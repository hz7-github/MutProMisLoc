"""Download article PDFs from Step 01 XLSX/CSV metadata.

Run: python scripts/02_download_papers.py --limit 3
Resume: repeat the same command without --limit.
Retry unsuccessful records: add --retry-failed.
Sources: Europe PMC discovery, official PMC cloud data, public publisher links.
The expanded preset is the default; --preset original selects the alternative paths.
Outputs: data/02_download_papers/<preset>/ without timestamp or batch folders.
"""
from __future__ import annotations

import argparse
import contextlib
import csv
import hashlib
import html
import json
import os
import re
import sqlite3
import signal
import threading
import time
import unicodedata
import uuid
import xml.etree.ElementTree as ET
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
from datetime import datetime, timezone
from difflib import SequenceMatcher
from pathlib import Path
from urllib.parse import parse_qs, quote, urlencode, urljoin, urlparse

import fitz
import requests
from bs4 import BeautifulSoup

ROOT = Path(__file__).resolve().parents[1]
INPUT_DIR = ROOT / 'data/01_get_pubmed_metadata'
DEFAULT_OUTPUT = ROOT / 'data/02_download_papers'
S3 = 'https://pmc-oa-opendata.s3.amazonaws.com/'
EPMC = 'https://www.ebi.ac.uk/europepmc/webservices/rest/search'
LOCAL = threading.local()
GATE = threading.Lock()
HOST_LOCKS = {}
HOST_LAST = {}
HOST_PAUSE = {}
STOP = threading.Event()


class DownloadCancelled(BaseException):
    """Do not let source fallback handlers swallow cancellation."""


def check_stop():
    if STOP.is_set():
        raise DownloadCancelled()


@contextlib.contextmanager
def interruptible_lock(lock):
    while not lock.acquire(timeout=0.2):
        check_stop()
    try:
        check_stop()
        yield
    finally:
        lock.release()


def request_stop(signum, frame):
    STOP.set()
FIELDS = ['PMID', 'PMCID', 'DOI', 'Title', 'status', 'source', 'source_url',
          'pdf_path', 'bytes', 'pages', 'sha256', 'identity_check', 'license',
          'article_version', 'is_manuscript', 'publication_alert', 'error', 'updated_at']


def now():
    return datetime.now(timezone.utc).isoformat()


def norm(text):
    text = html.unescape(html.unescape(str(text or '')))
    text = re.sub(r'<[^>]*>', '', text)
    return ''.join(c for c in unicodedata.normalize('NFKD', text).lower() if c.isalnum())


def title_match(a, b):
    return bool(norm(a) and norm(b)) and SequenceMatcher(None, norm(a), norm(b)).ratio() >= .88


def session():
    if not hasattr(LOCAL, 'session'):
        LOCAL.session = requests.Session()
        LOCAL.session.headers['User-Agent'] = 'MutProMisLoc-Downloader/1.0 (research fulltext retrieval)'
    return LOCAL.session


def fetch(url, attempts, *, max_bytes=100_000_000):
    """Bounded direct request, one request per publisher host, no challenge bypass."""
    check_stop()
    host = urlparse(url).netloc.lower()
    if urlparse(url).scheme not in ('http', 'https'):
        raise ValueError('Unsupported URL scheme')
    with GATE:
        lock = HOST_LOCKS.setdefault(host, threading.Lock())
    event = {'url': url, 'time': now()}
    attempts.append(event)
    # The public data bucket supports concurrent file retrieval. Serialize only
    # request starts there; holding this lock for a whole PDF would serialize all workers.
    cloud = host.endswith('s3.amazonaws.com')
    with interruptible_lock(lock):
        if HOST_PAUSE.get(host, 0) > time.monotonic():
            event['error'] = 'HOST_RETRY_AFTER'
            raise RuntimeError('HOST_RETRY_AFTER')
        STOP.wait(max(0, HOST_LAST.get(host, 0) + (.15 if host.endswith('s3.amazonaws.com') else 1) - time.monotonic()))
        check_stop()
        HOST_LAST[host] = time.monotonic()
    with contextlib.nullcontext() if cloud else interruptible_lock(lock):
        try:
            check_stop()
            with session().get(url, timeout=(8, 20), stream=True) as response:
                event.update(http_status=response.status_code, final_url=response.url)
                if response.status_code == 429:
                    delay = response.headers.get('Retry-After', '60')
                    HOST_PAUSE[host] = time.monotonic() + (int(delay) if delay.isdigit() else 60)
                response.raise_for_status()
                chunks = []
                size = 0
                started = time.monotonic()
                for chunk in response.iter_content(131072):
                    check_stop()
                    size += len(chunk)
                    pdf_stream = (chunks[0] if chunks else chunk).startswith(b'%PDF-')
                    if size > max_bytes or time.monotonic() - started > (600 if pdf_stream else 90):
                        raise ValueError('Response size/time limit exceeded')
                    chunks.append(chunk)
                body = b''.join(chunks)
                event['bytes'] = len(body)
                return body, response.url, response.headers.get('Content-Type', '')
        except Exception as exc:
            event['error'] = type(exc).__name__ + ': ' + str(exc)[:350]
            raise


def pdf_check(data, row=None, trusted_identity=False):
    if not data.startswith(b'%PDF-'):
        prefix = data[:4096].lstrip().lower()
        if any(tag in prefix for tag in (b'<!doctype html', b'<html', b'<head', b'<body')):
            raise ValueError('PDF_HTML_RESPONSE: HTML returned instead of a PDF; page purpose not verified')
        raise ValueError('PDF_MISSING_HEADER: response does not start with a PDF signature')
    if b'%%EOF' not in data[-8192:]:
        raise ValueError('PDF_MISSING_EOF: PDF end marker missing; response may be truncated')
    try:
        doc = fitz.open(stream=data, filetype='pdf')
    except Exception as exc:
        raise ValueError('PDF_PARSE_ERROR: unable to open PDF') from exc
    with doc:
        if doc.is_encrypted:
            raise ValueError('PDF_ENCRYPTED: PDF requires a password')
        if doc.page_count < 1:
            raise ValueError('PDF_EMPTY: PDF contains no pages')
        pages = doc.page_count
        # First pages only: a DOI in the references must not validate a different paper.
        try:
            text = '\n'.join(doc[i].get_text() for i in range(min(2, pages)))
        except Exception as exc:
            raise ValueError('PDF_PARSE_ERROR: unable to read PDF pages') from exc
    if row and not trusted_identity:
        front = norm(text[:18000])
        title = norm(row['Title'])
        aliases = [row.get(k, '') for k in ('DOI', 'DOI_before_independent_audit', 'DOI_canonical_alias')]
        has_doi = any(norm(d) in front for d in aliases if d.strip())
        words = set(re.findall(r'[a-zA-Z]{4,}', row['Title'].lower()))
        text_words = set(re.findall(r'[a-zA-Z]{4,}', text[:10000].lower()))
        title_coverage = len(words & text_words) / max(len(words), 1)
        if not ((title and title in front) or (has_doi and title_coverage >= .5)):
            raise ValueError('PDF_IDENTITY_UNCONFIRMED: title/DOI not found on first pages')
    return pages


def s3_url(value):
    if not value.startswith('s3://pmc-oa-opendata/'):
        raise ValueError('Unexpected cloud bucket')
    return S3 + value.split('s3://pmc-oa-opendata/', 1)[1].split('?', 1)[0]


def cloud_metadata(pmcid, out, attempts):
    # Exact PMCID + dot prevents matching e.g. PMC1234 when requesting PMC123.
    url = S3 + '?list-type=2&prefix=' + quote('metadata/' + pmcid + '.') + '&max-keys=1000'
    data, _, _ = fetch(url, attempts, max_bytes=2_000_000)
    root = ET.fromstring(data)
    keys = [e.text for e in root.findall('.//{*}Key') if e.text and e.text.endswith('.json')]
    if root.findtext('{*}IsTruncated') == 'true':
        raise ValueError('Unexpected metadata listing truncation')
    results = []
    for key in keys:
        content, _, _ = fetch(S3 + quote(key, safe='/'), attempts, max_bytes=1_000_000)
        results.append(json.loads(content))
    return results


def save_pdf(data, row, result, out, source, url, identity, md5=''):
    if md5 and hashlib.md5(data).hexdigest() != md5:
        raise ValueError('PDF_CHECKSUM_MISMATCH: PMC PDF MD5 differs from official metadata')
    pages = pdf_check(data, row, trusted_identity=identity)
    path = out / (row['PMID'] + '.pdf')
    temp = path.with_suffix('.pdf.part')
    temp.write_bytes(data)
    os.replace(temp, path)
    result.update(status='SUCCESS', source=source, source_url=url, pdf_path=str(path),
                  bytes=len(data), pages=pages, sha256=hashlib.sha256(data).hexdigest(),
                  identity_check='PMC_METADATA_PMID_AND_PMCID' if identity else 'PDF_FIRST_PAGES_TITLE_DOI', error='')


def publisher_pdf_links(body, base, row):
    soup = BeautifulSoup(body, 'html.parser')
    title_tag = soup.find('meta', attrs={'name':re.compile(r'^citation_title$', re.I)})
    doi_tag = soup.find('meta', attrs={'name':re.compile(r'^citation_doi$', re.I)})
    page_title = title_tag.get('content', '') if title_tag else ''
    page_doi = doi_tag.get('content', '') if doi_tag else ''
    expected = {norm(row.get(k, '')) for k in ('DOI','DOI_before_independent_audit','DOI_canonical_alias')} - {''}
    # Do not follow PDFs on search results, unrelated redirect pages, or generic portals.
    if not (title_match(page_title, row['Title']) or (page_doi and norm(page_doi) in expected)
            or (norm(row['Title']) and norm(row['Title']) in norm(soup.get_text(' ', strip=True)[:16000]))):
        return []
    links = []
    for meta in soup.find_all('meta', attrs={'name':re.compile(r'^citation_pdf_url$', re.I)}):
        if meta.get('content'):
            links.append(urljoin(base, meta['content']))
    for a in soup.find_all('a', href=True):
        href, text = a['href'], a.get_text(' ', strip=True).lower()
        if ('pdf' in href.lower() or 'pdf' in text) and not re.search(r'supp|appendix|figure|citation|reference', href+' '+text, re.I):
            links.append(urljoin(base, href))
    return list(dict.fromkeys(u for u in links if urlparse(u).scheme in ('http','https')))[:3]


def process(row, epmc, out, existing):
    attempts = []
    result = {'PMID':row['PMID'], 'PMCID':epmc.get('pmcid') or row.get('PMCID', ''),
        'DOI':row.get('DOI',''), 'Title':row['Title'], 'status':'NOT_FOUND', 'source':'', 'source_url':'',
        'pdf_path':'', 'bytes':0, 'pages':0, 'sha256':'', 'identity_check':'', 'license':'',
        'article_version':'', 'is_manuscript':'', 'publication_alert':row.get('publication_alert',''),
        'error':'', 'updated_at':now(), 'attempts':attempts}
    for path in existing.get(row['PMID'], []):
        try:
            pages = pdf_check(path.read_bytes())
            result.update(status='SUCCESS_EXISTING', source='LOCAL', pdf_path=str(path), pages=pages,
                          bytes=path.stat().st_size, identity_check='EXISTING_FILENAME_PMID_PDF_PARSE')
            return result
        except Exception:
            continue
    pmcid = result['PMCID']
    if re.fullmatch(r'PMC\d+', pmcid):
        try:
            versions = cloud_metadata(pmcid, out, attempts)
            # Prefer the published article, not blindly the highest PMC version.
            versions.sort(key=lambda x:(bool(x.get('is_manuscript')), -int(x.get('version', 0))))
            for meta in versions:
                if str(meta.get('pmid')) != row['PMID'] or meta.get('pmcid') != pmcid:
                    attempts.append({'error':'PMC_METADATA_IDENTITY_MISMATCH', 'metadata':meta})
                    continue
                if meta.get('is_retracted'):
                    result['publication_alert'] = 'PMC_RETRACTED'
                if not meta.get('pdf_url'):
                    continue
                try:
                    url = s3_url(meta['pdf_url'])
                    data, final, _ = fetch(url, attempts)
                    md5 = parse_qs(urlparse(meta['pdf_url']).query).get('md5',[''])[0]
                    save_pdf(data, row, result, out, 'PMC_CLOUD', final, True, md5)
                    result.update(license=meta.get('license_code',''), article_version=meta.get('version',''),
                                  is_manuscript=meta.get('is_manuscript',''))
                    return result
                except Exception as exc:
                    attempts.append({'url':meta['pdf_url'], 'error':str(exc)[:400]})
        except Exception as exc:
            attempts.append({'stage':'PMC_CLOUD', 'error':str(exc)[:400]})
    seen = set()
    for url in publisher_candidates(row, epmc):
        if url in seen:
            continue
        seen.add(url)
        try:
            data, final, kind = fetch(url, attempts)
            if data.startswith(b'%PDF-'):
                save_pdf(data, row, result, out, 'PUBLISHER', final, False)
                return result
            for pdf_url in publisher_pdf_links(data, final, row):
                if pdf_url in seen:
                    continue
                seen.add(pdf_url)
                try:
                    pdf, pdf_final, _ = fetch(pdf_url, attempts)
                    save_pdf(pdf, row, result, out, 'PUBLISHER', pdf_final, False)
                    return result
                except Exception as exc:
                    attempts.append({'url':pdf_url, 'error':str(exc)[:400]})
        except Exception as exc:
            attempts.append({'url':url,'error':str(exc)[:400]})
    result['status'], result['error'] = classify_failure(attempts)
    return result


def classify_failure(attempts):
    """Summarize mixed failures while preserving every reason in the checkpoint."""
    errors = [str(item['error']) for item in attempts if item.get('error')]
    invalid = ('PDF_HTML_RESPONSE', 'PDF_MISSING_HEADER', 'PDF_MISSING_EOF',
               'PDF_ENCRYPTED', 'PDF_EMPTY', 'PDF_PARSE_ERROR', 'PDF_CHECKSUM_MISMATCH',
               'Not a complete PDF', 'Encrypted or empty PDF', 'PMC PDF MD5 mismatch')
    if any('IDENTITY' in error for error in errors):
        status = 'IDENTITY_REVIEW'
    elif any(any(code in error for code in invalid) for error in errors):
        status = 'INVALID_PDF'
    elif any(item.get('http_status') in (401, 403, 429) for item in attempts):
        status = 'ACCESS_BLOCKED'
    elif errors or any(int(item.get('http_status') or 0) >= 400 for item in attempts):
        status = 'REQUEST_FAILED'
    else:
        status = 'NO_PUBLIC_PDF_LINK'
    return status, ' | '.join(dict.fromkeys(errors))[-2000:] or 'No public PDF link found among attempted sources'


def publisher_candidates(row, epmc):
    urls = []
    for link in epmc.get('fullTextUrlList', {}).get('fullTextUrl', []):
        if link.get('availabilityCode') in ('OA', 'F'):
            url = link.get('url', '')
            if urlparse(url).scheme in ('http', 'https') and not any(
                    host in urlparse(url).netloc for host in ('europepmc.org', 'ncbi.nlm.nih.gov')):
                urls.append(url)
    if row.get('DOI'):
        urls.append('https://doi.org/' + quote(row['DOI'], safe='/():'))
    return list(dict.fromkeys(urls))[:6]


def download_one(row, out):
    check_stop()
    events = []
    record = {}
    try:
        query = urlencode({'query': 'EXT_ID:' + row['PMID'] + ' AND SRC:MED',
                           'format': 'json', 'resultType': 'core', 'pageSize': 2})
        body, _, _ = fetch(EPMC + '?' + query, events, max_bytes=2_000_000)
        data = json.loads(body)
        for item in data.get('resultList', {}).get('result', []):
            if str(item.get('id')) == row['PMID'] and item.get('source') == 'MED':
                record = item
                break
    except Exception as exc:
        events.append({'stage': 'EUROPE_PMC', 'error': str(exc)[:400]})
    result = process(row, record, out, {})
    result['attempts'] = events + result.get('attempts', [])
    if not result['status'].startswith('SUCCESS'):
        result['status'], result['error'] = classify_failure(result['attempts'])
    return result


def load_input(path):
    if path.suffix.lower() == '.csv':
        with path.open(encoding='utf-8-sig', newline='') as stream:
            reader = csv.DictReader(stream)
            headers = reader.fieldnames or []
            records = list(reader)
    elif path.suffix.lower() == '.xlsx':
        import openpyxl
        book = openpyxl.load_workbook(path, read_only=True, data_only=True)
        try:
            sheet = book['Papers'] if 'Papers' in book.sheetnames else book.active
            values = sheet.iter_rows(values_only=True)
            headers = list(next(values, []))
            records = [dict(zip(headers, values)) for values in values]
        finally:
            book.close()
    else:
        raise ValueError('Input must be an XLSX or CSV file')
    if not {'PMID', 'Title'} <= set(headers):
        raise ValueError('Input requires PMID and Title columns')
    rows = []
    seen = set()
    for record in records:
        if not any(v is not None and str(v).strip() for v in record.values()):
            continue
        row = {k: str(record.get(k) or '').strip() for k in ('PMID', 'PMCID', 'DOI', 'Title')}
        if not re.fullmatch(r'[1-9][0-9]*', row['PMID']) or not row['Title']:
            raise ValueError('Every row requires a numeric PMID and a nonempty Title')
        if row['PMID'] in seen:
            raise ValueError('Duplicate PMID: ' + row['PMID'])
        if row['PMCID'] and not re.fullmatch(r'PMC[0-9]+', row['PMCID']):
            raise ValueError('Invalid PMCID for PMID ' + row['PMID'])
        row['DOI'] = re.sub(r'^https?://(?:dx\.)?doi\.org/|^doi:\s*', '', row['DOI'], flags=re.I)
        seen.add(row['PMID'])
        rows.append(row)
    if not rows:
        raise ValueError('Input contains no papers')
    return rows


def fingerprint(row):
    return hashlib.sha256(json.dumps(row, sort_keys=True).encode()).hexdigest()


def export_status(out, rows, results):
    target = out / 'download_status.csv'
    temp = out / ('download_status.' + uuid.uuid4().hex + '.tmp')
    try:
        with temp.open('w', encoding='utf-8-sig', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=FIELDS, extrasaction='ignore')
            writer.writeheader()
            for row in rows:
                writer.writerow(results.get(row['PMID'], dict(row, status='PENDING')))
        os.replace(temp, target)
    except PermissionError:
        print('Close download_status.csv to refresh it; progress is safe in the checkpoint.', flush=True)
    finally:
        if temp.exists():
            temp.unlink()


@contextlib.contextmanager
def output_lock(out):
    with (out / '.download.lock').open('a+b') as handle:
        handle.seek(0)
        handle.write(b'0')
        handle.flush()
        handle.seek(0)
        if os.name == 'nt':
            import msvcrt
            acquire = lambda: msvcrt.locking(handle.fileno(), msvcrt.LK_NBLCK, 1)
            release = lambda: msvcrt.locking(handle.fileno(), msvcrt.LK_UNLCK, 1)
        else:
            import fcntl
            acquire = lambda: fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
            release = lambda: fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
        try:
            acquire()
        except OSError as exc:
            raise RuntimeError('Another downloader is using this output directory') from exc
        try:
            yield
        finally:
            handle.seek(0)
            release()


def run(args, rows):
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    with output_lock(out), contextlib.closing(sqlite3.connect(out / '.download_state.sqlite3')) as db:
        db.execute('CREATE TABLE IF NOT EXISTS records (pmid TEXT PRIMARY KEY, fingerprint TEXT, result TEXT)')
        db.commit()
        results = {}
        todo = []
        for row in rows:
            if STOP.is_set():
                return
            saved = db.execute('SELECT fingerprint, result FROM records WHERE pmid=?', (row['PMID'],)).fetchone()
            old = json.loads(saved[1]) if saved and saved[0] == fingerprint(row) else {}
            if old.get('status') == 'NO_PUBLIC_PDF':
                # Reclassify legacy evidence without inventing a specific response type.
                events = old.get('attempts', [])
                if not events and old.get('error') != 'No verified public PDF link found':
                    events = [{'error': old.get('error', '')}]
                old['status'], old['error'] = classify_failure(events)
                db.execute('UPDATE records SET result=? WHERE pmid=?', (json.dumps(old), row['PMID']))
                db.commit()
            path = out / (row['PMID'] + '.pdf')
            if path.exists():
                try:
                    content = path.read_bytes()
                    trusted = bool(old.get('sha256') and old['sha256'] == hashlib.sha256(content).hexdigest())
                    pages = pdf_check(content, row, trusted_identity=trusted)
                    result = dict(old or row, status='SUCCESS_EXISTING', pdf_path=str(path),
                                  pages=pages, bytes=len(content), sha256=hashlib.sha256(content).hexdigest(),
                                  updated_at=now(), error='')
                    if not trusted:
                        result.update(source='LOCAL', identity_check='PDF_FIRST_PAGES_TITLE_DOI')
                    results[row['PMID']] = result
                    db.execute('INSERT OR REPLACE INTO records VALUES (?, ?, ?)',
                               (row['PMID'], fingerprint(row), json.dumps(result)))
                    db.commit()
                    continue
                except (ValueError, RuntimeError, OSError):
                    pass
            if old and not old.get('status', '').startswith('SUCCESS') and not args.retry_failed:
                results[row['PMID']] = old
                continue
            todo.append(row)
        if args.limit is not None:
            todo = todo[:args.limit]
        print('Input:', len(rows), 'queued:', len(todo), 'workers:', args.workers, flush=True)
        export_status(out, rows, results)
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            futures = {pool.submit(download_one, row, out): row for row in todo}
            pending = set(futures)
            number = 0
            stopping = False
            while pending:
                if STOP.is_set() and not stopping:
                    stopping = True
                    print('Stopping: cancelling queued tasks; waiting for active requests to return. Completed results are saved.', flush=True)
                    for future in pending:
                        future.cancel()
                # Remove cancelled futures explicitly: wait() may not notify for
                # futures cancelled before a worker has picked them up.
                pending = {f for f in pending if not f.cancelled()}
                if not pending:
                    break
                done, pending = wait(pending, timeout=0.2, return_when=FIRST_COMPLETED)
                for future in done:
                    row = futures[future]
                    try:
                        result = future.result()
                    except DownloadCancelled:
                        continue
                    except Exception as exc:
                        result = dict(row, status='INTERNAL_ERROR', error=type(exc).__name__ + ': ' + str(exc))
                    result['updated_at'] = now()
                    db.execute('INSERT OR REPLACE INTO records VALUES (?, ?, ?)',
                               (row['PMID'], fingerprint(row), json.dumps(result)))
                    db.commit()
                    results[row['PMID']] = result
                    export_status(out, rows, results)
                    number += 1
                    print(f"[{number}/{len(todo)}] PMID={row['PMID']} {result['status']}", flush=True)
        print('Status:', dict(Counter(results.get(r['PMID'], {}).get('status', 'PENDING') for r in rows)))
        print('Output:', out)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--preset', choices=('original', 'expanded'), default='expanded',
                        help='Input/output preset (default: expanded)')
    parser.add_argument('--input', type=Path, help='Override papers_<preset>.xlsx with an XLSX or CSV file')
    parser.add_argument('--output-dir', '--output', dest='output_dir', type=Path,
                        help='Override the exact output directory; no preset suffix is appended')
    parser.add_argument('--workers', type=int, choices=range(1, 9), default=3)
    parser.add_argument('--limit', type=int, help='Maximum new download attempts; 0 initializes status only')
    parser.add_argument('--retry-failed', action='store_true', help='Retry previous unsuccessful records')
    parser.add_argument('--plan-only', action='store_true', help='Validate input and show paths without network or writes')
    args = parser.parse_args(argv)
    if args.limit is not None and args.limit < 0:
        parser.error('--limit must be non-negative')
    args.input = args.input or INPUT_DIR / f'papers_{args.preset}.xlsx'
    args.output_dir = args.output_dir or DEFAULT_OUTPUT / args.preset
    return args


def main():
    args = parse_args()
    rows = load_input(args.input)
    if args.plan_only:
        print('Preset:', args.preset)
        print('Papers:', len(rows))
        print('Input:', args.input.resolve())
        print('Output:', args.output_dir.resolve())
        print('Sources: Europe PMC discovery, PMC cloud, public publisher links')
        return
    STOP.clear()
    previous = signal.signal(signal.SIGINT, request_stop)
    try:
        run(args, rows)
    finally:
        signal.signal(signal.SIGINT, previous)
    if STOP.is_set():
        print('Stopped. Repeat the command to resume unfinished papers.', flush=True)
        raise SystemExit(130)


if __name__ == '__main__':
    main()
