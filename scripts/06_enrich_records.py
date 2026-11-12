"""Enrich verified records and export rows with complete target fields.

Requires Accession_new, Reviewed and Wild Sequence_new. Database location is
required only where the existing no-change applicability rule permits it.
Offline by default; incomplete records remain in work/ and audit/.
ClinVar optionally fills clinical classifications; --online allows new queries.
Pipeline placement: Step 05 supplies source-verified records and normalized
variants. This step resolves protein identity/reference sequences, checks export
eligibility, then annotates clinical classifications before the final export.
The main workflow and Demo both use this entry point.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
import os
import re
import sys
import threading
import time
from collections import Counter
from concurrent.futures import Future, ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from http.client import HTTPException
from pathlib import Path
from types import SimpleNamespace
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts._internal import identity_recovery as recovery
from scripts._internal import legacy_enrichment as enrich
display = enrich.N.DISPLAY

VERSION = 'complete-records-enrichment-1.5.1'
DATA_SHEETS = ('Records', 'Completed Records')
LEGACY_SHEETS = ('A_candidates_26col', 'A_records_26col', 'Completed A Records')
COLUMNS = enrich.N.BASE_COLUMNS[:6] + [enrich.N.DETAIL] + enrich.N.BASE_COLUMNS[6:]


def load_records_input(path, sheet=None):
    """Require export provenance or a recognized data worksheet."""
    sources = {}
    if path.suffix.lower() == '.json':
        data = enrich.read_json(path)
        entries = data.get('records') if isinstance(data, dict) else data
        if not isinstance(entries, list):
            raise ValueError('Expected a records array.')
        if isinstance(data, dict) and data.get('tier') == 'A':
            proof = data.get('eligibility', {})
            rows = entries
            if any(proof.get(r.get('id'), {}).get('stage_status') != 'STAGE_A' for r in rows):
                raise ValueError('Every exported record needs an eligibility entry from source review.')
            sources = data.get('source_members', {})
        elif all(isinstance(r, dict) and 'stage_status' in r for r in entries):
            rows = [r['row'] for r in entries if r['stage_status'] == 'STAGE_A'
                    and not r.get('validation_errors') and not r.get('source_context_flags')
                    and r.get('canonical_event_record_id', r['id']) == r['id']]
        else:
            raise ValueError('Use Step 05 records.json, a stage-status ledger, or an exported data worksheet.')
    elif path.suffix.lower() == '.xlsx':
        from openpyxl import load_workbook
        book = load_workbook(path, read_only=True, data_only=True)
        try:
            accepted = DATA_SHEETS + LEGACY_SHEETS
            selected = sheet or next((name for name in accepted if name in book.sheetnames), None)
            if selected not in accepted or selected not in book.sheetnames:
                raise ValueError('Input workbook must contain an exported data worksheet: ' + ', '.join(DATA_SHEETS))
        finally:
            book.close()
        _, rows = enrich.load_rows(path, selected)
    else:
        raise ValueError('Use exported JSON or XLSX; an unlabelled CSV lacks source-review provenance.')
    if path.suffix.lower() == '.xlsx' and path.with_suffix('.json').is_file():
        adjacent = enrich.read_json(path.with_suffix('.json'))
        if isinstance(adjacent, dict) and adjacent.get('display_originals'):
            companion = {r['id']: r for r in adjacent.get('records', [])}
            if set(companion) != {r['id'] for r in rows} or any(
                    enrich.text(row.get(h)) != enrich.text(companion[row['id']].get(h))
                    for row in rows for h in COLUMNS):
                raise ValueError('Adjacent JSON does not match the displayed workbook; use the correct records.json.')
            data = adjacent
            sources = adjacent.get('source_members', sources)
        else:
            data = None
    elif path.suffix.lower() != '.json':
        data = None
    if isinstance(data, dict):
        originals = data.get('display_originals', {})
        rows = [display.restore_originals(row, originals.get(row.get('id'), {})) for row in rows]
    seen = set()
    for row in rows:
        if not isinstance(row, dict) or not set(COLUMNS) <= set(row):
            raise ValueError('A record is missing required 26-column headers.')
        row['id'], row['PMID'] = enrich.text(row['id']), enrich.text(row['PMID'])
        if not row['id'] or row['id'] in seen or not re.fullmatch(r'[1-9]\d*', row['PMID']):
            raise ValueError('Empty/duplicate ID or invalid PMID: ' + row['id'])
        if enrich.text(row.get('Wild Sequence_new')).startswith('[Full sequence:'):
            raise ValueError('Use JSON for long sequences; the workbook cell contains a reference marker.')
        seen.add(row['id'])
    return [{h: r.get(h) for h in COLUMNS} for r in rows], sources


def completeness(row):
    missing, errors = [], []
    def absent(value):
        return enrich.missing(value) or enrich.text(value).casefold() in {'unresolved', 'undetermined', 'not_applicable', 'not applicable', 'pending'}
    for field in ['Accession_new', 'Reviewed', 'Wild Sequence_new']:
        if absent(row.get(field)):
            missing.append(field)
    location_required = enrich.database_location_applicable(row)
    if location_required and absent(row.get('Subcellular_Location_Values')):
        missing.append('Subcellular_Location_Values')
    if 'Accession_new' not in missing and not enrich.ACCESSION.fullmatch(enrich.text(row['Accession_new'])):
        errors.append('INVALID_ACCESSION')
    if 'Reviewed' not in missing and enrich.text(row['Reviewed']).casefold() not in {'reviewed', 'unreviewed'}:
        errors.append('INVALID_REVIEWED_STATUS')
    if 'Wild Sequence_new' not in missing:
        sequence = enrich.text(row['Wild Sequence_new'])
        if not re.fullmatch('[ACDEFGHIKLMNPQRSTVWYBXZJUO]+', sequence):
            errors.append('INVALID_REFERENCE_SEQUENCE')
        elif enrich.sequence_check(row['Mutation Information'], sequence) in {'REFERENCE_RESIDUE_MISMATCH', 'POSITION_OUT_OF_RANGE'}:
            errors.append('REFERENCE_MISMATCH')
    return {'id': row['id'], 'PMID': row['PMID'], 'missing_fields': missing,
            'validation_errors': errors, 'database_location_required': location_required,
            'exported': not missing and not errors}


def prepare_sources(rows, embedded, args):
    source_root = (getattr(args, 'source_root', None) or ROOT).resolve()
    sources = enrich.HistoricalSources(args.source_audit, args.original_records, args.dedup_map)
    sources.load(rows)
    if args.identity_input:
        identities = enrich.read_json(args.identity_input)
        if len({r['id'] for r in identities}) != len(identities):
            raise ValueError('Duplicate identity-input IDs.')
        identities = {r['id']: r for r in identities}
    else:
        identities = {}
    for row in rows:
        members = embedded.get(row['id'])
        if not members and row['id'] in identities:
            item = identities[row['id']]
            if str(item['PMID']) != row['PMID']:
                raise ValueError('Identity input PMID differs: ' + row['id'])
            members = [{'id': row['id'], 'record': item['original_record'],
                        'source_path': item.get('source_path'), 'source_sha256': item['source_sha256']}]
            embedded[row['id']] = members
        if not members:
            continue
        if not isinstance(members, list) or not members or len({r['id'] for r in members}) != len(members):
            raise ValueError('Invalid source-member mapping: ' + row['id'])
        sources.groups[row['id']] = [r['id'] for r in members]
        for member in members:
            raw = member['record']
            if raw.get('pmid') and str(raw['pmid']) != row['PMID']:
                raise ValueError('Source member crosses papers: ' + row['id'])
            if member.get('source_path'):
                path = (source_root / member['source_path']).resolve()
                if not path.is_relative_to(source_root) or enrich.digest(path) != member.get('source_sha256'):
                    raise ValueError('Source member file changed: ' + row['id'])
                sources.hashes[str(path)] = member['source_sha256']
            sources.records[member['id']] = raw
            sources.locators[member['id']] = 'A source member: ' + member['id']
    return sources


def build_recovery_workset(rows, audit, sources, embedded, input_path):
    by_id = {r['id']: r for r in audit}; targets = []
    for row in rows:
        a = by_id[row['id']]
        if not enrich.missing(row.get('Wild Sequence_new')):
            continue
        if 'SEQUENCE_WITHHELD_ISOFORM_CONTEXT' in a['notes']:
            cohort = 'isoform'
        elif a['uniprot_status'] != 'UNIQUE':
            cohort = 'identity'
        else:
            continue
        item = {'row': row, 'prior_audit': a, 'cohort': cohort,
                'originals': [{'id': mid, 'record': sources.records.get(mid, {})} for mid in sources.members(row['id'])]}
        paths = {(m.get('source_path'), m.get('source_sha256')) for m in embedded.get(row['id'], []) if m.get('source_path')}
        if len(paths) > 1:
            raise ValueError('Multiple source bindings for one A record: ' + row['id'])
        if paths:
            item['source_path'], item['source_sha256'] = next(iter(paths))
        targets.append(item)
    return {'input_sha256': enrich.digest(input_path), 'records': targets}


CLINICAL_FIELD = 'Variant Clinical Classification'
NUCLEOTIDE_HGVS = re.compile(
    r'(?<![A-Za-z0-9_])N[CMGRWT]_\d+\.\d+:[cgmn]\.[A-Za-z0-9_+*?>=()\-]+')


def ncbi_api_key(env_file=None):
    """Return a credential without logging it; an explicit file takes priority."""
    def from_file(path):
        if not path.is_file():
            return None
        values = []
        for line in path.read_text(encoding='utf-8-sig').splitlines():
            match = re.match(r'^\s*(?:export\s+)?NCBI_API_KEY\s*=\s*(.*?)\s*$', line)
            if match:
                value = match[1]
                if value[:1] in {'"', "'"} and value[-1:] == value[:1]:
                    value = value[1:-1]
                else:
                    value = re.split(r'\s+#', value, maxsplit=1)[0].strip()
                values.append(value)
        if len(values) > 1:
            raise ValueError('Duplicate NCBI_API_KEY entries in an environment file.')
        return values[0] if values else None
    if env_file is not None:
        if not env_file.is_file():
            raise ValueError('--env-file does not exist.')
        value = from_file(env_file)
        if value is not None:
            return value
    if os.environ.get('NCBI_API_KEY', '').strip():
        return os.environ['NCBI_API_KEY'].strip()
    for path in (ROOT / '.env', Path.home() / '.env'):
        value = from_file(path)
        if value:
            return value
    return ''


class ClinVarClient:
    """POST-only NCBI requests, credential-free caching and bounded retries.

    API: https://www.ncbi.nlm.nih.gov/clinvar/docs/maintenance_use/
    """
    BASE = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'

    def __init__(self, args, remaining):
        self.key = ncbi_api_key(args.env_file)
        self.cache = args.clinvar_cache
        self.online, self.limit = args.online, remaining
        self.timeout, self.retries = args.timeout, args.retries
        self.interval = max(args.interval, 0.11 if self.key else 0.35)
        self.max_age = args.cache_max_age_days * 86400
        self.calls = self.cache_hits = 0
        self.last_request = 0.0
        self.auth_failed = False
        self.memo = {}
        self.pending = {}
        self.pending_lock = threading.Lock()
        self.request_lock = threading.Lock()

    def request(self, endpoint, **params):
        token = json.dumps([endpoint, params], sort_keys=True)
        with self.pending_lock:
            future = self.pending.get(token)
            owner = future is None
            if owner:
                future = Future()
                self.pending[token] = future
            else:
                self.cache_hits += 1
        if not owner:
            return future.result()
        try:
            result = self._request(endpoint, **params)
        except Exception as exc:
            future.set_exception(exc)
            raise
        future.set_result(result)
        return result

    def _request(self, endpoint, **params):
        public = {'db': 'clinvar', **params}
        token = hashlib.sha256(json.dumps([endpoint, public], sort_keys=True).encode()).hexdigest()
        if token in self.memo:
            self.cache_hits += 1
            return self.memo[token]
        path = self.cache / (token + '.json')
        if path.is_file():
            try:
                saved = enrich.read_json(path)
                payload = saved['payload']
                age = time.time() - datetime.fromisoformat(saved['retrieved_utc']).timestamp()
                valid = (saved.get('endpoint') == endpoint and saved.get('params') == public
                         and hashlib.sha256(payload.encode()).hexdigest() == saved.get('sha256')
                         and 0 <= age <= self.max_age)
                if valid:
                    self.cache_hits += 1
                    self.memo[token] = (payload, 'OK')
                    return payload, 'OK'
            except (ValueError, KeyError, TypeError, AttributeError, OSError):
                pass
        if not self.online:
            return None, 'OFFLINE_CACHE_MISS_OR_STALE'
        if self.auth_failed:
            return None, 'AUTHENTICATION_ERROR'
        private = {**public, 'tool': 'MutProMisLoc_Step06'}
        if self.key:
            private['api_key'] = self.key
        for attempt in range(self.retries + 1):
            with self.request_lock:
                if self.auth_failed:
                    return None, 'AUTHENTICATION_ERROR'
                if self.limit is not None and self.calls >= self.limit:
                    return None, 'REQUEST_LIMIT'
                time.sleep(max(0, self.interval - (time.monotonic() - self.last_request)))
                self.last_request = time.monotonic()
                self.calls += 1
            req = Request(self.BASE + endpoint, data=urlencode(private).encode(),
                          headers={'User-Agent': 'MutProMisLoc-ClinVar/1.0'})
            try:
                with urlopen(req, timeout=self.timeout) as response:
                    raw = response.read()
            except HTTPError as exc:
                body = exc.read().lower()
                if b'api key' in body and b'invalid' in body:
                    self.auth_failed = True
                    return None, 'AUTHENTICATION_ERROR'
                if (exc.code == 429 or exc.code >= 500) and attempt < self.retries:
                    time.sleep(min(4, 2 ** attempt))
                    continue
                return None, 'HTTP_' + str(exc.code)
            except (URLError, TimeoutError, OSError, HTTPException):
                if attempt < self.retries:
                    time.sleep(min(4, 2 ** attempt))
                    continue
                return None, 'CONNECTION_ERROR'
            if self.key and self.key.encode() in raw:
                return None, 'RESPONSE_SUPPRESSED_CREDENTIAL'
            try:
                payload = raw.decode('utf-8')
                if endpoint == 'esearch.fcgi':
                    obj = json.loads(payload)
                    search = obj.get('esearchresult', {})
                    errors = search.get('errorlist') or {}
                    warnings = search.get('warninglist') or {}
                    # ESearch returns errorlist with empty arrays on successful
                    # requests. An unmatched phrase at count=0 is a valid miss.
                    if ('error' in obj or 'ERROR' in obj or 'ERROR' in search
                            or 'count' not in search or 'idlist' not in search
                            or errors.get('fieldsnotfound') or warnings.get('phrasesignored')
                            or (str(search.get('count')) != '0' and
                                (any(errors.values()) or warnings.get('quotedphrasesnotfound')))):
                        return None, 'API_ERROR'
                elif params.get('retmode') == 'text' and params.get('rettype') == 'fasta':
                    if not payload.startswith('>'):
                        return None, 'INVALID_FASTA_RESPONSE'
                else:
                    xml = ET.fromstring(payload)
                    if any(n.tag.upper() == 'ERROR' for n in xml.iter()):
                        return None, 'API_ERROR'
            except (ValueError, TypeError, AttributeError, ET.ParseError):
                return None, 'INVALID_RESPONSE'
            enrich.atomic_json(path, {'endpoint': endpoint, 'params': public,
                'retrieved_utc': datetime.now(timezone.utc).isoformat(), 'payload': payload,
                'sha256': hashlib.sha256(payload.encode()).hexdigest()})
            self.memo[token] = (payload, 'OK')
            return payload, 'OK'
        return None, 'CONNECTION_ERROR'


def parse_clinvar_vcv(payload):
    """Use the aggregate germline label, preserving mixed/conflicting labels."""
    records = []
    for archive in ET.fromstring(payload).iter('VariationArchive'):
        allele = archive.find('./ClassifiedRecord/SimpleAllele')
        is_single = allele is not None
        if allele is None:
            allele = ET.Element('UnsupportedAllele')
        germline = archive.find('./ClassifiedRecord/Classifications/GermlineClassification')
        records.append({
            'variation_id': archive.get('VariationID'), 'accession': archive.get('Accession'),
            'version': archive.get('Version'),
            'record_type': 'SINGLE_ALLELE' if is_single else 'HAPLOTYPE_OR_GENOTYPE',
            'genes': sorted({n.get('Symbol') for n in allele.findall('./GeneList/Gene') if n.get('Symbol')}),
            'nucleotide_hgvs': sorted({n.text for n in allele.findall('./HGVSlist/HGVS/NucleotideExpression/Expression') if n.text}),
            'protein_hgvs': sorted({n.text for n in allele.findall('./HGVSlist/HGVS/ProteinExpression/Expression') if n.text}),
            'classification': germline.findtext('Description') if germline is not None else None,
            'review_status': germline.findtext('ReviewStatus') if germline is not None else None,
            'last_evaluated': germline.get('DateLastEvaluated') if germline is not None else None,
            'grch38_locations': [dict(n.attrib) for n in allele.findall('./Location/SequenceLocation')
                                if n.get('Assembly') == 'GRCh38'],
            'url': 'https://www.ncbi.nlm.nih.gov/clinvar/variation/' + str(archive.get('VariationID')) + '/',
        })
    return records


def protein_change(value):
    """Normalize predicted parentheses only, preserving positions and operations."""
    return enrich.text(value).split(':', 1)[-1].replace('(', '').replace(')', '')


def clinical_label_missing(value):
    return enrich.missing(value) or enrich.text(value).casefold() in {
        'not provided', 'not classified', 'no classification provided', 'not applicable'}


def clinvar_reference_match(row, candidate, protein, client):
    """Bind a protein-level hit to the table's existing reference protein."""
    accession = enrich.text(row.get('Accession_new'))
    expressions = candidate['protein_hgvs']
    expected = protein_change(protein)
    reference = enrich.text(row.get('Wild Sequence_new'))
    if reference and enrich.sequence_check(protein, reference) in {'REFERENCE_RESIDUE_MISMATCH', 'POSITION_OUT_OF_RANGE'}:
        return 'REFERENCE_RESIDUE_MISMATCH', None
    direct = [h for h in expressions if h.split(':', 1)[0] == accession]
    if direct:
        return ('UNIPROT_ACCESSION_MATCH', accession) if any(
            protein_change(h) == expected for h in direct) else ('REFERENCE_ACCESSION_CONFLICT', None)
    if not reference:
        return 'REFERENCE_UNAVAILABLE', None
    refs = sorted({h.split(':', 1)[0] for h in expressions
                   if protein_change(h) == expected and re.match(r'^(?:NP|XP)_\d+\.\d+:', h)})
    if not refs:
        return 'REFERENCE_ACCESSION_UNRESOLVED', None
    if len(refs) > 30:
        return 'REFERENCE_CANDIDATES_TRUNCATED', None
    payload, status = client.request('efetch.fcgi', db='protein', id=','.join(refs),
                                     rettype='fasta', retmode='text')
    if status != 'OK':
        return status, None
    seen = set()
    for block in payload.split('>'):
        if not block.strip():
            continue
        header, *lines = block.splitlines()
        identifiers = set(re.findall(r'(?:NP|XP)_\d+\.\d+', header))
        linked = identifiers & set(refs)
        seen.update(linked)
        sequence = ''.join(lines).replace(' ', '').strip().rstrip('*')
        if linked and sequence == reference:
            return 'REFSEQ_FULL_SEQUENCE_MATCH', sorted(linked)[0]
    return ('REFERENCE_SEQUENCE_MISMATCH' if seen == set(refs) else 'REFERENCE_RESPONSE_INCOMPLETE'), None


def clinvar_report(client, audit, changes):
    return {'enabled': client is not None, 'http_requests': client.calls if client else 0,
            'cache_hits': client.cache_hits if client else 0, 'filled_records': len(changes),
            'queried_records': sum(bool(r.get('lookup_performed')) for r in audit),
            'records_with_search_hits': sum(r.get('search_matches', 0) > 0 for r in audit),
            'status_counts': dict(Counter(r['status'] for r in audit)),
            'api_key_configured': bool(client and client.key)}


def enrich_clinvar(rows, original_rows, sources, args, remaining):
    """Match gene plus complete protein HGVS, verifying reference and consensus.

    A protein-level annotation is identified as such in the audit; it does not
    establish a nucleotide allele. Exact nucleotide matches remain supported.
    """
    if len(rows) > 1 and not getattr(args, '_clinvar_serial', False) and not args.skip_clinvar:
        client = ClinVarClient(args, remaining)
        worker_args = copy.copy(args)
        worker_args._clinvar_serial = True
        worker_args._clinvar_client = client
        workers = max(1, getattr(args, 'clinvar_workers', 4))
        originals = {r['id']: r for r in original_rows}
        audit, changes = [None] * len(rows), []
        def resolve(index, row):
            local = dict(row)
            try:
                items, updates, _ = enrich_clinvar([local], [originals[row['id']]], sources, worker_args, remaining)
            except (OSError, ValueError, TypeError, KeyError, IndexError, HTTPException) as exc:
                # A malformed/partial public response must not abort the entire
                # table. Keep a credential-free error type for a targeted retry.
                return index, dict(row), {'id': row['id'], 'PMID': row['PMID'],
                    'status': 'LOOKUP_PROCESSING_ERROR', 'lookup_performed': True,
                    'error_type': type(exc).__name__, 'candidates': []}, []
            return index, local, items[0], updates
        with ThreadPoolExecutor(max_workers=workers) as pool:
            tasks = [pool.submit(resolve, index, row) for index, row in enumerate(rows)]
            for completed, task in enumerate(as_completed(tasks), 1):
                index, local, item, updates = task.result()
                rows[index][CLINICAL_FIELD] = local.get(CLINICAL_FIELD)
                audit[index] = item
                changes.extend(updates)
                if completed == 1 or completed % 25 == 0 or completed == len(rows):
                    progress = {'completed': completed, 'total': len(rows), 'filled_records': len(changes),
                                'http_requests': client.calls, 'cache_hits': client.cache_hits,
                                'status_counts': dict(Counter(r['status'] for r in audit if r))}
                    progress_path = getattr(args, '_clinvar_progress', None)
                    if progress_path:
                        enrich.atomic_json(progress_path, progress)
                    print(json.dumps({'clinvar_progress': progress}), flush=True)
        return audit, changes, clinvar_report(client, audit, changes)
    audit, changes = [], []
    client = None if args.skip_clinvar else getattr(args, '_clinvar_client', None)
    if client is None and not args.skip_clinvar:
        client = ClinVarClient(args, remaining)
    original = {r['id']: r for r in original_rows}
    for row in rows:
        item = {'id': row['id'], 'PMID': row['PMID'], 'status': '',
                'lookup_performed': False, 'candidates': []}
        audit.append(item)
        if client is None:
            item['status'] = 'DISABLED'
            continue
        if not enrich.missing(row.get(CLINICAL_FIELD)):
            item['status'] = 'EXISTING_VALUE_PRESERVED'
            continue
        gene = enrich.text(row.get('Gene Name'))
        if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', gene):
            item['status'] = 'UNRESOLVED_GENE'
            continue
        texts = [enrich.text(original.get(row['id'], {}).get('Mutation Information')),
                 enrich.text(row.get('Mutation Information'))]
        for member in sources.members(row['id']):
            raw = sources.records.get(member, {})
            for field in ('mutation_information', 'mutation_identity'):
                if isinstance(raw.get(field), str):
                    texts.append(raw[field])
        hgvs = sorted({m[0] for value in texts for m in NUCLEOTIDE_HGVS.finditer(value)})
        formatted, _ = display.format_record(row)
        protein = enrich.text(formatted.get('Mutation Information'))
        if ',' in protein or len({protein_change(p) for p in enrich.N.P_RE.findall(' '.join(texts))}) > 1:
            item['status'] = 'COMPLEX_VARIANT_REQUIRES_REVIEW'
            continue
        if hgvs:
            variants = hgvs
        elif enrich.N.P_RE.fullmatch(protein):
            variants = [protein_change(protein)]
        else:
            item['status'] = 'MISSING_PRECISE_VARIANT'
            continue
        term = '"' + gene + '"[gene] AND (' + ' OR '.join('"' + v + '"[varnam]' for v in variants) + ')'
        item.update(query=term, nucleotide_hgvs=hgvs, lookup_performed=True)
        payload, status = client.request('esearch.fcgi', term=term, retmode='json', retmax='500')
        if status != 'OK':
            item['status'] = status
            continue
        result = json.loads(payload)['esearchresult']
        try:
            count, ids = int(result['count']), result['idlist']
            if count < 0 or not isinstance(ids, list) or any(not isinstance(v, str) or not v.isdigit() for v in ids):
                raise ValueError('Invalid search identifiers')
        except (ValueError, TypeError, KeyError):
            item['status'] = 'INVALID_SEARCH_RESPONSE'
            continue
        if count != len(ids) or len(ids) != len(set(ids)):
            item['status'] = 'TRUNCATED_SEARCH'
            item['search_matches'] = count
            continue
        item['search_matches'] = count
        if not ids:
            item['status'] = 'NO_MATCH'
            continue
        candidates = []
        for start in range(0, len(ids), 40):
            payload, status = client.request('efetch.fcgi', id=','.join(ids[start:start+40]), rettype='vcv',
                                             retmode='xml', is_variationid='true')
            if status != 'OK':
                break
            candidates.extend(parse_clinvar_vcv(payload))
        if status != 'OK':
            item['status'] = status
            continue
        item['candidates'] = candidates
        if {v['variation_id'] for v in candidates} != set(ids):
            item['status'] = 'RETURNED_ID_MISMATCH_OR_UNSUPPORTED_RECORD'
            continue
        matched = [v for v in candidates if gene.casefold() in {g.casefold() for g in v['genes']}
                   and v.get('record_type', 'SINGLE_ALLELE') == 'SINGLE_ALLELE'
                   and (set(hgvs) <= set(v['nucleotide_hgvs']) if hgvs else
                        any(protein_change(h) == protein_change(protein) for h in v['protein_hgvs']))]
        item['exact_variant_candidates'] = len(matched)
        if not matched or (hgvs and len(matched) > 1):
            item['status'] = 'AMBIGUOUS_MATCH' if matched else 'NO_EXACT_VARIANT_MATCH'
            continue
        if hgvs and protein and not any(protein_change(h) == protein_change(protein)
                                       for h in matched[0]['protein_hgvs']):
            item['status'] = 'PROTEIN_CONSEQUENCE_MISMATCH'
            continue
        if not hgvs:
            verified, unresolved, reference_checks = [], [], []
            for candidate in matched:
                reference_status, reference_accession = clinvar_reference_match(row, candidate, protein, client)
                reference_checks.append({'variation_id': candidate['variation_id'],
                    'status': reference_status, 'accession': reference_accession})
                if reference_status in {'UNIPROT_ACCESSION_MATCH', 'REFSEQ_FULL_SEQUENCE_MATCH'}:
                    verified.append(candidate)
                elif reference_status not in {'REFERENCE_ACCESSION_CONFLICT', 'REFERENCE_RESIDUE_MISMATCH', 'REFERENCE_SEQUENCE_MISMATCH'}:
                    unresolved.append(reference_status)
            item.update(reference_checks=reference_checks, match_level='PROTEIN',
                        reference_verified_candidates=len(verified))
            if unresolved or not verified:
                item['status'] = unresolved[0] if unresolved else reference_checks[0]['status']
                continue
            matched = verified
        else:
            item['match_level'] = 'NUCLEOTIDE'
        labels = [enrich.text(v.get('classification')) for v in matched]
        if any(clinical_label_missing(label) for label in labels):
            item['status'] = 'NO_GERMLINE_CLASSIFICATION'
            continue
        if len({label.casefold() for label in labels}) != 1:
            item['status'] = 'AMBIGUOUS_CLASSIFICATION'
            continue
        label = labels[0]
        variation_ids = [v['variation_id'] for v in matched]
        urls = [v['url'] for v in matched]
        before = row.get(CLINICAL_FIELD)
        row[CLINICAL_FIELD] = label
        item.update(status='MATCHED_EXACT_ALLELE' if hgvs else
                    'MATCHED_PROTEIN_VARIANT' if len(matched) == 1 else 'MATCHED_PROTEIN_CONSENSUS',
                    variation_ids=variation_ids, classification=label)
        changes.append({'id': row['id'], 'PMID': row['PMID'], 'field': CLINICAL_FIELD,
                        'before': before, 'after': label,
                        'rule': 'CLINVAR_EXACT_GENE_VERSIONED_NUCLEOTIDE_HGVS' if hgvs else 'CLINVAR_EXACT_GENE_PROTEIN_HGVS_REFERENCE_VERIFIED',
                        'match_level': item['match_level'], 'variation_ids': variation_ids, 'urls': urls})
    return audit, changes, clinvar_report(client, audit, changes)


def refresh_clinvar(args):
    """Refresh clinical annotations of a previously completed JSON export.

    Preserve previously checked rows and display fields byte-for-byte as values;
    a displayed shorthand is not a new sequence-validation input.
    """
    src, out = args.input.resolve(), args.output_dir.resolve()
    if src.suffix.lower() != '.json':
        raise ValueError('--clinvar-only requires the completed_records.json companion file.')
    if out == src.parent or src.is_relative_to(out):
        raise ValueError('Output must be separate from the input and its directory.')
    if out.exists() and any(out.iterdir()):
        if not getattr(args, 'resume', False) or any(p.name != 'work' for p in out.iterdir()):
            raise ValueError('Use a new output directory; --resume accepts only an unfinished ClinVar work directory.')
    input_hash = enrich.digest(src)
    data = enrich.read_json(src)
    query_rows, embedded = load_records_input(src, args.sheet)
    if not isinstance(data, dict) or any(
            data.get('eligibility', {}).get(r['id'], {}).get('enrichment_complete') is not True
            for r in query_rows):
        raise ValueError('--clinvar-only requires a previous completed export with completeness provenance.')
    rows = [{h: r.get(h) for h in COLUMNS} for r in data['records']]
    sources = prepare_sources(query_rows, embedded, args)
    args._clinvar_progress = out / 'work/clinvar_progress.json'
    clinical_audit, changes, clinical_report = enrich_clinvar(rows, query_rows, sources, args, args.max_requests)
    if enrich.digest(src) != input_hash or any(enrich.digest(Path(name)) != expected
                                            for name, expected in sources.hashes.items()):
        raise ValueError('Input or source file changed during ClinVar refresh.')
    # Preserve every unrelated value and every original record, including fields
    # that were intentionally left blank by the earlier formatting step.
    if any(before.get(h) != after.get(h) for before, after in zip(data['records'], rows)
           for h in COLUMNS if h != CLINICAL_FIELD):
        raise ValueError('ClinVar refresh changed an unrelated field.')
    out.mkdir(parents=True, exist_ok=True)
    audit_dir = out / 'audit'; audit_dir.mkdir()
    enrich.atomic_json(audit_dir / 'clinvar.json', clinical_audit)
    enrich.atomic_json(audit_dir / 'field_changes.json', changes)
    result = copy.deepcopy(data)
    result.update(schema_version=VERSION, records=rows,
                  clinical_refresh={'input_sha256': input_hash, 'clinvar': clinical_report})
    enrich.atomic_json(out / 'completed_records.json', result)
    enrich.write_csv(out / 'completed_records.csv', COLUMNS, rows)
    workbook_report = enrich.export_workbook(out, COLUMNS, rows, [], [], [],
        data_only=True, sheet_name='Records', filename='completed_records.xlsx')
    report = {'version': VERSION, 'mode': 'clinvar-only', 'input_sha256': input_hash,
        'input_records': len(rows), 'exported_complete_records': len(rows),
        'http_requests': clinical_report['http_requests'], 'clinvar': clinical_report,
        'columns': len(COLUMNS), 'worksheet_count': 1,
        'unrelated_fields_preserved': True, 'source_sha256': sources.hashes,
        **workbook_report}
    enrich.atomic_json(out / 'run_report.json', report)
    print(json.dumps({k: v for k, v in report.items() if k != 'source_sha256'}, ensure_ascii=False, indent=2))
    return report



def run(args):
    if getattr(args, 'clinvar_only', False):
        return refresh_clinvar(args)
    src = args.input.resolve(); out = args.output_dir.resolve()
    if out == src.parent or src.is_relative_to(out):
        raise ValueError('Output must be separate from the input and its directory.')
    if out.exists() and any(out.iterdir()):
        raise ValueError('Use a new output directory.')
    input_hash = enrich.digest(src)
    rows, embedded = load_records_input(src, args.sheet)
    sources = prepare_sources(rows, embedded, args)
    work = out / 'work'; work.mkdir(parents=True)
    client = enrich.UniProtClient(args.uniprot_cache, args.online, args.max_requests,
        args.timeout, args.retries, args.interval, args.cache_max_age_days, args.include_unreviewed)
    try:
        updated, changes, audit, details = enrich.process(rows, sources, client)
    finally:
        client.close()
    first = work / 'database_enrichment'; first.mkdir()
    enrich.atomic_json(first / 'enriched_records.json', {'columns': COLUMNS, 'records': updated})
    enrich.atomic_json(first / 'record_audit.json', audit)
    enrich.atomic_json(first / 'other_source_details.json', details)
    enrich.atomic_json(first / 'field_changes.json', changes)
    workset = build_recovery_workset(updated, audit, sources, embedded, first / 'enriched_records.json')
    workset_path = work / 'recovery_workset.json'; enrich.atomic_json(workset_path, workset)
    recovery_calls = 0
    if workset['records']:
        recovered = work / 'identity_recovery'
        recovery.run(SimpleNamespace(input_dir=first, workset=workset_path,
            hgnc=args.hgnc, old_cache=args.uniprot_cache, cache_dir=args.reference_cache,
            output_dir=recovered, online=args.online, xlsx=False, force=False,
            max_requests=None if args.max_requests is None else max(0, args.max_requests - client.calls)))
        updated = enrich.read_json(recovered / 'enriched_records.json')['records']
        changes += enrich.read_json(recovered / 'field_changes.json')
        audit = enrich.read_json(recovered / 'record_audit.json')
        recovery_calls = enrich.read_json(recovered / 'run_report.json')['http_requests']
    checks = [completeness(row) for row in updated]
    audit_by_id = {r['id']: r for r in audit}
    for check in checks:
        item = audit_by_id[check['id']]
        if 'CONFLICT' in item['uniprot_status']:
            check['validation_errors'].append('IDENTITY_CONFLICT')
        if any(field + ':EXISTING_VALUE_DIFFERS' in item['notes'] for field in enrich.ENRICH_FIELDS):
            check['validation_errors'].append('EXISTING_DATABASE_VALUE_DIFFERS')
        check['exported'] = not check['missing_fields'] and not check['validation_errors']
    selected_ids = {r['id'] for r in checks if r['exported']}
    selected = [r for r in updated if r['id'] in selected_ids]
    args._clinvar_progress = work / 'clinvar_progress.json'
    remaining = None if args.max_requests is None else max(0, args.max_requests - client.calls - recovery_calls)
    clinical_audit, clinical_changes, clinical_report = enrich_clinvar(selected, rows, sources, args, remaining)
    changes += clinical_changes
    # Completeness and sequence checks above always see the machine/original
    # descriptions. Compact strings are a final presentation step only.
    selected, display_audit = display.format_records(selected)
    if enrich.digest(src) != input_hash:
        raise ValueError('Input changed during enrichment.')
    for name, expected in sources.hashes.items():
        if enrich.digest(Path(name)) != expected:
            raise ValueError('Historical source changed: ' + name)
    audit_dir = out / 'audit'; audit_dir.mkdir()
    enrich.atomic_json(audit_dir / 'completeness.json', checks)
    enrich.atomic_json(audit_dir / 'record_display.json', display_audit)
    enrich.atomic_json(audit_dir / 'clinvar.json', clinical_audit)
    enrich.atomic_json(audit_dir / 'field_changes.json', changes)
    enrich.atomic_json(audit_dir / 'record_audit.json', audit)
    enrich.atomic_json(audit_dir / 'all_records.json', {'columns': COLUMNS, 'records': updated})
    enrich.write_csv(audit_dir / 'completeness.csv', ['id', 'PMID', 'missing_fields', 'validation_errors', 'database_location_required', 'exported'], checks)
    result = {'schema_version': VERSION, 'tier': 'A', 'columns': COLUMNS, 'records': selected,
              'display_originals': {r['id']: r['original_values'] for r in display_audit},
              'eligibility': {r['id']: {'stage_status': 'STAGE_A', 'enrichment_complete': True} for r in selected},
              'source_members': {rid: embedded[rid] for rid in selected_ids if rid in embedded}}
    enrich.atomic_json(out / 'completed_records.json', result)
    enrich.write_csv(out / 'completed_records.csv', COLUMNS, selected)
    target = out / 'completed_records.xlsx'
    workbook_report = enrich.export_workbook(out, COLUMNS, selected, [], [], [],
        data_only=True, sheet_name='Records', filename=target.name)
    report = {'version': VERSION, 'input_sha256': input_hash, 'input_records': len(rows),
              'exported_complete_records': len(selected), 'incomplete_records': len(rows) - len(selected),
              'required_fields': ['Accession_new', 'Reviewed', 'Wild Sequence_new'],
              'database_location_rule': 'Required only for applicable explicit no-change records',
              'missing_field_counts': dict(Counter(f for r in checks for f in r['missing_fields'])),
              'validation_error_counts': dict(Counter(e for r in checks for e in r['validation_errors'])),
              'http_requests': client.calls + recovery_calls + clinical_report['http_requests'],
              'clinvar': clinical_report, 'exported_tiers': ['A'],
              'workbook_sha256': enrich.digest(target), 'source_sha256': sources.hashes,
              'columns': len(COLUMNS), 'worksheet_count': 1,
              **workbook_report}
    enrich.atomic_json(out / 'run_report.json', report)
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return report


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input', type=Path, required=True, help='Step 05 records.json (recommended) or exported data worksheet')
    p.add_argument('--sheet', choices=DATA_SHEETS)
    p.add_argument('--output-dir', type=Path, required=True)
    p.add_argument('--original-records', type=Path)
    p.add_argument('--identity-input', type=Path, help='Original source workset when input is an exported workbook')
    p.add_argument('--source-audit', type=Path)
    p.add_argument('--source-root', type=Path, help='Root of the original project containing source-bound paper files')
    p.add_argument('--dedup-map', type=Path)
    p.add_argument('--uniprot-cache', type=Path, default=ROOT / 'data/06_completed_records/uniprot_cache')
    p.add_argument('--reference-cache', type=Path, default=ROOT / 'data/06_completed_records/reference/cache')
    p.add_argument('--hgnc', type=Path, default=ROOT / 'data/06_completed_records/reference/hgnc_complete_set.json')
    p.add_argument('--clinvar-cache', type=Path, default=ROOT / 'data/06_completed_records/clinvar_cache')
    p.add_argument('--env-file', type=Path, help='NCBI_API_KEY file; otherwise environment, project .env and home .env are checked')
    p.add_argument('--skip-clinvar', action='store_true', help='Skip optional ClinVar classification lookup')
    p.add_argument('--clinvar-only', action='store_true', help='Refresh only clinical classifications in a previously completed JSON export')
    p.add_argument('--clinvar-protein-candidates', action='store_true', help=argparse.SUPPRESS)
    p.add_argument('--clinvar-workers', type=int, default=4, help='Concurrent ClinVar lookups sharing one rate limiter and request budget (default: 4)')
    p.add_argument('--resume', action='store_true', help='Resume an unfinished --clinvar-only run from its request cache; completed exports cannot be overwritten')
    mode = p.add_mutually_exclusive_group()
    mode.add_argument('--online', action='store_true')
    mode.add_argument('--offline', action='store_true')
    p.add_argument('--include-unreviewed', action='store_true')
    p.add_argument('--max-requests', type=int, help='HTTP attempt limit shared by database enrichment, identity recovery and ClinVar')
    p.add_argument('--timeout', type=float, default=30)
    p.add_argument('--retries', type=int, default=2)
    p.add_argument('--interval', type=float, default=0.4)
    p.add_argument('--cache-max-age-days', type=float, default=30)
    return p


def main():
    p = parser(); args = p.parse_args()
    try:
        if args.max_requests is not None and args.max_requests < 0:
            raise ValueError('--max-requests must be nonnegative.')
        if args.retries < 0:
            raise ValueError('--retries must be nonnegative.')
        if not 1 <= args.clinvar_workers <= 32:
            raise ValueError('--clinvar-workers must be between 1 and 32.')
        if args.resume and not args.clinvar_only:
            raise ValueError('--resume requires --clinvar-only.')
        for field in ('timeout', 'cache_max_age_days', 'interval'):
            value = getattr(args, field)
            if not math.isfinite(value) or value < 0 or (field == 'timeout' and value == 0):
                raise ValueError('--' + field.replace('_', '-') + ' has an invalid value.')
        run(args)
    except (OSError, ValueError, KeyError, IndexError) as exc:
        p.exit(2, f'Error: {exc}\n')


if __name__ == '__main__':
    main()
