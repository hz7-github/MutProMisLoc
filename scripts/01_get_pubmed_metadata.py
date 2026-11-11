#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""01 Retrieve PubMed metadata without downloading full text. Python >= 3.10.

Install: python -m pip install openpyxl
Run: python 01_get_pubmed_metadata.py
Change cutoff: python 01_get_pubmed_metadata.py --end-date 2026-11-03
Custom query: python 01_get_pubmed_metadata.py --query "protein mislocalization mutation"
Expanded query is the default: python 01_get_pubmed_metadata.py
Custom output name: python 01_get_pubmed_metadata.py --output-name papers_custom
Preview without network access: python 01_get_pubmed_metadata.py --plan-only
Each run performs a fresh search; no checkpoint or resume files are written.
Default output: data/01_get_pubmed_metadata/ under the project root, without run subdirectories.
Sample run: --limit 3 (explicitly reported as a subset, not a complete export).
Offline self-test: --self-test

The original query and dates are embedded from 01_pubmed_search_strategy.json;
no legacy scripts or configuration files are required.
Outputs: <output-name>.xlsx and <output-name>.csv only (default: papers_expanded).
Use --output-name without an extension or directory. Search details are included in the
Excel workbook. Metadata is held in memory during the run. No PDF or full-text
download functionality is included. An API key is optional. Environment variables
NCBI_API_KEY and NCBI_EMAIL are supported; API keys are not saved in logs or output.
Searches exceeding 9999 results stop rather than silently truncate the result set.
Documentation: https://www.ncbi.nlm.nih.gov/books/NBK25499/
"""
import argparse
import csv
import json
import os
from pathlib import Path
import re
import tempfile
import time
from datetime import date, datetime, timezone
from urllib.parse import urlencode
from urllib.request import Request, urlopen
import xml.etree.ElementTree as ET

DEFAULT_QUERY = "protein mislocalization mutation"
EXPANDED_QUERY = (
    "protein[Title/Abstract] AND mutation[Title/Abstract] AND "
    "(trafficking[Title/Abstract] OR redistribution[Title/Abstract] OR "
    "retention[Title/Abstract] OR translocation[Title/Abstract])"
)
SEARCH_PRESETS = {"original": DEFAULT_QUERY, "expanded": EXPANDED_QUERY}


def build_query(term, start_date, end_date):
    return f'({term}) AND ("{start_date}"[Date - Publication] : "{end_date}"[Date - Publication])'


HEADERS = ["PMID", "DOI", "PMCID", "Title", "Abstract", "First Author", "Authors",
           "Journal", "Publication Date", "Publication Types", "Language", "Keywords", "MeSH Terms", "PubMed URL", "DOI URL"]


def txt(node, path=None):
    node = node.find(path) if node is not None and path else node
    return "" if node is None else "".join(node.itertext()).strip()


def parse_articles(root):
    result = {}
    for entry in list(root):
        if entry.tag == "PubmedArticle":
            citation = entry.find("MedlineCitation")
            article = citation.find("Article") if citation is not None else None
            id_nodes = entry.findall("PubmedData/ArticleIdList/ArticleId")
            journal = txt(article, "Journal/Title")
            pubdate = article.find("Journal/JournalIssue/PubDate") if article is not None else None
            keywords = [txt(n) for n in citation.findall("KeywordList/Keyword")] if citation is not None else []
            mesh = [txt(n) for n in citation.findall("MeshHeadingList/MeshHeading/DescriptorName")] if citation is not None else []
        elif entry.tag == "PubmedBookArticle":
            citation = article = entry.find("BookDocument")
            id_nodes = article.findall("ArticleIdList/ArticleId") if article is not None else []
            journal, pubdate = txt(article, "Book/BookTitle"), article.find("Book/PubDate") if article is not None else None
            keywords, mesh = [], []
        else:
            continue
        pmid = txt(citation, "PMID")
        if not pmid or article is None:
            continue
        identifiers = {n.get("IdType", "").lower(): txt(n) for n in id_nodes}
        authors = []
        for n in article.findall("AuthorList/Author"):
            name = txt(n, "CollectiveName") or " ".join(filter(None, [txt(n, "ForeName") or txt(n, "Initials"), txt(n, "LastName")]))
            if name:
                authors.append(name)
        abstract = "\n".join((n.get("Label", "") + ": " if n.get("Label") else "") + txt(n) for n in article.findall("Abstract/AbstractText"))
        doi = identifiers.get("doi", "")
        row = [pmid, doi, identifiers.get("pmc", ""), txt(article, "ArticleTitle"), abstract,
               authors[0] if authors else "", "; ".join(authors), journal,
               " ".join(txt(n) for n in pubdate) if pubdate is not None else "",
               "; ".join(txt(n) for n in article.findall("PublicationTypeList/PublicationType")),
               "; ".join(txt(n) for n in article.findall("Language")), "; ".join(keywords), "; ".join(mesh),
               "https://pubmed.ncbi.nlm.nih.gov/" + pmid + "/", "https://doi.org/" + doi if doi else ""]
        if pmid in result:
            raise ValueError("Duplicate PMID in the returned metadata")
        result[pmid] = row
    return result


class Client:
    def __init__(self, timeout, retries):
        self.timeout, self.retries, self.last = timeout, retries, 0
        self.key, self.email = os.environ.get("NCBI_API_KEY", ""), os.environ.get("NCBI_EMAIL", "")

    def request(self, endpoint, **params):
        params.update(db="pubmed", retmode="xml", tool="MutProML_Metadata_01")
        if self.key:
            params["api_key"] = self.key
        if self.email:
            params["email"] = self.email
        for attempt in range(self.retries + 1):
            time.sleep(max(0, 0.36 - (time.monotonic() - self.last)))
            self.last = time.monotonic()
            try:
                url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/" + endpoint + "?" + urlencode(params)
                with urlopen(Request(url, headers={"User-Agent": "MutProML-Metadata/1.0"}), timeout=self.timeout) as f:
                    root = ET.fromstring(f.read())
                if root.find(".//ERROR") is not None or root.find(".//ErrorList") is not None:
                    raise ValueError("NCBI reported a search or record error; check the parameters")
                return root
            except (OSError, ET.ParseError):
                if attempt == self.retries:
                    raise RuntimeError(f"Request failed for {endpoint}; rerun to start a fresh search (request URLs are omitted to protect credentials)") from None
                time.sleep(2 ** attempt)


def validate_output_name(name):
    if (not name or name != name.strip() or name in {".", ".."}
            or name.endswith((".", " ")) or re.search(r'[<>:"/\\|?*\x00-\x1f]', name)
            or name.lower().endswith((".xlsx", ".csv"))
            or name.split(".")[0].upper() in {"CON", "PRN", "AUX", "NUL", *[f"COM{i}" for i in range(1, 10)], *[f"LPT{i}" for i in range(1, 10)]}):
        raise ValueError("output-name must be a valid filename stem without a directory or .xlsx/.csv extension")
    return name


def export(out, rows, report, output_name="papers"):
    output_name = validate_output_name(output_name)
    from openpyxl import Workbook, load_workbook
    from openpyxl.styles import Font, PatternFill, Alignment
    def safe(v):
        return re.sub(r"[\x00-\x08\x0b\x0c\x0e-\x1f]", lambda m: "\\u" + format(ord(m[0]), "04x"), str(v))
    # Reject silent Excel truncation.
    if any(len(safe(v)) > 32767 for r in rows for v in r):
        raise ValueError("Metadata exceeds Excel's cell length limit; export stopped without truncating data")
    wb = Workbook()
    ws = wb.active
    ws.title = "Papers"
    ws.append(HEADERS)
    for row in rows:
        ws.append([safe(v) for v in row])
        for cell in ws[ws.max_row]:
            cell.data_type = "s"
    ws.freeze_panes = "D2"
    ws.auto_filter.ref = ws.dimensions
    for c in ws[1]:
        c.font = Font(bold=True, color="FFFFFF")
        c.fill = PatternFill("solid", fgColor="203E5B")
    for col, width in zip("ABCDEFGHIJKLMNO", [14, 35, 18, 80, 100, 25, 65, 40, 22, 40, 16, 50, 60, 50, 55]):
        ws.column_dimensions[col].width = width
    for row in ws:
        for c in row:
            c.alignment = Alignment(vertical="top", wrap_text=True)
    audit = wb.create_sheet("Search Details")
    audit.append(["Field", "Value"])
    for k, v in report.items():
        audit.append([k, json.dumps(v, ensure_ascii=False) if isinstance(v, (list, dict)) else v])
    audit.column_dimensions["A"].width = 32
    audit.column_dimensions["B"].width = 100
    # Temporary files are removed on failure and never remain as run outputs.
    with tempfile.TemporaryDirectory(prefix=".export_", dir=out) as staging:
        temp = Path(staging) / "papers.xlsx"
        wb.save(temp)
        check = load_workbook(temp, read_only=True)
        try:
            actual = list(check["Papers"].values)
        finally:
            check.close()
        expected = [[None if safe(v) == "" else safe(v) for v in r] for r in rows]
        assert [list(r) for r in actual[1:]] == expected, "XLSX readback mismatch"
        csv_temp = Path(staging) / "papers.csv"
        with csv_temp.open("w", encoding="utf-8-sig", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(HEADERS)
            writer.writerows(rows)
        os.replace(temp, out / f"{output_name}.xlsx")
        os.replace(csv_temp, out / f"{output_name}.csv")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    selection = p.add_mutually_exclusive_group()
    selection.add_argument("--preset", choices=SEARCH_PRESETS, default="expanded", help="Built-in search strategy; default: expanded")
    selection.add_argument("--query", help="Custom query instead of a built-in preset")
    p.add_argument("--start-date", default="1900-01-01")
    p.add_argument("--end-date", default="2025-08-01")
    p.add_argument("--output-dir", type=Path, help="Override data/01_get_pubmed_metadata/ (no automatic subdirectory)")
    p.add_argument("--output-name", type=validate_output_name,
                   help="Output filename stem without extension (default: papers_<preset>)")
    p.add_argument("--limit", type=int, help="Limit exported records for testing; the report identifies the subset")
    p.add_argument("--timeout", type=float, default=30)
    p.add_argument("--retries", type=int, default=2)
    p.add_argument("--plan-only", action="store_true")
    p.add_argument("--self-test", action="store_true")
    a = p.parse_args()
    if a.self_test:
        root = ET.fromstring('<PubmedArticleSet><PubmedArticle><MedlineCitation><PMID>1</PMID><Article><ArticleTitle>A <i>gene</i></ArticleTitle><Abstract><AbstractText Label="BACKGROUND">abc</AbstractText></Abstract><AuthorList><Author><LastName>Smith</LastName><ForeName>Jane</ForeName></Author></AuthorList></Article></MedlineCitation><PubmedData><ArticleIdList><ArticleId IdType="doi">own</ArticleId></ArticleIdList><ReferenceList><ArticleId IdType="doi">wrong</ArticleId></ReferenceList></PubmedData></PubmedArticle><PubmedBookArticle><BookDocument><PMID>2</PMID><ArticleTitle>Book chapter</ArticleTitle><Book><BookTitle>Book</BookTitle></Book></BookDocument></PubmedBookArticle></PubmedArticleSet>')
        rows = parse_articles(root)
        assert len(rows) == 2 and rows['1'][1] == 'own' and rows['1'][3] == 'A gene' and rows['1'][5] == 'Jane Smith'
        assert rows['1'][4] == 'BACKGROUND: abc' and rows['2'][7] == 'Book'
        with tempfile.TemporaryDirectory() as td:
            export(Path(td), list(rows.values()), {"test": True})
            assert {f.name for f in Path(td).iterdir()} == {"papers.xlsx", "papers.csv"}
        with tempfile.TemporaryDirectory() as td:
            export(Path(td), list(rows.values()), {"test": True}, "papers_expanded")
            assert {f.name for f in Path(td).iterdir()} == {"papers_expanded.xlsx", "papers_expanded.csv"}
        for invalid in ("../papers", "folder/papers", "papers.xlsx", "CON", ""):
            try:
                validate_output_name(invalid)
            except ValueError:
                pass
            else:
                raise AssertionError(f"Invalid filename accepted: {invalid!r}")
        print("Self-test passed: first record, book records, inline title tags, abstract sections, reference DOI isolation, and XLSX readback")
        return
    term = a.query if a.query is not None else SEARCH_PRESETS[a.preset or "expanded"]
    preset = "custom" if a.query is not None else a.preset or "expanded"
    a.output_name = a.output_name or f"papers_{preset}"
    if date.fromisoformat(a.start_date) > date.fromisoformat(a.end_date) or not term.strip():
        p.error("Invalid date range or empty query")
    if (a.limit is not None and a.limit < 1) or a.timeout <= 0 or not 0 <= a.retries <= 5:
        p.error("limit and timeout must be positive; retries must be between 0 and 5")
    query = build_query(term, a.start_date, a.end_date)
    if a.plan_only:
        print(query)
        print("Retrieve bibliographic metadata only; do not download full text.")
        return
    import openpyxl  # Check dependencies before making network requests.
    out = a.output_dir or Path(__file__).resolve().parents[1] / "data" / "01_get_pubmed_metadata"
    out.mkdir(parents=True, exist_ok=True)
    if any((out / f"{a.output_name}.{ext}").exists() for ext in ("xlsx", "csv")):
        raise ValueError("Output files already exist; choose a new --output-dir or --output-name")
    client = Client(a.timeout, a.retries)
    root = client.request("esearch.fcgi", term=query, retmax=9999, sort="pub_date")
    count = int(root.findtext("Count", "-1"))
    ids = [txt(n) for n in root.findall("IdList/Id")]
    if count > 9999:
        raise ValueError("More than 9999 results; narrow the date range and search in segments. No incomplete result set was exported")
    if count < 0 or len(ids) != count or len(set(ids)) != count:
        raise ValueError("Incomplete search ID set; stopping rather than silently omitting records")
    snap = dict(query=query, query_translation=root.findtext("QueryTranslation", ""), reported_count=count,
                pmids=ids, searched_at=datetime.now(timezone.utc).isoformat(), warnings=[txt(n) for n in root.findall("WarningList/*")])
    ids = snap['pmids'][:a.limit] if a.limit else snap['pmids']
    records = {}
    missing = [i for i in ids if i not in records]
    for start in range(0, len(missing), 100):
        batch = missing[start:start + 100]
        fetched = parse_articles(client.request("efetch.fcgi", id=','.join(batch)))
        if set(fetched) - set(batch):
            raise ValueError("Unexpected PMID returned; stopping to prevent incorrect metadata assignment")
        records.update(fetched)
        print(f"Metadata retrieved: {sum(i in records for i in ids)}/{len(ids)}", flush=True)
    pending = [i for i in ids if i not in records]
    report = dict(search_preset=preset, query=query, query_translation=snap.get('query_translation'), search_count=snap['reported_count'],
                  selected_count=len(ids), exported_count=sum(i in records for i in ids), missing_pmids=pending,
                  sample_limit=a.limit, complete_search_export=not pending and len(ids)==len(snap['pmids']),
                  full_text_downloaded=False, generated_at=datetime.now(timezone.utc).isoformat())
    export(out, [records[i] for i in ids if i in records], report, a.output_name)
    print(f"Exported: {out / (a.output_name + '.xlsx')}\nMissing metadata: {len(pending)}; full-text downloads: 0")
    if pending:
        raise RuntimeError(f"Some PMIDs lack metadata; see Search Details in {a.output_name}.xlsx. Rerun with a new output directory for a fresh search")


if __name__ == '__main__':
    try:
        main()
    except (OSError, ValueError, RuntimeError, ImportError) as exc:
        print(f"Failed: {exc}")
        raise SystemExit(2)
