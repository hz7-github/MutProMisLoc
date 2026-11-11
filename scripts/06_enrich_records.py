"""Enrich verified records and export rows with complete target fields.

Requires Accession_new, Reviewed and Wild Sequence_new. Database location is
required only where the existing no-change applicability rule permits it.
Offline by default; incomplete records remain in work/ and audit/.
"""
from __future__ import annotations

import argparse
import copy
import json
import re
import sys
from collections import Counter
from pathlib import Path
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts._internal import identity_recovery as recovery
from scripts._internal import legacy_enrichment as enrich

VERSION = 'complete-records-enrichment-1.1.0'
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
                path = (ROOT / member['source_path']).resolve()
                if not path.is_relative_to(ROOT) or enrich.digest(path) != member.get('source_sha256'):
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


def run(args):
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
    if enrich.digest(src) != input_hash:
        raise ValueError('Input changed during enrichment.')
    for name, expected in sources.hashes.items():
        if enrich.digest(Path(name)) != expected:
            raise ValueError('Historical source changed: ' + name)
    audit_dir = out / 'audit'; audit_dir.mkdir()
    enrich.atomic_json(audit_dir / 'completeness.json', checks)
    enrich.atomic_json(audit_dir / 'field_changes.json', changes)
    enrich.atomic_json(audit_dir / 'record_audit.json', audit)
    enrich.atomic_json(audit_dir / 'all_records.json', {'columns': COLUMNS, 'records': updated})
    enrich.write_csv(audit_dir / 'completeness.csv', ['id', 'PMID', 'missing_fields', 'validation_errors', 'database_location_required', 'exported'], checks)
    result = {'schema_version': VERSION, 'tier': 'A', 'columns': COLUMNS, 'records': selected,
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
              'http_requests': client.calls + recovery_calls, 'exported_tiers': ['A'],
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
    p.add_argument('--dedup-map', type=Path)
    p.add_argument('--uniprot-cache', type=Path, default=ROOT / 'data/06_completed_records/uniprot_cache')
    p.add_argument('--reference-cache', type=Path, default=ROOT / 'data/06_completed_records/reference/cache')
    p.add_argument('--hgnc', type=Path, default=ROOT / 'data/06_completed_records/reference/hgnc_complete_set.json')
    mode = p.add_mutually_exclusive_group()
    mode.add_argument('--online', action='store_true')
    mode.add_argument('--offline', action='store_true')
    p.add_argument('--include-unreviewed', action='store_true')
    p.add_argument('--max-requests', type=int, help='HTTP attempt limit shared by database enrichment and identity recovery')
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
        run(args)
    except (OSError, ValueError, KeyError, IndexError) as exc:
        p.exit(2, f'Error: {exc}\n')


if __name__ == '__main__':
    main()
