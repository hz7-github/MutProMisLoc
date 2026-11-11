"""Export an offline, traceable stage workbook in the existing 26-column layout.

This is a candidate snapshot, not a final release. No network calls are made.
Source wording and uncertainty are retained; frozen releases are never modified.
"""
from __future__ import annotations
import argparse
import collections
import csv
import hashlib
import importlib.util
import json
from io import BytesIO
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
VERSION = 'records-summary-1.2.0'


def module(name, filename):
    spec = importlib.util.spec_from_file_location(name, ROOT / 'scripts' / filename)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def read(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def path_label(path):
    path = path.resolve()
    return path.relative_to(ROOT).as_posix() if path.is_relative_to(ROOT) else path.as_posix()


def known(value):
    if value is None or str(value).strip().casefold() in {'', 'unknown', 'not_reported', 'not reported', 'null', 'n/a', 'na', '*', 'undetermined'}:
        return None
    return str(value).strip()


def stage_row(item, original, title, normalizer):
    review = item['review']; annotation = review.get('annotation', {})
    row = {h: None for h in normalizer.BASE_COLUMNS + [normalizer.DETAIL]}
    before = known(annotation.get('normal_localization'))
    after = known(annotation.get('mutated_localization'))
    outcome = review.get('outcome_class')
    if outcome == 'normal_experimentally_confirmed':
        phenotype, provenance = 'No change', 'EXPERIMENTAL_NO_CHANGE'
    elif outcome == 'observation_only':
        phenotype, provenance = 'Observation only', 'EXPERIMENTAL_OBSERVATION_ONLY'
    elif outcome == 'mislocalization':
        phenotype, provenance = (before or '*') + '>' + (after or '*'), 'EXPERIMENTAL_MISLOCALIZATION'
    else:
        raise ValueError('Unresolved outcome cannot enter an A/B stage table.')
    mutation, mutation_note = normalizer.normalize_mutation(original.get('mutation_information'), original.get('mutation_type'))
    mutation_type, type_note = normalizer.normalize_type(original.get('mutation_information'), original.get('mutation_type'))
    cell, _, cell_notes = normalizer.normalize_cell(annotation.get('cell_line'))
    details = [known(annotation.get('effect_description'))]
    for key, label in [('construct_context', 'Construct context'), ('treatment', 'Reported condition'),
                       ('assay_locator', 'Source locator'), ('quantitative_result', 'Reported measurement')]:
        if known(annotation.get(key)):
            details.append(label + ': ' + annotation[key])
    if review.get('limitations'):
        details.append('Limitations: ' + '; '.join(review['limitations']))
    row.update({'id': item['id'], 'PMID': str(item['PMID']), 'Gene Name': original.get('gene_name'),
                'Mutation Type': mutation_type, 'Mutation Information': mutation,
                'Mislocalization Phenotype': phenotype, normalizer.BASE_COLUMNS[6]: before,
                normalizer.BASE_COLUMNS[7]: after, 'Cell Line': known(cell),
                'Disease Information': known(original.get('disease_information')), 'Title': title,
                'Experimental_Method': known(annotation.get('experimental_method')),
                'Tag_Strategy': known(annotation.get('tag_strategy')),
                'Expression_Regime': known(annotation.get('expression_regime')),
                'Outcome_Provenance': provenance,
                normalizer.DETAIL: '\n'.join(x for x in details if x)})
    # Metadata enums not adjudicated by the current schema stay blank. Do not
    # infer quantitative extent, tag position or evidence type from keywords.
    return row, {'mutation_note': mutation_note, 'type_note': type_note,
                 'cell_notes': cell_notes, 'localization_wording_preserved': True}


def classify(item, context_flags):
    if item.get('pending_adjudication'):
        return 'PENDING_ADJUDICATION'
    if item['review'].get('material_required'):
        return 'PENDING_SOURCE_MATERIAL'
    if item['validation_errors']:
        return 'PENDING_VALIDATION'
    if context_flags:
        return 'PENDING_SOURCE_CONTEXT'
    if item['review']['tier'] == 'C':
        return 'C_' + item['review'].get('c_disposition', 'UNRESOLVED')
    if item.get('canonical_event_record_id', item['id']) != item['id']:
        return 'LINKED_DUPLICATE_RETAINED_IN_AUDIT'
    return 'STAGE_' + item['review']['tier']


def decision_label(status):
    return {
        'STAGE_A': 'A类候选（待最终质量核查）',
        'STAGE_B': 'B类条件性观察候选',
        'C_OUT_OF_SCOPE': '模型建议排除（保留审计）',
        'C_UNRESOLVED': '待补科学证据',
        'PENDING_VALIDATION': '待修复证据引用或校验问题',
        'PENDING_SOURCE_CONTEXT': '待核原文语境',
        'PENDING_SOURCE_MATERIAL': '待补必要原文图表或材料',
        'PENDING_ADJUDICATION': '尚未完成模型核验',
        'LINKED_DUPLICATE_RETAINED_IN_AUDIT': '关联重复项（保留审计）',
    }.get(status, status)


def run(args):
    from openpyxl import Workbook, load_workbook
    n = module('stage_normalizer', '_internal/record_cleaning.py')
    root = (getattr(args, 'pipeline_root', None) or ROOT / 'data/12_final_pipeline').resolve()
    out = args.output_dir.resolve()
    if out.exists() and any(p.name != 'work' or not getattr(args, 'allow_work_dir', False) for p in out.iterdir()):
        raise ValueError('Use a new versioned output directory.')
    identity_input = getattr(args, 'identity_input', None)
    goal_path = root / 'goal_manifest.json'
    goal = read(goal_path) if goal_path.exists() else {'baseline_files': {}}
    if not identity_input and not goal_path.exists():
        raise ValueError('Legacy export requires its goal manifest; dataset export requires --identity-input.')
    for name, expected in goal['baseline_files'].items():
        if digest(ROOT / name) != expected:
            raise ValueError('Frozen baseline changed: ' + name)
    audit_path = args.audit.resolve()
    audit = read(audit_path); items = []; inputs = {path_label(audit_path): digest(audit_path)}
    for cohort in audit['cohorts']:
        path = root / cohort['source']
        if digest(path) != cohort['sha256']:
            raise ValueError('Cohort changed since the audited snapshot: ' + cohort['cohort'])
        inputs[path_label(path)] = digest(path)
        for item in read(path):
            items.append({**item, 'cohort': cohort['cohort'], 'ledger_path': path_label(path)})
    if len({r['id'] for r in items}) != len(items):
        raise ValueError('Cohort IDs overlap.')
    wanted = {r['id'] for r in items}; originals = {}; counts = collections.Counter()
    if identity_input:
        identity_input = identity_input.resolve()
        inputs[path_label(identity_input)] = digest(identity_input)
        identity_rows = read(identity_input)
        identities = {r['id']: r for r in identity_rows}
        if len(identities) != len(identity_rows) or not wanted <= set(identities):
            raise ValueError('Duplicate or missing dataset identities.')
        checked_sources = {}
        for item in items:
            identity = identities[item['id']]
            if str(item['PMID']) != str(identity['PMID']) or item['source_sha256'] != identity['source_sha256']:
                raise ValueError('Adjudication and input source binding differ: ' + item['id'])
            if identity.get('source_path'):
                source = (ROOT / identity['source_path']).resolve()
                if not source.is_relative_to(ROOT):
                    raise ValueError('Source must remain within the project.')
                if source not in checked_sources:
                    checked_sources[source] = digest(source)
                if checked_sources[source] != identity['source_sha256']:
                    raise ValueError('Source changed after adjudication: ' + item['id'])
        for rid, identity in identities.items():
            originals[rid] = identity['original_record']
            if rid not in wanted:
                items.append({'id': rid, 'PMID': str(identity['PMID']), 'source_sha256': identity['source_sha256'],
                              'source_path': identity.get('source_path'), 'pending_adjudication': True,
                              'validation_errors': [], 'cohort': 'not_yet_adjudicated', 'ledger_path': None,
                              'review': {'tier': None, 'annotation': {}, 'limitations': [], 'evidence': []}})
        wanted = set(identities)
    else:
        original_path = ROOT / 'data/04_extract_deepseek/expanded/extracted_records.json'
        inputs[path_label(original_path)] = digest(original_path)
        for original in read(original_path)['records']:
            pmid = str(original['pmid']); counts[pmid] += 1
            rid = f'{pmid}:{counts[pmid]:04d}'
            if rid in wanted:
                originals[rid] = original
        for relative in ['recovered_source_adjudication_input_v1.json', 'recovered_zero_output_adjudication_input_v1.json']:
            path = root / relative; inputs[path_label(path)] = digest(path)
            for row in read(path):
                if row['id'] in wanted:
                    if row['id'] in originals:
                        raise ValueError('New and original IDs overlap.')
                    originals[row['id']] = row['original_record']
    if set(originals) != wanted:
        raise ValueError('Missing original identity records.')
    metadata = {}; metadata_path = (getattr(args, 'metadata', None) or ROOT / 'data/01_get_pubmed_metadata/papers_expanded_available.csv').resolve()
    inputs[path_label(metadata_path)] = digest(metadata_path)
    with metadata_path.open(encoding='utf-8-sig', newline='') as stream:
        for row in csv.DictReader(stream):
            if row['PMID'] in metadata and metadata[row['PMID']] != row['Title']:
                raise ValueError('Multiple cached titles for one PMID.')
            metadata[row['PMID']] = row['Title']
    flags = collections.defaultdict(list)
    flag_path = getattr(args, 'context_flags', None)
    if not flag_path and not identity_input:
        flag_path = root / 'remaining_saved_pass_v1/source_context_audit_flags.json'
    if flag_path:
        flag_path = flag_path.resolve()
        inputs[path_label(flag_path)] = digest(flag_path)
        for row in read(flag_path)['records']:
            flags[row['id']].append(row['reason'])
    if set(flags) - wanted:
        raise ValueError('Context flags reference records outside this export.')
    columns = n.BASE_COLUMNS[:6] + [n.DETAIL] + n.BASE_COLUMNS[6:]
    assert len(columns) == 26
    tables = {'Records': []}
    audit_rows = []; canonical = []; status_counts = collections.Counter()
    for item in items:
        rid = item['id']; review = item['review']; original = originals[rid]
        context = flags[rid] + [r['reason'] for r in item.get('source_audit_flags', [])]
        status = classify(item, context); status_counts[status] += 1
        event_notes = [r['reason'] for r in item.get('event_review_decisions', [])]
        refs = [{'claim': e.get('claim'), 'source_id': e.get('source_id')} for e in review.get('evidence', []) if isinstance(e, dict)]
        annotation = review.get('annotation', {})
        output_row, normalization = None, None
        if status == 'STAGE_A':
            output_row, normalization = stage_row(item, original, metadata.get(str(item['PMID'])), n)
            tables['Records'].append([output_row[h] for h in columns])
            assert all(output_row[h] is None for h in n.BLANK_COLUMNS)
        audit_rows.append([rid, str(item['PMID']), review['tier'], status, decision_label(status),
                           review.get('reason'), original.get('gene_name'),
                           original.get('mutation_information'), review.get('outcome_class'), review.get('entity_scope'),
                           '; '.join(review.get('limitations', [])), '; '.join(item['validation_errors']),
                           '; '.join(context), item.get('canonical_event_record_id', rid), '; '.join(event_notes),
                           annotation.get('assay_locator'), annotation.get('treatment'), annotation.get('construct_context'),
                           item['ledger_path'], item['source_sha256'],
                           '; '.join(f"{e['claim']}:{e['source_id']}" for e in refs),
                           'https://pubmed.ncbi.nlm.nih.gov/' + str(item['PMID']) + '/'])
        canonical.append({'id': rid, 'PMID': str(item['PMID']), 'stage_status': status,
                          'row': output_row, 'normalization': normalization,
                          'source_ledger': item['ledger_path'], 'source_sha256': item['source_sha256'],
                          'model_tier': review['tier'], 'model_reason': review.get('reason'),
                          'decision_label': decision_label(status),
                          'limitations': review.get('limitations'), 'validation_errors': item['validation_errors'],
                          'source_context_flags': context, 'event_notes': event_notes,
                          'canonical_event_record_id': item.get('canonical_event_record_id', rid),
                          'covered_by_event_record_ids': item.get('covered_by_event_record_ids',
                                                                   [item.get('canonical_event_record_id', rid)]),
                          'evidence_references': refs, 'final_admission': False})
    baseline_path = getattr(args, 'baseline', None)
    if not baseline_path and not identity_input:
        baseline_path = ROOT / 'data/10_targeted_followup/expanded/final_records.json'
    baseline = []
    if baseline_path:
        baseline_path = baseline_path.resolve()
        baseline = read(baseline_path)['records']; inputs[path_label(baseline_path)] = digest(baseline_path)
    audit_headers = ['id', 'PMID', 'Model tier', 'Stage status', '核验去向', 'Model reason', 'Gene Name', 'Original mutation',
                     'Outcome class', 'Entity scope', 'Limitations', 'Validation errors', 'Source context holds',
                     'Canonical event id', 'Event review notes', 'Assay locator', 'Reported condition',
                     'Construct context', 'Source ledger', 'Source SHA256', 'Evidence references', 'PubMed URL']
    out.mkdir(parents=True, exist_ok=True)
    book = Workbook(write_only=True); expected = {}
    for name, rows in tables.items():
        expected[name] = n.write_sheet(book, name, columns, rows, plain=True)
    workbook_path = out / 'records.xlsx'; temporary = out / 'records.tmp.xlsx'
    book.save(temporary)
    check = load_workbook(BytesIO(temporary.read_bytes()), read_only=True, data_only=False)
    for name, value in expected.items():
        if n.data_digest(check[name].values) != value:
            raise ValueError('Workbook readback mismatch: ' + name)
    for name in tables:
        if len(list(check[name].iter_rows(min_row=1, max_row=1, values_only=True))[0]) != 26:
            raise ValueError('Wrong stage column count.')
    check.close(); temporary.replace(workbook_path)
    audit_dir = out / 'audit'; audit_dir.mkdir()
    (audit_dir / 'stage_records.json').write_text(json.dumps(canonical, ensure_ascii=False, indent=2), encoding='utf-8', newline='\n')
    with (audit_dir / 'record_audit.csv').open('w', encoding='utf-8-sig', newline='') as stream:
        writer = csv.writer(stream); writer.writerow(audit_headers)
        writer.writerows([["'" + v if isinstance(v, str) and v.lstrip().startswith(('=', '+', '-', '@')) else v for v in row] for row in audit_rows])
    selected = [r for r in canonical if r['stage_status'] == 'STAGE_A']
    identity_by_id = {r['id']: r for r in identity_rows} if identity_input else {}
    source_members = {}
    for entry in selected:
        members = [r for r in items if r.get('canonical_event_record_id', r['id']) == entry['id']]
        source_members[entry['id']] = [dict(id=r['id'], record=originals[r['id']],
            source_path=identity_by_id.get(r['id'], {}).get('source_path', r.get('source_path')),
            source_sha256=r['source_sha256']) for r in members]
    payload = {'schema_version': VERSION, 'tier': 'A', 'columns': columns,
               'records': [r['row'] for r in selected], 'source_members': source_members,
               'eligibility': {r['id']: {'stage_status': r['stage_status'], 'source_sha256': r['source_sha256']} for r in selected}}
    (out / 'records.json').write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding='utf-8', newline='\n')
    with (out / 'records.csv').open('w', encoding='utf-8-sig', newline='') as stream:
        writer = csv.writer(stream); writer.writerow(columns)
        writer.writerows([["'" + v if isinstance(v, str) and v.lstrip().startswith(('=', '+', '-', '@')) else v for v in row] for row in tables['Records']])
    report = {'version': VERSION, 'records_accounted': len(items), 'status_counts': dict(status_counts),
              'exported_records': len(tables['Records']), 'exported_tiers': ['A'],
              'B_exported_rows': 0, 'C_exported_rows': 0,
              'columns': 26, 'frozen_baseline_rows': len(baseline), 'input_sha256': inputs,
              'workbook_sha256': digest(workbook_path), 'script_sha256': digest(Path(__file__)),
              'new_api_calls': 0, 'readback_verified': True, 'release_complete': False,
              'worksheet_count': 1, 'data_only': True, 'cell_fills': False,
              'human_verified': False, 'unresolved_event_duplicates_may_remain': True,
              'pending_adjudication_records': status_counts.get('PENDING_ADJUDICATION', 0)}
    (out / 'export_manifest.json').write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8', newline='\n')
    for name, value in goal['baseline_files'].items():
        assert digest(ROOT / name) == value
    if baseline_path:
        assert digest(baseline_path) == inputs[path_label(baseline_path)]
    print(json.dumps({k: v for k, v in report.items() if k != 'input_sha256'}, indent=2))


