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
VERSION = 'records-summary-1.4.0'


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




def _build_metadata_rules():
    """Recover optional metadata from source-matched saved review evidence only.

    Rules classify verified quotations, not free model prose or raw extraction enums.
    Unknown/unreported information stays blank. These are evidence types, not scores.
    """
    import re

    FIELDS=('Evidence_Type','Quantitative_Reported','Redistribution_Extent','No_Change_Basis','Cell_System_Type')
    SPATIAL=r'\b(?:immunofluorescen\w*|immunostain\w*|confocal|fluorescen\w* microscop\w*|electron microscop\w*|fluorescen\w* imag\w*)\b'
    BIOCHEM=r'\b(?:subcellular fractionation|cell fractionation|membrane fractionation|surface biotinylation|surface accessibility|protease protection)\b'
    LOCATION=r'\b(?:locali[sz]\w*|redistribut\w*|nucle(?:ar|us)|cytoplas\w*|cytosol\w*|membrane|cell surface|colocali[sz]\w*|compartments?)\b'


    def recover(review,units,parser,normalizer):
        values=dict.fromkeys(FIELDS);verified={};proofs={};unmatched=0
        for ref in review.get('evidence',[]):
            if not isinstance(ref,dict):continue
            quote=ref.get('quote');unit=units.get(ref.get('source_id'),{}).get('text','')
            if not isinstance(quote,str) or not quote.strip():continue
            q=parser.normalize_quote_text(quote)
            if not q or q not in parser.normalize_quote_text(unit):unmatched+=1;continue
            claim=ref.get('claim','')
            for prefix in ('annotation.','checks.'):
                if claim.startswith(prefix):claim=claim[len(prefix):]
            verified.setdefault(claim,[]).append(quote)
            proofs.setdefault(claim,[]).append({'claim':claim,'source_id':ref.get('source_id')})
        used={};rules={}
        def put(field,value,rule,claims):
            values[field]=value;rules[field]=rule
            used[field]=[p for claim in claims for p in proofs.get(claim,[])]
        def quotes(*claims):return '\n'.join(q for claim in claims for q in verified.get(claim,[]))
        def has(pattern,text):return bool(re.search(pattern,text,re.I))
        checks=review.get('checks',{});annotation=review.get('annotation',{})
        method=quotes('experimental_method')
        direct=checks.get('direct_localization_measurement')=='SUPPORTED' and bool(verified.get('direct_localization_measurement'))
        # Do not classify a method stated only as explicitly excluded/not used.
        method_parts=[p for p in re.split(r'[.;\n]',method) if not has(r'\b(?:not used|did not use|without using|no imaging|not performed)\b',p)]
        method='\n'.join(method_parts)
        if direct and has(SPATIAL,method):
            put('Evidence_Type','SPATIAL_DIRECT','verified direct-localization check and imaging method quote',('experimental_method','direct_localization_measurement'))
        elif direct and has(BIOCHEM,method):
            put('Evidence_Type','BIOCHEMICAL_LOCALIZATION','verified direct-localization check and biochemical localization method quote',('experimental_method','direct_localization_measurement'))
        quantitative=quotes('quantitative_result')
        if quantitative and has(LOCATION,quantitative):
            negative=has(r'\b(?:not quantified|no quantification|without quantification|quantification was not|not quantitatively)\b',quantitative)
            if negative:
                put('Quantitative_Reported','NO','explicit statement that localization was not quantified',('quantitative_result',))
            elif has(r'\d|\bquantif\w*\b',quantitative) and has(r'%|\b(?:percent\w*|ratio\w*|quantif\w*|intensit\w*|fraction of cells|proportion\w*)\b',quantitative):
                put('Quantitative_Reported','YES','verified localization measurement; abundance/sample-size alone is insufficient',('quantitative_result',))
        outcome=review.get('outcome_class')
        comparison=checks.get('comparator_binding')=='SUPPORTED' and bool(verified.get('comparator_binding'))
        outcome_supported=checks.get('outcome_supported')=='SUPPORTED' and bool(verified.get('outcome_supported'))
        if outcome=='normal_experimentally_confirmed' and comparison and direct and outcome_supported:
            put('Redistribution_Extent','NO_CHANGE','source-supported normal outcome with direct matched comparison',('outcome_supported','comparator_binding','direct_localization_measurement'))
            comparison_quote=quotes('outcome_supported','effect_description','comparator_binding')
            if values['Quantitative_Reported']=='YES' and has(r'\bwild[- ]type\b|\bWT\b|\bmutant\w*\b',quantitative):
                put('No_Change_Basis','QUANTITATIVE_COMPARISON','source-verified localization quantification and matched comparison',('quantitative_result','comparator_binding','outcome_supported'))
            elif has(r'no (?:significant )?(?:change|difference)[^.\n]{0,60}(?:locali[sz]|distribut)|(?:locali[sz]\w*|distribut\w*)[^.\n]{0,60}(?:unchanged|unaffected|identical)|same (?:subcellular )?(?:locali[sz]ation|distribution)',comparison_quote):
                put('No_Change_Basis','EXPLICIT_NO_CHANGE_STATEMENT','explicit source statement bound to the confirmed outcome',('outcome_supported','effect_description','comparator_binding'))
            elif values['Evidence_Type']=='SPATIAL_DIRECT':
                put('No_Change_Basis','QUALITATIVE_COMPARISON','verified imaging and matched normal-outcome comparison',('experimental_method','comparator_binding','outcome_supported'))
        elif outcome=='mislocalization' and outcome_supported:
            result=quotes('effect_description','outcome_supported')
            if has(LOCATION,result):
                if has(r'\b(?:partially|partial)\s+(?:redistribut\w*|relocat\w*|mislocali[sz]\w*|retained|retention)|\bpartial loss of (?:nuclear|membrane|surface|locali[sz])',result):
                    put('Redistribution_Extent','PARTIAL','explicit partial localization outcome quote',('effect_description','outcome_supported'))
                elif has(r'\b(?:completely|complete)\s+(?:redistribut\w*|relocat\w*|mislocali[sz]\w*|retained|retention)|\bcomplete exclusion from (?:the )?(?:nucleus|membrane)',result):
                    put('Redistribution_Extent','COMPLETE','explicit complete localization outcome quote',('effect_description','outcome_supported'))
        # Annotated cell/system must be backed by a source-matched cell-line quote.
        cell=annotation.get('cell_line') or '';cell_quote=quotes('cell_line')
        if cell_quote and cell and str(cell).casefold() not in {'not_reported','unknown'}:
            compact=lambda s:re.sub('[^a-z0-9]','',s.casefold())
            matched=compact(cell) in compact(cell_quote)
            if matched:
                for pattern,value in [(r'\borganoid\w*\b','ORGANOID_3D'),
                                      (r'\bstem[- ]cell[- ]derived\b|\biPSC[- ]derived\b','STEM_CELL_DERIVED'),
                                      (r'\bprimary (?:cells|neurons|fibroblasts|hepatocytes|cultures)\b','PRIMARY_CELLS'),
                                      (r'\btissue sections?\b|\btissue samples?\b','TISSUE'),
                                      (r'\bcell[- ]free\b','CELL_FREE')]:
                    if has(pattern,cell):
                        put('Cell_System_Type',value,'explicit system description in source-bound cell annotation',('cell_line',));break
            if values['Cell_System_Type'] is None:
                canonical,ids,notes=normalizer.normalize_cell(cell)
                if ids:
                    supported=False
                    for name,accession,aliases in normalizer.CELL_LINES:
                        if accession not in ids:continue
                        for alias in [name,*aliases]:
                            pattern=r'(?<![A-Za-z0-9])'+re.escape(alias).replace(r'\-',r'[- ]?')+r'(?![A-Za-z0-9])'
                            if has(pattern,cell_quote):supported=True
                    if supported:put('Cell_System_Type','CELL_LINE','recognized cell-line identity also present in the verified quotation',('cell_line',))
        return values,{'rules':rules,'evidence':used,'unmatched_quotes_ignored':unmatched,
                       'unfilled_fields':[k for k,v in values.items() if v is None],
                       'no_change_basis_applicable':outcome=='normal_experimentally_confirmed'}

    from types import SimpleNamespace
    return SimpleNamespace(recover=recover, FIELDS=FIELDS)

EVIDENCE_METADATA = _build_metadata_rules()


def known(value):
    if value is None or str(value).strip().casefold() in {'', 'unknown', 'not_reported', 'not reported', 'null', 'n/a', 'na', '*', 'undetermined'}:
        return None
    return str(value).strip()


def stage_row(item, original, title, normalizer, source_units=None, source_parser=None):
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
    metadata = EVIDENCE_METADATA
    if source_units is not None and source_parser is not None:
        recovered, metadata_audit = metadata.recover(review, source_units, source_parser, normalizer)
        row.update(recovered)
    else:
        metadata_audit = {'unfilled_fields': list(metadata.FIELDS), 'reason': 'SOURCE_UNITS_UNAVAILABLE'}
    display = normalizer.DISPLAY
    row, display_audit = display.format_record(row)
    return row, {'mutation_note': mutation_note, 'type_note': type_note,
                 'cell_notes': cell_notes, 'localization_wording_preserved': True,
                 'evidence_metadata': metadata_audit, 'record_display': display_audit}


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
    metadata_parser = module('stage_metadata_source_units', '_internal/source_evidence.py')
    metadata_units = {}
    for item in items:
        rid = item['id']; review = item['review']; original = originals[rid]
        context = flags[rid] + [r['reason'] for r in item.get('source_audit_flags', [])]
        status = classify(item, context); status_counts[status] += 1
        event_notes = [r['reason'] for r in item.get('event_review_decisions', [])]
        refs = [{'claim': e.get('claim'), 'source_id': e.get('source_id')} for e in review.get('evidence', []) if isinstance(e, dict)]
        annotation = review.get('annotation', {})
        output_row, normalization = None, None
        if status == 'STAGE_A':
            source = (ROOT / item.get('source_path', f"data/03_mineru_pdf_to_md/expanded/clean/{item['PMID']}.md")).resolve()
            source_units = None
            if source.is_relative_to(ROOT) and source.is_file() and digest(source) == item['source_sha256']:
                if source not in metadata_units:
                    metadata_units[source], _ = metadata_parser.complete_text_units(source.read_text(encoding='utf-8-sig'), metadata_parser)
                source_units = metadata_units[source]
            output_row, normalization = stage_row(item, original, metadata.get(str(item['PMID'])), n, source_units, metadata_parser)
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
    (audit_dir / 'evidence_metadata.json').write_text(json.dumps([
        {'id': r['id'], **r['normalization']['evidence_metadata']}
        for r in canonical if r['normalization']], ensure_ascii=False, indent=2), encoding='utf-8')
    (audit_dir / 'stage_records.json').write_text(json.dumps(canonical, ensure_ascii=False, indent=2), encoding='utf-8', newline='\n')
    with (audit_dir / 'record_audit.csv').open('w', encoding='utf-8-sig', newline='') as stream:
        writer = csv.writer(stream); writer.writerow(audit_headers)
        writer.writerows([["'" + v if isinstance(v, str) and v.lstrip().startswith(('=', '+', '-', '@')) else v for v in row] for row in audit_rows])
    selected = [r for r in canonical if r['stage_status'] == 'STAGE_A']
    (audit_dir / 'record_display.json').write_text(json.dumps([
        r['normalization']['record_display'] for r in selected], ensure_ascii=False, indent=2), encoding='utf-8')
    identity_by_id = {r['id']: r for r in identity_rows} if identity_input else {}
    source_members = {}
    for entry in selected:
        members = [r for r in items if r.get('canonical_event_record_id', r['id']) == entry['id']]
        source_members[entry['id']] = [dict(id=r['id'], record=originals[r['id']],
            source_path=identity_by_id.get(r['id'], {}).get('source_path', r.get('source_path')),
            source_sha256=r['source_sha256']) for r in members]
    payload = {'schema_version': VERSION, 'tier': 'A', 'columns': columns,
               'display_originals': {r['id']: r['normalization']['record_display']['original_values'] for r in selected},
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


