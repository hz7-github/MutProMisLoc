"""Prepare reproducible source-adjudication inputs without network requests."""
from __future__ import annotations
import argparse
import collections
import copy
import csv
import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location('workset_judge', ROOT / 'scripts/_internal/source_adjudication.py')
judge = importlib.util.module_from_spec(spec)
spec.loader.exec_module(judge)


def read(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def prepare_repairs(original, repairs, notes):
    indexed = {r['id']: r for r in original}
    if len(indexed) != len(original):
        raise ValueError('Duplicate original record IDs.')
    output = []
    for row in repairs:
        if not row['validation_errors']:
            continue
        item = copy.deepcopy(indexed[row['id']])
        if item['source_sha256'] != row['source_sha256']:
            raise ValueError('Repair source hash differs from original input.')
        item['remaining_errors'] = row['validation_errors']
        item['adjudication_feedback'] = {
            'previous_proposal': row['review'],
            'validation_errors': row['validation_errors'],
            'citation_suggestions': row.get('citation_repair_suggestions', []),
            'source_context_inspection': notes.get(row['id']),
            'instruction': 'Verify against source; repair evidence or narrow the observation honestly. Do not force acceptance.'}
        output.append(item)
    return output


def prepare_dataset(extracted, md_dir, paper_ids=None):
    """Prepare every record in the selected paper scope, without scientific filtering."""
    counts, sources, output = collections.Counter(), {}, []
    for record in extracted:
        pmid = str(record['pmid'])
        counts[pmid] += 1
        if paper_ids is not None and pmid not in paper_ids:
            continue
        if pmid not in sources:
            path = (md_dir / (pmid + '.md')).resolve()
            if not path.is_relative_to(ROOT):
                raise ValueError('Dataset Markdown must remain within the project.')
            sources[pmid] = (path.relative_to(ROOT).as_posix(), judge.sha(path))
        source_path, source_hash = sources[pmid]
        output.append({'id': f'{pmid}:{counts[pmid]:04d}', 'PMID': pmid,
                       'original_record': record, 'revised_review': {'corrected_fields': record},
                       'remaining_errors': ['Raw extraction requires source adjudication; missingness and uncertainty retained.'],
                       'source_sha256': source_hash, 'source_path': source_path,
                       'origin': 'DATASET_EXTRACTION'})
    return output


def run(args):
    root = ROOT / 'data/12_final_pipeline'
    bindings = {}
    def bound(path):
        bindings[str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)] = judge.sha(path)
        return read(path)
    if args.mode == 'dataset':
        extracted = bound(args.extraction_input)['records']
        paper_ids = None
        if args.paper_list:
            bindings[args.paper_list.resolve().as_posix()] = judge.sha(args.paper_list)
            with args.paper_list.open(encoding='utf-8-sig', newline='') as stream:
                reader = csv.DictReader(stream)
                if 'PMID' not in (reader.fieldnames or []):
                    raise ValueError('Paper scope CSV requires a PMID column.')
                ids = [r['PMID'].strip() for r in reader]
                if len(ids) != len(set(ids)) or not all(p.isdigit() for p in ids):
                    raise ValueError('Paper scope requires unique numeric PMIDs.')
                paper_ids = set(ids)
        output = prepare_dataset(extracted, args.md_dir.resolve(), paper_ids)
    elif args.mode == 'repair':
        original = bound(args.base_input)
        repairs = bound(args.repair_input)
        notes = {r['id']: r for r in bound(args.context_notes)['records']} if args.context_notes else {}
        output = prepare_repairs(original, repairs, notes)
    elif args.mode == 'subset':
        output = bound(args.base_input)
    elif args.mode == 'extraction':
        extracted = bound(args.extraction_input)['records']
        overrides = bound(args.source_manifest)
        output, counts = [], collections.Counter()
        for record in extracted:
            pmid = str(record['pmid'])
            counts[pmid] += 1
            info = overrides[pmid]
            source = (args.source_manifest.parent / info['markdown_file']).resolve()
            if not source.is_relative_to(ROOT) or judge.sha(source) != info['markdown_sha256']:
                raise ValueError('Recovered source path or hash mismatch: ' + pmid)
            output.append({'id': f'{pmid}:R{counts[pmid]:04d}', 'PMID': pmid,
                           'original_record': record, 'revised_review': {'corrected_fields': record},
                           'remaining_errors': ['Recovered-source extraction requires source adjudication.'],
                           'source_sha256': info['markdown_sha256'],
                           'source_path': source.relative_to(ROOT).as_posix(),
                           'origin': 'RECOVERED_SOURCE_NEW_CANDIDATE'})
    else:
        worksets = bound(root / 'recovery_worksets/record_worksets.json')
        selected = {r['id']: r for r in worksets if set(args.queue) & set(r['queues'])}
        if not selected:
            raise ValueError('No records match the requested queues.')
        original = bound(ROOT / 'data/04_extract_deepseek' / args.preset / 'extracted_records.json')['records']
        counts, source_hashes, output = collections.Counter(), {}, []
        for record in original:
            pmid = str(record['pmid'])
            counts[pmid] += 1
            rid = f'{pmid}:{counts[pmid]:04d}'
            if rid not in selected:
                continue
            if pmid not in source_hashes:
                source = ROOT / 'data/03_mineru_pdf_to_md' / args.preset / 'clean' / (pmid + '.md')
                source_hashes[pmid] = judge.sha(source)
            output.append({'id': rid, 'PMID': pmid, 'original_record': record,
                           'revised_review': {'corrected_fields': record},
                           'remaining_errors': [selected[rid]['initial_reasons']],
                           'source_sha256': source_hashes[pmid],
                           'recovery_queues': selected[rid]['queues']})
        if set(selected) != {r['id'] for r in output}:
            raise ValueError('Original record IDs and recovery queues do not agree.')
    exclusions = set()
    for path in args.exclude_ids_from:
        entries = bound(path)
        if isinstance(entries, dict):
            entries = entries['records']
        exclusions.update(r['id'] for r in entries)
    excluded_count = sum(r['id'] in exclusions for r in output)
    output = [r for r in output if r['id'] not in exclusions]
    if not output or len({r['id'] for r in output}) != len(output):
        raise ValueError('Empty or duplicate workset.')
    if args.output.exists() and read(args.output) != output:
        raise ValueError('Existing workset differs; choose a new versioned output path.')
    judge.save(args.output, output)
    report = {'mode': args.mode, 'records': len(output), 'papers': len({r['PMID'] for r in output}),
              'queue_filters': args.queue, 'input_sha256': bindings,
              'excluded_existing_ids': excluded_count,
              'reported_barrier_counts': dict(collections.Counter(str(reason) for r in output for reason in r.get('remaining_errors', []))),
              'output_sha256': judge.sha(args.output), 'script_sha256': judge.sha(Path(__file__)),
              'api_calls': 0, 'scientific_decisions_changed': False}
    judge.save(args.output.with_suffix('.manifest.json'), report)
    print(json.dumps(report, indent=2))


