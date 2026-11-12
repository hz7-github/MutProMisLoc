"""Clean records, adjudicate source evidence and export data tables.

Offline by default. Use --online for new model requests. All intermediate
records and decisions remain in work/ and audit/, outside the final tables.
Mutation/localization display and source-backed metadata recovery are integrated
into this step through the existing runtime modules.
Step 06 handles ClinVar classification after resolving protein references.
"""
from __future__ import annotations

import argparse
import asyncio
import json
import sys
from pathlib import Path
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts._internal import adjudication_workset as worksets
from scripts._internal import adjudication_merge as merger
from scripts._internal import evidence_reconciliation as evidence
from scripts._internal import record_cleaning as cleaning
from scripts._internal import source_adjudication as judge
from scripts._internal import stage_export as exporter


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--preset', choices=['original', 'expanded'], default='expanded')
    p.add_argument('--input', type=Path, help='Raw Step 04 extraction JSON')
    p.add_argument('--metadata', type=Path, help='CSV with PMID and Title')
    p.add_argument('--md-dir', type=Path)
    p.add_argument('--workset', '--identity-input', dest='workset', type=Path)
    p.add_argument('--base', type=Path, help='Saved first-round adjudications; skips new initial requests')
    p.add_argument('--followup', type=Path, action='append', default=[])
    p.add_argument('--context-flags', type=Path)
    p.add_argument('--event-decisions', type=Path)
    p.add_argument('--audit', type=Path, help='Re-export a saved, hash-bound final audit snapshot')
    p.add_argument('--pipeline-root', type=Path, help='Cohort root for --audit')
    p.add_argument('--baseline', type=Path, help='Reference for event decisions; excluded from the exported table')
    p.add_argument('--output-dir', type=Path)
    p.add_argument('--online', action='store_true', help='Allow new paid model requests')
    p.add_argument('--max-requests', type=int, help='Maximum newly selected request batches across all rounds')
    p.add_argument('--repair-rounds', type=int, default=1)
    p.add_argument('--response-dir', type=Path, help='Reuse an existing compatible first-round response cache')
    p.add_argument('--env-file', type=Path)
    p.add_argument('--api-key', nargs='?', const=True)
    p.add_argument('--model', default='deepseek-flash')
    p.add_argument('--base-url', default='https://api.deepseek.com')
    p.add_argument('--thinking', choices=['enabled', 'disabled'], default='enabled')
    p.add_argument('--concurrency', type=int, choices=range(1, 9), default=4)
    p.add_argument('--output-tokens', type=int, default=64000)
    p.add_argument('--policy-profile', choices=['legacy', judge.OPTIONAL_CELL_PROFILE], default=judge.OPTIONAL_CELL_PROFILE)
    return p


def run(args):
    if args.max_requests is not None and args.max_requests < 0 or args.repair_rounds < 0:
        raise ValueError('Request and repair limits must be nonnegative.')
    if args.max_requests and not args.online:
        raise ValueError('--max-requests above zero requires --online.')
    if args.audit and (args.base or args.followup or args.event_decisions or args.online):
        raise ValueError('--audit replays an existing snapshot; use --base for new merges.')
    if args.followup and not args.base:
        raise ValueError('--followup requires --base.')
    out = (args.output_dir or ROOT / 'data/05_records' / args.preset).resolve()
    if out.exists() and any(out.iterdir()):
        raise ValueError('Use a new output directory; saved inputs and caches are reusable.')
    metadata = args.metadata or ROOT / 'data/01_get_pubmed_metadata' / f'papers_{args.preset}.csv'
    md_dir = args.md_dir or ROOT / 'data/03_mineru_pdf_to_md' / args.preset / 'clean'
    if args.audit:
        exporter.run(SimpleNamespace(audit=args.audit, output_dir=out,
            pipeline_root=args.pipeline_root, identity_input=args.workset,
            metadata=metadata, baseline=args.baseline, context_flags=args.context_flags))
        return exporter.read(out / 'export_manifest.json')
    if args.base and not args.workset:
        raise ValueError('--base requires --workset to bind original records and source files.')
    work = out / 'work'; work.mkdir(parents=True)
    if args.workset:
        identities = exporter.read(args.workset)
    else:
        extraction = args.input or ROOT / 'data/04_extract_deepseek' / args.preset / 'extracted_records.json'
        cleaning.unified_clean(SimpleNamespace(input=extraction, metadata=metadata,
            md_dir=md_dir, output_dir=work / 'cleaning', force=False))
        identities = worksets.prepare_dataset(exporter.read(extraction)['records'], md_dir.resolve())
    if not identities or len({r['id'] for r in identities}) != len(identities):
        raise ValueError('A nonempty workset with unique IDs is required.')
    for item in identities:
        path = (ROOT / item['source_path']).resolve()
        if not path.is_relative_to(ROOT) or judge.sha(path) != item['source_sha256']:
            raise ValueError('Workset source changed or is outside the project: ' + item['id'])
    workset = work / 'workset.json'; judge.save(workset, identities)
    by_id = {r['id']: r for r in identities}
    budget = args.max_requests if args.online else 0
    calls = 0

    def adjudicate(path, directory):
        nonlocal budget, calls
        options = vars(args).copy()
        options.update(input=path, output_dir=directory, max_requests=budget,
            run_name=directory.name, plan_only=False, retry_invalid=True,
            retry_existing_only=False, retry_output_tokens=args.output_tokens)
        asyncio.run(judge.run(SimpleNamespace(**options)))
        report = exporter.read(directory / 'adjudication_report.json')
        # Count selected batches, including failed attempts, against the run budget.
        used = report['new_successful_requests'] + report['new_failed_requests']
        calls += used
        if budget is not None:
            budget = max(0, budget - used)
        return directory / 'adjudications.json'

    def reconcile(path, label):
        rows = exporter.read(path)
        for row in rows:
            identity = by_id.get(row['id'])
            if not identity or str(row['PMID']) != str(identity['PMID']) or row['source_sha256'] != identity['source_sha256']:
                raise ValueError('Adjudication does not match the workset: ' + row['id'])
            row['source_path'] = identity['source_path']
        bound = work / f'{label}_source_bound.json'; judge.save(bound, rows)
        directory = work / label
        evidence.run(SimpleNamespace(input=bound, output_dir=directory, preset=args.preset))
        return directory / 'reconciled_adjudications.json'

    initial = args.base or adjudicate(workset, args.response_dir or work / 'initial_responses')
    versions = [reconcile(initial, 'initial_checked')]
    versions.extend(reconcile(path, f'followup_{i}_checked') for i, path in enumerate(args.followup, 1))
    for index in range(args.repair_rounds):
        current = merger.merge_versions([(str(p), exporter.read(p)) for p in versions], {}, allow_empty=True)
        repairs = [r for r in current if r['validation_errors']]
        if not repairs or budget == 0:
            break
        repair_input = work / f'repair_{index + 1}_workset.json'
        judge.save(repair_input, worksets.prepare_repairs(identities, repairs, {}))
        response = adjudicate(repair_input, work / f'repair_{index + 1}_responses')
        versions.append(reconcile(response, f'repair_{index + 1}_checked'))
    merged_dir = work / 'merged'
    merger.run(SimpleNamespace(base=versions[0], followup=versions[1:],
        audit_flags=[args.context_flags] if args.context_flags else [],
        output_dir=merged_dir, identity_input=workset, baseline=args.baseline,
        event_decisions=args.event_decisions, preset=args.preset, md_dir=md_dir,
        allow_empty=True))
    ledger = merged_dir / 'merged_adjudications.json'
    snapshot = work / 'export_snapshot.json'
    judge.save(snapshot, {'cohorts': [{'cohort': 'merged', 'source': str(ledger.resolve()), 'sha256': judge.sha(ledger)}]})
    exporter.run(SimpleNamespace(audit=snapshot, output_dir=out,
        pipeline_root=out, identity_input=workset, metadata=metadata,
        baseline=args.baseline, context_flags=None, allow_work_dir=True))
    report = exporter.read(out / 'export_manifest.json')
    report.update(pipeline_version='records-pipeline-1.1.0', new_model_request_batches=calls)
    judge.save(out / 'export_manifest.json', report)
    return report


def main():
    p = parser(); args = p.parse_args()
    try:
        run(args)
    except (OSError, ValueError, KeyError) as exc:
        p.exit(2, f'Error: {exc}\n')


if __name__ == '__main__':
    main()
