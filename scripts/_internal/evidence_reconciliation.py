"""Reconcile unambiguous claim names and inventory evidence gaps offline.

Exact quotations found elsewhere are suggestions, never automatic scientific
approval. Original responses and scientific annotations remain unchanged.
"""
from __future__ import annotations
import argparse
import collections
import copy
import json
from pathlib import Path
import importlib.util

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location('source_judge', ROOT / 'scripts/_internal/source_adjudication.py')
judge = importlib.util.module_from_spec(spec)
spec.loader.exec_module(judge)
POLICY = 'offline-evidence-reconciliation-1.1.0'
ANNOTATIONS = {'normal_localization', 'mutated_localization', 'cell_line',
               'experimental_method', 'effect_description', 'assay_locator',
               'treatment', 'construct_context', 'tag_strategy',
               'expression_regime', 'quantitative_result'}


def reconcile(review, units, normalizer, policy_profile='legacy'):
    revised = copy.deepcopy(review)
    changes, suggestions = [], []
    for index, ref in enumerate(revised.get('evidence', [])):
        if not isinstance(ref, dict):
            continue
        claim = ref.get('claim', '')
        for prefix, allowed in [('checks.', set(judge.CHECKS)), ('annotation.', ANNOTATIONS)]:
            if isinstance(claim, str) and claim.startswith(prefix) and claim[len(prefix):] in allowed:
                ref['claim'] = claim[len(prefix):]
                changes.append({'evidence_index': index, 'old_claim': claim,
                                'new_claim': ref['claim'], 'rule': 'explicit_schema_namespace'})
                break
        quote = ref.get('quote')
        if not isinstance(quote, str):
            continue
        q = normalizer.normalize_quote_text(quote)
        cited = units.get(ref.get('source_id'), {}).get('text', '')
        if q and q not in normalizer.normalize_quote_text(cited):
            matches = [key for key, unit in units.items()
                       if q in normalizer.normalize_quote_text(unit['text'])]
            suggestions.append({'evidence_index': index, 'claim': ref.get('claim'),
                                'cited_source_id': ref.get('source_id'),
                                'exact_match_source_ids': matches,
                                'action': 'CHECK_EXPERIMENT_BINDING' if matches else 'RECHECK_QUOTATION_AGAINST_SOURCE',
                                'automatically_applied': False})
    errors = judge.validate(revised, units, normalizer, policy_profile)
    # A concrete measured entity is essential even for a conditional observation.
    if revised.get('tier') in {'A', 'B'} and revised.get('checks', {}).get('observed_entity_binding') not in {'SUPPORTED', 'LIMITED'}:
        errors.append('observed_entity_binding_unresolved')
    return revised, changes, suggestions, list(dict.fromkeys(errors))


def run(args):
    inputs = args.input.resolve()
    if args.output_dir.resolve() == inputs.parent:
        raise ValueError('Use a separate output directory to preserve original adjudications.')
    rows = json.loads(inputs.read_text(encoding='utf-8'))
    pipeline = parser = normalizer = judge.module('reconcile_evidence', '_internal/source_evidence.py')
    cache, source_hashes, output = {}, {}, []
    for row in rows:
        pmid = row['PMID']
        if pmid not in cache:
            source = (ROOT / row['source_path']).resolve() if row.get('source_path') else ROOT / 'data/03_mineru_pdf_to_md' / args.preset / 'clean' / (pmid + '.md')
            if not source.resolve().is_relative_to(ROOT):
                raise ValueError('Source must remain within the project.')
            source_hashes[pmid] = judge.sha(source)
            cache[pmid], _ = pipeline.complete_text_units(source.read_text(encoding='utf-8-sig'), parser)
        if source_hashes[pmid] != row['source_sha256']:
            raise ValueError('Source hash changed: ' + pmid)
        review, changes, suggestions, errors = reconcile(row['review'], cache[pmid], normalizer,
                                                       row.get('policy_profile', 'legacy'))
        item = copy.deepcopy(row)
        item.update(review=review, reconciliation_changes=changes,
                    citation_repair_suggestions=suggestions,
                    previous_validation_errors=row['validation_errors'], validation_errors=errors,
                    status='REQUIRES_EVIDENCE_REPAIR' if errors else 'MODEL_SUPPORTED_' + review['tier'],
                    source_audit_complete=False)
        item['evidence_match_audit']=[{'evidence_index':i,'claim':ref.get('claim'),
                                      'source_id':ref.get('source_id'),
                                      'match_mode':judge.quote_match_mode(ref['quote'],cache[pmid][ref['source_id']]['text'],normalizer)}
                                     for i,ref in enumerate(review.get('evidence',[]))
                                     if isinstance(ref,dict) and isinstance(ref.get('quote'),str) and ref.get('source_id') in cache[pmid]]
        output.append(item)
    report = {'policy': POLICY, 'validation_policy': judge.VALIDATION_POLICY, 'records': len(output), 'input_sha256': judge.sha(inputs),
              'script_sha256': judge.sha(Path(__file__)),
              'status_counts': dict(collections.Counter(x['status'] for x in output)),
              'records_with_schema_repairs': sum(bool(x['reconciliation_changes']) for x in output),
              'newly_validation_complete': [x['id'] for x in output if x['previous_validation_errors'] and not x['validation_errors']],
              'records_with_exact_quote_relink_candidates': sum(any(s['exact_match_source_ids'] for s in x['citation_repair_suggestions']) for x in output),
              'remaining_error_counts': dict(collections.Counter(e for x in output for e in x['validation_errors'])),
              'new_api_calls': 0, 'release_complete': False, 'human_verified': False}
    judge.save(args.output_dir / 'reconciled_adjudications.json', output)
    judge.save(args.output_dir / 'evidence_repair_workset.json', [x for x in output if x['validation_errors']])
    judge.save(args.output_dir / 'reconciliation_report.json', report)
    print(json.dumps(report, indent=2))


