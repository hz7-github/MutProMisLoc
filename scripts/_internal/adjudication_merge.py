"""Merge source adjudication revisions with explicit provenance and audit holds.

Latest results supersede earlier proposals even when they become unresolved or C.
This is a review ledger, not a final scientific release or an event deduplicator.
"""
from __future__ import annotations
import argparse
import collections
import copy
import importlib.util
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location('merge_judge', ROOT / 'scripts/_internal/source_adjudication.py')
judge = importlib.util.module_from_spec(spec)
spec.loader.exec_module(judge)


def merge_versions(versions, audit_flags, allow_empty=False):
    if not versions or (not versions[0][1] and not allow_empty):
        raise ValueError('A nonempty base cohort is required.')
    if not versions[0][1] and any(rows for _, rows in versions[1:]):
        raise ValueError('A nonempty follow-up requires a nonempty base cohort.')
    merged, history, order = {}, collections.defaultdict(list), []
    for label, rows in versions:
        ids = [r['id'] for r in rows]
        if len(ids) != len(set(ids)):
            raise ValueError('Duplicate IDs in input: ' + label)
        for row in rows:
            rid = row['id']
            if rid in merged:
                if merged[rid]['source_sha256'] != row['source_sha256'] or merged[rid]['PMID'] != row['PMID']:
                    raise ValueError('Source binding changed: ' + rid)
            elif merged and label != versions[0][0]:
                raise ValueError('Follow-up contains an ID outside the base cohort: ' + rid)
            if rid not in merged:
                order.append(rid)
            history[rid].append({'input': label, 'tier': row['review']['tier'],
                                 'validation_errors': row['validation_errors'],
                                 'response_sha256': row.get('response_sha256'),
                                 'review': copy.deepcopy(row['review'])})
            merged[rid] = copy.deepcopy(row)
    output = []
    for rid in order:
        row = merged[rid]
        row['revision_history'] = history[rid]
        row['source_audit_flags'] = audit_flags.get(rid, row.get('source_audit_flags', []))
        row['release_eligible'] = False
        row['ledger_status'] = ('EVIDENCE_REPAIR_REQUIRED' if row['validation_errors'] else
                                'SOURCE_CONTEXT_AUDIT_REQUIRED' if row['source_audit_flags'] else
                                'PROPOSAL_READY_FOR_QUALITY_AUDIT')
        output.append(row)
    return output


def event_collision_groups(entries):
    """Find possible repeated variants; never collapse conditions automatically."""
    grouped = collections.defaultdict(list)
    for entry in entries:
        key = (str(entry['PMID']), entry['variant_key'])
        if entry['variant_key']:
            grouped[key].append(entry)
    groups = []
    for (pmid, variant), members in grouped.items():
        if len(members) < 2 or not any(m['cohort'] == 'current' for m in members):
            continue
        contexts = {json.dumps(m['experimental_context'], sort_keys=True, ensure_ascii=False) for m in members}
        genes = {m['gene_name'].strip().casefold() for m in members}
        groups.append({'PMID': pmid, 'variant_key': variant, 'members': members,
                       'comparison': 'IDENTICAL_CONTEXT_CANDIDATE' if len(contexts) == 1 and len(genes) == 1 else 'CONTEXT_OR_IDENTITY_COMPARISON_REQUIRED',
                       'automatically_merged': False,
                       'rule': 'Same paper and normalized variant is a candidate key, not proof of the same experimental event.'})
    return groups


def build_event_audit(output, identity_rows, baseline):
    normalizer = judge.module('event_variant_normalizer', '_internal/record_cleaning.py')
    identities = {r['id']: r['original_record'] for r in identity_rows}
    entries = []
    def variant_key(value, kind):
        converted, _ = normalizer.normalize_mutation(value, kind)
        text = str(converted or '').strip()
        if text.casefold() in {'unknown', 'not_reported', 'none', ''}:
            return ''
        return re.sub(r'\s+', ' ', text).casefold()
    for row in output:
        if row['review']['tier'] not in {'A', 'B'}:
            continue
        original = identities[row['id']]
        entries.append({'id': row['id'], 'PMID': row['PMID'], 'cohort': 'current',
                        'gene_name': str(original.get('gene_name') or ''),
                        'variant_key': variant_key(original.get('mutation_information'), original.get('mutation_type')),
                        'experimental_context': row['review'].get('annotation', {}),
                        'tier': row['review']['tier']})
    for row in baseline:
        value = row['row']
        entries.append({'id': value['id'], 'PMID': str(value['PMID']), 'cohort': 'frozen_baseline',
                        'gene_name': str(value.get('Gene Name') or ''),
                        'variant_key': variant_key(value.get('Mutation Information'), value.get('Mutation Type')),
                        'experimental_context': {'frozen_row': value, 'prior_admission': row.get('review', {}).get('admission')},
                        'tier': 'UNREASSESSED_BASELINE'})
    return {'groups': event_collision_groups(entries), 'entries_examined': len(entries),
            'automatic_merges': 0, 'deduplication_complete': False,
            'limitations': 'Candidate blocking can miss aliases and complex variant equivalences. Different hosts, tags, treatments, times, comparators and assay endpoints must remain separate unless source-supported equivalence is established.'}


def attach_event_decisions(rows, decisions, baseline):
    """Attach source-reviewed event links while retaining every original row."""
    known = {r['id']: r['PMID'] for r in rows}
    known.update({r['row']['id']: str(r['row']['PMID']) for r in baseline})
    links, notes, coverage = {}, collections.defaultdict(list), {}
    for decision in decisions:
        ids = decision['record_ids']; pmid = str(decision['PMID'])
        if any(rid not in known or str(known[rid]) != pmid for rid in ids):
            raise ValueError('Event decision crosses papers or references an unknown record.')
        if decision['action'] not in {'MERGE_IN_FINAL_EVENT', 'KEEP_SEPARATE', 'COVERED_BY_EXISTING_EVENTS'}:
            raise ValueError('Unknown event decision action.')
        if not decision.get('evidence') or not decision.get('reason'):
            raise ValueError('Event decisions require source evidence and a reason.')
        if decision['action'] == 'MERGE_IN_FINAL_EVENT':
            target = decision['canonical_event_record_id']
            if target not in ids:
                raise ValueError('Canonical event must belong to the reviewed group.')
            for rid in ids:
                if rid in links and links[rid] != target:
                    raise ValueError('Conflicting event decisions: ' + rid)
                links[rid] = target
        elif decision['action'] == 'COVERED_BY_EXISTING_EVENTS':
            covering = decision.get('covering_record_ids', [])
            redundant = decision.get('redundant_record_ids', [])
            target = decision.get('canonical_event_record_id')
            if (len(covering) < 2 or not redundant or target not in covering
                    or len(set(covering)) != len(covering) or len(set(redundant)) != len(redundant)
                    or set(covering) & set(redundant) or set(covering) | set(redundant) != set(ids)):
                raise ValueError('Composite coverage requires distinct covering and redundant records.')
            for rid in redundant:
                if rid in links and links[rid] != target:
                    raise ValueError('Conflicting composite event decisions: ' + rid)
                links[rid] = target
                coverage[rid] = covering
        for rid in ids:
            notes[rid].append(decision)
    if any(links.get(target, target) != target for target in links.values()):
        raise ValueError('Event links must point directly to retained canonical records.')
    for covering in coverage.values():
        if any(links.get(rid, rid) != rid for rid in covering):
            raise ValueError('Composite coverage references a redundant record.')
    for row in rows:
        row['canonical_event_record_id'] = links.get(row['id'], row['id'])
        row['covered_by_event_record_ids'] = coverage.get(row['id'], [row['canonical_event_record_id']])
        row['event_review_decisions'] = notes.get(row['id'], [])
        row['event_link_is_final_admission'] = False
    return rows


def run(args):
    paths = [args.base] + args.followup
    if (args.output_dir / 'merged_adjudications.json').resolve() in {p.resolve() for p in paths}:
        raise ValueError('Output must not overwrite an input ledger.')
    versions = [(str(p), json.loads(p.read_text(encoding='utf-8'))) for p in paths]
    flags = collections.defaultdict(list)
    for path in args.audit_flags:
        data = json.loads(path.read_text(encoding='utf-8'))
        for row in data['records']:
            flags[row['id']].append({'source': str(path), **row})
    output = merge_versions(versions, flags, getattr(args, 'allow_empty', False))
    baseline = json.loads(args.baseline.read_text(encoding='utf-8'))['records'] if args.baseline else []
    if args.event_decisions:
        decisions = json.loads(args.event_decisions.read_text(encoding='utf-8'))['decisions']
        pipeline = parser = judge.module('event_source_units', '_internal/source_evidence.py')
        checked = {}
        for decision in decisions:
            pmid = str(decision['PMID'])
            explicit = {row['source_path'] for row in output if str(row['PMID']) == pmid and row.get('source_path')}
            if len(explicit) > 1:
                raise ValueError('Ambiguous source paths for event decision: ' + pmid)
            if explicit:
                path = (ROOT / next(iter(explicit))).resolve()
            else:
                directory = getattr(args, 'md_dir', None) or ROOT / 'data/03_mineru_pdf_to_md' / args.preset / 'clean'
                path = (directory / (pmid + '.md')).resolve()
            if not path.is_relative_to(ROOT):
                raise ValueError('Event source must remain within the project.')
            if judge.sha(path) != decision['source_sha256']:
                raise ValueError('Event decision source changed: ' + pmid)
            if pmid not in checked:
                checked[pmid], _ = pipeline.complete_text_units(path.read_text(encoding='utf-8-sig'), parser)
            for ref in decision['evidence']:
                if not ref.get('quote') or ref['quote'] not in checked[pmid].get(ref['source_id'], {}).get('text', ''):
                    raise ValueError('Event decision quote not matched: ' + pmid)
        output = attach_event_decisions(output, decisions, baseline)
    unknown = set(flags) - {r['id'] for r in output}
    if unknown:
        raise ValueError('Audit flags reference unknown IDs: ' + ', '.join(sorted(unknown)))
    transitions = collections.Counter()
    changes = []
    for row in output:
        history = row['revision_history']
        if len(history) > 1:
            old, new = history[0]['review'], history[-1]['review']
            transitions[old['tier'] + '->' + new['tier']] += 1
            changes.append({'id': row['id'], 'changed_review_fields': [k for k in sorted(set(old) | set(new)) if old.get(k) != new.get(k)],
                            'previous_validation_errors': history[0]['validation_errors'],
                            'current_validation_errors': row['validation_errors']})
    report = {'records': len(output), 'revised_records': len(changes),
              'tier_counts': dict(collections.Counter(r['review']['tier'] for r in output)),
              'ledger_status_counts': dict(collections.Counter(r['ledger_status'] for r in output)),
              'tier_transitions': dict(transitions),
              'input_sha256': {str(p): judge.sha(p) for p in paths + args.audit_flags},
              'script_sha256': judge.sha(Path(__file__)),
              'event_decisions_sha256': judge.sha(args.event_decisions) if args.event_decisions else None,
              'linked_duplicate_records_retained_in_ledger': sum(r.get('canonical_event_record_id', r['id']) != r['id'] for r in output),
              'release_complete': False, 'human_verified': False, 'event_deduplication_complete': False}
    judge.save(args.output_dir / 'merged_adjudications.json', output)
    judge.save(args.output_dir / 'revision_changes.json', changes)
    judge.save(args.output_dir / 'merge_report.json', report)
    if args.identity_input:
        audit = build_event_audit(output, json.loads(args.identity_input.read_text(encoding='utf-8')),
                                  baseline)
        audit['identity_input_sha256'] = judge.sha(args.identity_input)
        audit['baseline_sha256'] = judge.sha(args.baseline) if args.baseline else None
        audit['ledger_sha256'] = judge.sha(args.output_dir / 'merged_adjudications.json')
        judge.save(args.output_dir / 'event_collision_audit.json', audit)
        report['event_collision_groups'] = len(audit['groups'])
        judge.save(args.output_dir / 'merge_report.json', report)
    print(json.dumps(report, indent=2))


