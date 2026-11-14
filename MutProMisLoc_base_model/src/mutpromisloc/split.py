"""All-versus-all homology graph with a coverage-of-longer-protein predicate."""
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedGroupKFold

from .data import read_fasta, load_samples, write_fasta, sha256, write_json, new_output


def is_homolog(hit):
    """Inclusive identity>=40%, coordinate coverage of the longer protein>=80%."""
    nident, alnlen = int(hit['nident']), int(hit['alnlen'])
    qlen, tlen = int(hit['qlen']), int(hit['tlen'])
    if alnlen <= 0 or min(qlen, tlen) <= 0 or 5 * nident < 2 * alnlen:
        return False
    qspan = abs(int(hit['qend']) - int(hit['qstart'])) + 1
    tspan = abs(int(hit['tend']) - int(hit['tstart'])) + 1
    if qlen > tlen:
        return 5 * qspan >= 4 * qlen
    if tlen > qlen:
        return 5 * tspan >= 4 * tlen
    return 5 * min(qspan, tspan) >= 4 * qlen


def component_groups(sequences, samples, edges):
    parent = {key: key for key in sequences}

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    def union(a, b):
        ra, rb = find(a), find(b)
        parent[max(ra, rb)] = min(ra, rb)

    exact = {}
    for key, seq in sequences.items():
        if seq in exact:
            union(key, exact[seq])
        exact[seq] = key
    protein_groups = {}
    for r in samples.to_dict('records'):
        union(r['wt_id'], r['mut_id'])
        name = r.get('protein_group', '')
        if pd.notna(name) and str(name).strip():
            if name in protein_groups:
                union(r['wt_id'], protein_groups[name])
            protein_groups[name] = r['wt_id']
    for a, b in edges:
        union(a, b)
    return {key: find(key) for key in parent}


def assign_partitions(meta, seed=2026):
    if meta.cluster_id.nunique() < 7:
        raise ValueError('Need at least 7 independent clusters for 4:1 split and five-fold development CV; collect more proteins')
    y, groups = meta.label.to_numpy(), meta.cluster_id.to_numpy()
    outer = StratifiedGroupKFold(5, shuffle=True, random_state=seed)
    development, test = next(outer.split(np.zeros(len(y)), y, groups))
    out = meta.copy().reset_index(drop=True)
    out['split'] = 'development'
    out.loc[test, 'split'] = 'test'
    out['cv_fold'] = -1
    if len(set(groups[development])) < 5:
        raise ValueError('Development partition has fewer than 5 independent clusters')
    inner = StratifiedGroupKFold(5, shuffle=True, random_state=seed)
    for fold, (a, b) in enumerate(inner.split(np.zeros(len(development)), y[development], groups[development])):
        if set(groups[development[a]]) & set(groups[development[b]]):
            raise AssertionError('CV cluster overlap')
        out.loc[development[b], 'cv_fold'] = fold
        for indices in [development[a], development[b]]:
            if set(y[indices]) != {0, 1}:
                raise ValueError('A development training/validation fold lacks a class; need more independently clustered labeled data')
    for indices in [development, test]:
        if set(y[indices]) != {0, 1}:
            raise ValueError('Development or test partition lacks a class; need more independently clustered data')
    if set(groups[development]) & set(groups[test]):
        raise AssertionError('Outer cluster overlap')
    return out


def prepare(fasta, samples, out, mmseqs='mmseqs', threads=8, seed=2026):
    seqs = read_fasta(fasta)
    meta = load_samples(samples, seqs)
    used = set(meta.wt_id) | set(meta.mut_id)
    seqs = {k: v for k, v in seqs.items() if k in used}
    out = new_output(out)
    # Canonical tool IDs avoid FASTA header interpretation and reuse exact duplicates.
    unique = list(dict.fromkeys(seqs.values()))
    seq_to_tool = {seq: f'p{i:07d}' for i, seq in enumerate(unique)}
    original_for_tool = {}
    for name, seq in seqs.items():
        original_for_tool.setdefault(seq_to_tool[seq], name)
    write_fasta(out / 'unique.fa', {seq_to_tool[s]: s for s in unique})
    database, hits_db = out / 'seqdb', out / 'alndb'
    version = subprocess.check_output([str(mmseqs), 'version'], text=True).strip()
    n = str(len(unique))
    columns = 'query,target,nident,alnlen,qlen,tlen,qstart,qend,tstart,tend'
    commands = [
        [str(mmseqs), 'createdb', str(out / 'unique.fa'), str(database), '--dbtype', '1', '--shuffle', '0', '-v', '1'],
        [str(mmseqs), 'search', str(database), str(database), str(hits_db), str(out / 'tmp'),
         '--exhaustive-search', '1', '--min-seq-id', '0', '-c', '0', '--seq-id-mode', '0',
         '--alignment-mode', '3', '-a', '1', '--max-seqs', n, '--max-accept', n, '--max-rejected', n,
         '-e', '1000000', '-s', '7.5', '--threads', str(threads), '-v', '1'],
        [str(mmseqs), 'convertalis', str(database), str(database), str(hits_db), str(out / 'alignments.tsv'),
         '--format-output', columns, '--threads', str(threads), '-v', '1']]
    for i, command in enumerate(commands):
        print('MMseqs2', ' '.join(command), flush=True)
        with (out / f'mmseqs_{i}.log').open('w') as log:
            subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
    hits = pd.read_csv(out / 'alignments.tsv', sep='\t', names=columns.split(','))
    if hits.empty:
        raise ValueError('MMseqs2 returned no alignments; see mmseqs_*.log')
    # convertalis nident needs a stored alignment backtrace; without -a it can
    # silently report zeros even for a full-length self alignment.
    if not hits.nident.gt(0).any():
        raise ValueError('MMseqs2 returned no identical-residue counts; require alignment backtrace (-a 1)')
    self_hits = hits[hits['query'] == hits['target']]
    if not self_hits.empty and self_hits.nident.eq(0).any():
        raise ValueError('MMseqs2 self alignments have zero identity counts; check alignment backtrace')
    pairs = []
    for h in hits.to_dict('records'):
        if h['query'] != h['target'] and is_homolog(h):
            a, b = original_for_tool[h['query']], original_for_tool[h['target']]
            pairs.append((a, b))
    groups = component_groups(seqs, meta, pairs)
    meta['cluster_id'] = meta.wt_id.map(groups)
    partitions = assign_partitions(meta, seed)
    partitions.to_csv(out / 'splits.csv', index=False)
    pd.DataFrame(pairs, columns=['sequence_a', 'sequence_b']).to_csv(out / 'homology_edges.csv', index=False)
    pd.DataFrame([{'sequence_id': k, 'tool_id': seq_to_tool[seqs[k]], 'cluster_id': groups[k]} for k in seqs]).to_csv(out / 'sequence_groups.csv', index=False)
    for name in ['development', 'test']:
        sub = partitions[partitions.split == name]
        ids = set(sub.wt_id) | set(sub.mut_id)
        write_fasta(out / f'{name}.fa', {key: seqs[key] for key in seqs if key in ids})
    summary = {}
    for name, sub in [('development', partitions[partitions.split == 'development']), ('test', partitions[partitions.split == 'test'])]:
        summary[name] = dict(samples=len(sub), positive=int(sub.label.sum()), negative=int((sub.label == 0).sum()), clusters=int(sub.cluster_id.nunique()))
    summary['test_fraction'] = summary['test']['samples'] / len(meta)
    summary['total_clusters'] = int(partitions.cluster_id.nunique())
    summary['homology_edges_crossing_partitions'] = 0
    # Independent edge and duplicate check against partition assignments.
    seq_partition = {}
    for r in partitions.to_dict('records'):
        for key in [r['wt_id'], r['mut_id']]:
            if key in seq_partition and seq_partition[key] != r['split']:
                raise AssertionError('Sequence ID crosses partitions')
            seq_partition[key] = r['split']
    if any(seq_partition[a] != seq_partition[b] for a,b in pairs):
        raise AssertionError('Homology edge crosses partitions')
    manifest = dict(fasta=str(Path(fasta).resolve()), samples=str(Path(samples).resolve()),
                    fasta_sha256=sha256(fasta), samples_sha256=sha256(samples), seed=seed,
                    split_sha256=sha256(out / 'splits.csv'), mmseqs_version=version, commands=commands,
                    identity='nident/alnlen >= 0.40 (including gaps in alignment denominator)',
                    coverage='(abs(end-start)+1)/length >= 0.80 on the longer protein; equal lengths require both',
                    grouping='connected components of observed homology edges + exact duplicates + WT/MUT pair + optional protein_group',
                    search_note='MMseqs2 exhaustive all-vs-all search; edges reflect the alignments found, not a proof that all possible local alignments were enumerated',
                    outer='first fold of seeded StratifiedGroupKFold(5); approximately 80/20 by rows',
                    inner='seeded StratifiedGroupKFold(5) within development only', summary=summary)
    write_json(out / 'manifest.json', manifest)
    write_json(out / 'summary.json', summary)
    from .release import export_existing
    export_existing(fasta, out / 'splits.csv', out / 'dataset')
    print(summary, flush=True)
