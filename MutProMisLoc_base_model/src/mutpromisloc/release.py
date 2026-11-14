"""Publish fixed FASTA/CSV development, test and five-fold data partitions."""
import json
from pathlib import Path

import pandas as pd

from .data import read_fasta, load_samples, write_fasta, sha256, write_json, new_output


def export_dataset(sequences, splits, out, provenance=None):
    """Materialize existing assignments, never calculate or change a split."""
    required = {'sample_id', 'wt_id', 'mut_id', 'position', 'label',
                'cluster_id', 'split', 'cv_fold'}
    if not required <= set(splits):
        raise ValueError(f'Missing split columns: {sorted(required - set(splits))}')
    if splits.empty or not splits.sample_id.is_unique or splits.sample_id.isna().any():
        raise ValueError('Sample IDs must be nonempty and unique')
    if splits.cluster_id.isna().any() or splits.cluster_id.astype(str).str.strip().eq('').any():
        raise ValueError('Each sample must have a cluster ID')
    if set(splits.split) != {'development', 'test'}:
        raise ValueError('Expected nonempty development and test partitions')
    if not splits.label.isin([0, 1]).all():
        raise ValueError('Dataset release requires binary labels')
    development = splits[splits.split == 'development']
    test = splits[splits.split == 'test']
    if set(development.cv_fold) != set(range(5)) or not test.cv_fold.eq(-1).all():
        raise ValueError('Expected five development folds (0..4) and test fold=-1')
    sequence_partitions, cluster_partitions = {}, {}
    for row in splits.itertuples(index=False):
        partition = (row.split, int(row.cv_fold))
        if row.cluster_id in cluster_partitions and cluster_partitions[row.cluster_id] != partition:
            raise ValueError('A cluster crosses test/development or validation folds')
        cluster_partitions[row.cluster_id] = partition
        for key in [row.wt_id, row.mut_id]:
            if key not in sequences:
                raise ValueError(f'Missing sequence ID: {key}')
            sequence = sequences[key]
            if sequence in sequence_partitions and sequence_partitions[sequence] != partition:
                raise ValueError('An identical sequence crosses test/development or validation folds')
            sequence_partitions[sequence] = partition
    subsets = {'development': development, 'test': test}
    for fold in range(5):
        subsets[f'folds/fold_{fold}/train'] = development[development.cv_fold != fold]
        subsets[f'folds/fold_{fold}/validation'] = development[development.cv_fold == fold]
    for name, table in subsets.items():
        if set(table.label) != {0, 1}:
            raise ValueError(f'{name} lacks one of the two classes')

    # Validate all boundaries before writing any files.
    out = new_output(out)
    all_ids = set(splits.wt_id) | set(splits.mut_id)
    write_fasta(out / 'proteins.fa', {k: s for k, s in sequences.items() if k in all_ids})
    sample_columns = [c for c in splits.columns if c not in {'split', 'cv_fold', 'cluster_id'}]
    splits[sample_columns].to_csv(out / 'samples.csv', index=False)
    splits.to_csv(out / 'splits.csv', index=False)
    counts = {}
    for name, table in subsets.items():
        base = out / name
        base.parent.mkdir(parents=True, exist_ok=True)
        table.to_csv(str(base) + '.csv', index=False)
        ids = set(table.wt_id) | set(table.mut_id)
        write_fasta(str(base) + '.fa', {k: s for k, s in sequences.items() if k in ids})
        counts[name] = dict(samples=len(table), positive=int(table.label.sum()),
                            negative=int(table.label.eq(0).sum()),
                            clusters=int(table.cluster_id.nunique()), sequence_records=len(ids))
    (out / 'README.md').write_text(
        '# Fixed protein mutation dataset partitions\n\n'
        'proteins.fa and samples.csv contain all used WT/MUT sequences and sample records. '
        'FASTA IDs are the first token in each header. Each CSV row is one WT/MUT pair; '
        'position is 1-based and label is the study-supplied binary 0/1 label. '
        'Consult the accompanying study for the biological definition and source of each class.\n\n'
        'development.fa / development.csv contain the approximately 80% development set. '
        'test.fa / test.csv contain the fixed approximately 20% held-out test set. '
        'Within folds/fold_0/ through folds/fold_4/, train.fa / train.csv and '
        'validation.fa / validation.csv reproduce each of the five development CV folds. '
        'The test set is shared by all folds and is not used for fitting or early stopping.\n\n'
        'Train and validation are disjoint within each fold, and neither overlaps the test set. '
        'Each development sample is validation exactly once and training in the other four folds; '
        'repetition of development records across different fold folders is intentional. '
        'Each subset FASTA includes both WT and MUT records, with shared IDs written once. '
        'Sequence-record counts can therefore differ from sample counts. Labels and pair '
        'relationships are in the CSV, not encoded in the FASTA headers.\n\n'
        'splits.csv is the authoritative assignment table (sample_id, cluster_id, split, cv_fold). '
        'Files are exported from these assignments without re-randomization. Cluster and exact '
        'sequence boundaries are checked before export. The grouping rule used in the project is '
        'sequence identity >=40% and alignment coverage of the longer protein >=80%, '
        'with connected components and WT/MUT pair links. Dataset manifests retain the source '
        'split hash; the original prepare manifest contains homology search details.\n\n'
        'dataset_manifest.json records per-subset sample/class/cluster counts and file hashes. '
        'Verify file integrity on Linux with: sha256sum -c SHA256SUMS\n\n'
        'For public release, include the relevant study methods, data origin, class definition '
        'and applicable data permissions/license alongside these files. The export itself '
        'does not assign a new data license. If using project data/example/ inputs, all '
        'sequences and labels are synthetic and have no experimental meaning.\n'
    )
    files = {str(p.relative_to(out)): dict(sha256=sha256(p), bytes=p.stat().st_size)
             for p in sorted(out.rglob('*')) if p.is_file()}
    manifest = dict(schema_version=1, split_policy='existing frozen assignments; no new split',
                    validation_folds=5, shared_independent_test_set=True,
                    partitions=counts, files=files, provenance=provenance or {})
    write_json(out / 'dataset_manifest.json', manifest)
    with (out / 'SHA256SUMS').open('w') as stream:
        for path in sorted(p for p in out.rglob('*') if p.is_file() and p.name != 'SHA256SUMS'):
            stream.write(f'{sha256(path)}  {path.relative_to(out).as_posix()}\n')
    print(f'Exported fixed training/validation/test FASTA + CSV files to {out}', flush=True)
    return manifest


def export_existing(fasta, splits, out):
    """Export a prepared dataset after checking its frozen source/split hashes."""
    source = Path(splits).parent / 'manifest.json'
    manifest = json.loads(source.read_text())
    if sha256(fasta) != manifest['fasta_sha256'] or sha256(splits) != manifest['split_sha256']:
        raise ValueError('FASTA or frozen split hash differs from prepare manifest')
    sequences = read_fasta(fasta)
    table = load_samples(splits, sequences)
    # load_samples parses input IDs as strings; restore numerical validation folds.
    table['cv_fold'] = pd.to_numeric(table.cv_fold, errors='raise')
    return export_dataset(sequences, table, out, provenance=dict(
        fasta_sha256=sha256(fasta), split_sha256=sha256(splits),
        prepare_manifest_sha256=sha256(source), seed=manifest.get('seed'),
        identity=manifest.get('identity'), coverage=manifest.get('coverage'),
        grouping=manifest.get('grouping')))
