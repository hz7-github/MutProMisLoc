"""FASTA/sample validation and aligned mutation-window preprocessing."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

AMINO_ACIDS = set('ACDEFGHIKLMNPQRSTVWYBXZOU')


def sha256(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def write_json(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False) + '\n')


def new_output(path):
    path = Path(path)
    if path.exists() and any(path.iterdir()):
        raise ValueError(f'Output directory is not empty: {path}. Choose a new directory.')
    path.mkdir(parents=True, exist_ok=True)
    return path


def read_fasta(path):
    sequences, key, pieces = {}, None, []

    def finish():
        if key is not None:
            seq = ''.join(pieces).upper()
            if not seq or set(seq) - AMINO_ACIDS:
                raise ValueError(f'Invalid protein sequence for {key}: no gaps, stop symbols or empty sequences allowed')
            sequences[key] = seq

    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            finish()
            if not line[1:].strip():
                raise ValueError('Empty FASTA header')
            key = line[1:].split()[0]
            if key in sequences:
                raise ValueError(f'Duplicate FASTA ID: {key}')
            pieces = []
        else:
            if key is None:
                raise ValueError('Sequence appears before first FASTA header')
            pieces.append(line)
    finish()
    if not sequences:
        raise ValueError('FASTA file is empty')
    return sequences


def write_fasta(path, sequences):
    with Path(path).open('w') as stream:
        for name, seq in sequences.items():
            stream.write(f'>{name}\n')
            for i in range(0, len(seq), 80):
                stream.write(seq[i:i + 80] + '\n')


def load_samples(path, sequences, require_labels=True):
    df = pd.read_csv(path, dtype=str, keep_default_na=False)
    required = {'sample_id', 'wt_id', 'mut_id', 'position'}
    if require_labels:
        required.add('label')
    if not required <= set(df):
        raise ValueError(f'Missing CSV columns: {sorted(required - set(df))}')
    if df.empty or not df.sample_id.is_unique or df.sample_id.str.strip().eq('').any():
        raise ValueError('sample_id must be nonempty and unique; CSV cannot be empty')
    if 'label' in df:
        if not df.label.isin(['0', '1']).all():
            raise ValueError('label must be 0 or 1; omit the column for unlabeled prediction')
        df['label'] = df.label.astype(int)
    seen_pairs = {}
    positions = []
    for r in df.to_dict('records'):
        if r['wt_id'] not in sequences or r['mut_id'] not in sequences:
            raise ValueError(f"{r['sample_id']}: wt_id/mut_id missing from FASTA")
        wt, mut = sequences[r['wt_id']], sequences[r['mut_id']]
        try:
            position = int(r['position'])
        except ValueError as exc:
            raise ValueError('position must be a 1-based integer') from exc
        if len(wt) != len(mut):
            raise ValueError(f"{r['sample_id']}: only equal-length single substitutions are supported")
        if not 1 <= position <= len(wt):
            raise ValueError(f"{r['sample_id']}: position outside sequence")
        changed = [i for i, (a,b) in enumerate(zip(wt, mut)) if a != b]
        if changed != [position - 1]:
            raise ValueError(f"{r['sample_id']}: sequences must differ only at the specified position")
        pair = (wt, mut)
        if pair in seen_pairs and 'label' in r and seen_pairs[pair] != r['label']:
            raise ValueError('Identical WT/MUT pair has conflicting labels')
        seen_pairs[pair] = r.get('label')
        positions.append(position)
    df['position'] = positions
    return df


def context_window(sequence, position, max_length=1022):
    if max_length < 1 or max_length > 1022 or not 1 <= position <= len(sequence):
        raise ValueError('Invalid maximum context length or mutation position')
    p = position - 1
    start = max(0, min(p - max_length // 2, len(sequence) - max_length))
    return sequence[start:start + max_length], p - start, start


def local_tokens(hidden, length, center, radius=30):
    hidden = np.asarray(hidden, dtype=np.float32)
    if hidden.ndim != 2 or len(hidden) < length + 1 or not 0 <= center < length or radius < 0:
        raise ValueError('Invalid token representation/window')
    residues = hidden[1:length + 1]  # Exclude BOS, EOS and padding.
    if not np.isfinite(residues).all():
        raise ValueError('Nonfinite ESM-2 representation')
    width = 2 * radius + 1
    out = np.zeros((width, residues.shape[1]), dtype=np.float32)
    mask = np.zeros(width, dtype=bool)
    left, right = max(0, center - radius), min(length, center + radius + 1)
    dest = left - center + radius
    out[dest:dest + right - left] = residues[left:right]
    mask[dest:dest + right - left] = True
    return out, mask
