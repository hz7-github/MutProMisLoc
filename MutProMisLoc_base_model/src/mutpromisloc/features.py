"""Frozen ESM-2 inference shared by training and new-pair prediction."""
import hashlib
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

from .data import (read_fasta, load_samples, context_window, local_tokens,
                   sha256, write_json, new_output)


def build_requests(sequences, samples, max_context=1022):
    requests, references, windows = {}, defaultdict(list), []
    for i, row in enumerate(samples.itertuples(index=False)):
        wt = sequences[row.wt_id]
        mut = sequences[row.mut_id]
        wt_context, center, start = context_window(wt, row.position, max_context)
        mut_context = mut[start:start + len(wt_context)]
        for role, sequence in [('wt', wt_context), ('mut', mut_context)]:
            key = hashlib.sha256(sequence.encode()).hexdigest()
            requests[key] = sequence
            references[key].append((i, role, center))
        windows.append(dict(sample_id=row.sample_id, window_start_0based=start,
                            window_length=len(wt_context), center_0based=center,
                            full_length=len(wt)))
    return requests, references, pd.DataFrame(windows)


def extract(fasta, samples, out, checkpoint, device='cpu', batch_size=2,
            radius=30, max_context=1022, splits=None):
    import torch
    import esm
    if radius < 0 or batch_size < 1:
        raise ValueError('radius must be >=0 and batch size >=1')
    seqs = read_fasta(fasta)
    meta = load_samples(samples, seqs, require_labels=False)
    out = new_output(out)
    split_manifest = None
    if splits:
        split = pd.read_csv(splits, dtype={'sample_id': str, 'wt_id': str, 'mut_id': str})
        if split.sample_id.tolist() != meta.sample_id.tolist():
            raise ValueError('Splits and CSV sample order differ')
        for col in ['wt_id', 'mut_id', 'position', 'label']:
            if col not in meta or not np.array_equal(split[col], meta[col]):
                raise ValueError(f'Split/CSV mismatch: {col}')
        split_manifest = Path(splits).parent / 'manifest.json'
        provenance = __import__('json').loads(split_manifest.read_text())
        if sha256(fasta) != provenance['fasta_sha256'] or sha256(samples) != provenance['samples_sha256'] or sha256(splits) != provenance['split_sha256']:
            raise ValueError('Input hashes differ from prepared splits')
    requests, refs, windows = build_requests(seqs, meta, max_context)
    n, width, dim = len(meta), 2 * radius + 1, 1280
    # Accumulate in FP32, store in FP16. WT-MUT direction follows the reference.
    delta = np.zeros((n, width, dim), dtype=np.float32)
    masks = np.zeros((n, width), dtype=bool)
    seen = np.zeros((n, 2), dtype=bool)
    torch.set_num_threads(2)
    model, alphabet = esm.pretrained.load_model_and_alphabet_local(str(Path(checkpoint).resolve()))
    if model.num_layers != 33 or model.embed_dim != 1280:
        raise ValueError('Use esm2_t33_650M_UR50D, 33 layers and 1280 dimensions')
    model.eval().to(device)
    convert = alphabet.get_batch_converter()
    keys = sorted(requests, key=lambda k: (len(requests[k]), k))
    started = time.monotonic()
    for start in range(0, len(keys), batch_size):
        batch = keys[start:start + batch_size]
        _, _, tokens = convert([(k, requests[k]) for k in batch])
        with torch.inference_mode():
            hidden = model(tokens.to(device), repr_layers=[33], return_contacts=False)['representations'][33]
        hidden = hidden.float().cpu().numpy()
        for h, key in zip(hidden, batch):
            for i, role, center in refs[key]:
                token, mask = local_tokens(h, len(requests[key]), center, radius)
                if seen[i].any() and not np.array_equal(masks[i], mask):
                    raise AssertionError('WT/MUT valid masks differ')
                delta[i] += token if role == 'wt' else -token
                masks[i] = mask
                seen[i, 0 if role == 'wt' else 1] = True
        if start == 0 or start % 100 == 0 or start + batch_size >= len(keys):
            print(f'ESM-2 {min(start + batch_size,len(keys))}/{len(keys)} unique context windows', flush=True)
    if not seen.all() or not np.isfinite(delta).all():
        raise ValueError('Incomplete or nonfinite features')
    np.savez_compressed(out / 'features.npz', sample_id=meta.sample_id.to_numpy(dtype=str),
                        difference=delta.astype(np.float16), mask=masks)
    windows.to_csv(out / 'windows.csv', index=False)
    preprocessing = dict(model='esm2_t33_650M_UR50D', layer=33, dimension=1280,
                         radius=radius, max_context=max_context, difference='WT-MUT',
                         checkpoint_sha256=sha256(checkpoint))
    manifest = dict(fasta=str(Path(fasta).resolve()), samples=str(Path(samples).resolve()),
                    fasta_sha256=sha256(fasta), samples_sha256=sha256(samples),
                    checkpoint=str(Path(checkpoint).resolve()), preprocessing=preprocessing,
                    samples_count=n, unique_context_windows=len(keys), shape=list(delta.shape),
                    features_sha256=sha256(out / 'features.npz'), wall_seconds=time.monotonic() - started,
                    split_path=str(Path(splits).resolve()) if splits else None,
                    split_sha256=sha256(splits) if splits else None,
                    split_manifest_sha256=sha256(split_manifest) if split_manifest else None,
                    torch=torch.__version__, fair_esm=esm.__version__, device=str(device))
    write_json(out / 'manifest.json', manifest)
    print(f'Saved {out / "features.npz"}', flush=True)
