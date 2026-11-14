"""Grouped development CV, final refit and separately protected test evaluation."""
import json
import os
os.environ.setdefault('CUBLAS_WORKSPACE_CONFIG', ':4096:8')
import random
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from sklearn.metrics import (roc_auc_score, average_precision_score, f1_score, accuracy_score,
                             matthews_corrcoef, precision_score, recall_score,
                             balanced_accuracy_score, confusion_matrix)

from .data import sha256, write_json, new_output
from .model import Classifier


def fit_rms(x, mask, rows):
    total, count = 0., 0
    for i in rows:
        valid = np.asarray(x[i][mask[i]], dtype=np.float64)
        total += np.square(valid).sum()
        count += valid.size
    if count == 0:
        raise ValueError('No valid training tokens')
    return max(float(np.sqrt(total / count)), 1e-8)


def select_threshold(y, scores):
    grid = np.unique(np.r_[0., np.quantile(scores, np.linspace(.02, .98, 97))])
    return max((matthews_corrcoef(y, scores >= t), -abs(t), float(t)) for t in grid)[2]


def metrics(y, scores, threshold=0.):
    labels = np.asarray(scores) >= threshold
    tn, fp, fn, tp = confusion_matrix(y, labels, labels=[0, 1]).ravel()
    return dict(roc_auc=float(roc_auc_score(y, scores)), ap=float(average_precision_score(y, scores)),
                f1=float(f1_score(y, labels, zero_division=0)), accuracy=float(accuracy_score(y, labels)),
                balanced_accuracy=float(balanced_accuracy_score(y, labels)),
                mcc=float(matthews_corrcoef(y, labels)), precision=float(precision_score(y, labels, zero_division=0)),
                recall=float(recall_score(y, labels, zero_division=0)), specificity=float(tn / (tn + fp)),
                TN=int(tn), FP=int(fp), FN=int(fn), TP=int(tp))


def seed_all(seed):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.use_deterministic_algorithms(True)
    torch.backends.cudnn.benchmark = False
    torch.backends.cudnn.deterministic = True
    torch.backends.cuda.matmul.allow_tf32 = False
    torch.backends.cudnn.allow_tf32 = False
    torch.backends.cuda.enable_flash_sdp(False)
    torch.backends.cuda.enable_mem_efficient_sdp(False)
    torch.set_num_threads(2)


def load_features(features):
    features = Path(features)
    manifest = json.loads((features.parent / 'manifest.json').read_text())
    if sha256(features) != manifest['features_sha256']:
        raise ValueError('Feature hash changed')
    with np.load(features, allow_pickle=False) as z:
        ids = z['sample_id']
        x, mask = z['difference'].astype(np.float32), z['mask']
    if x.ndim != 3 or mask.shape != x.shape[:2] or not np.isfinite(x).all() or not mask.any(axis=1).all():
        raise ValueError('Invalid feature/mask arrays')
    return ids, x, mask, manifest


def load_development(features, splits):
    ids, x, mask, manifest = load_features(features)
    meta = pd.read_csv(splits, dtype={'sample_id': str})
    if meta.sample_id.tolist() != ids.tolist() or manifest['split_sha256'] != sha256(splits):
        raise ValueError('Feature rows or split hash differ')
    if set(meta.split) != {'development', 'test'}:
        raise ValueError('Expected development/test partitions')
    dev, test = meta.split.eq('development').to_numpy(), meta.split.eq('test').to_numpy()
    if set(meta.loc[dev, 'cluster_id']) & set(meta.loc[test, 'cluster_id']):
        raise ValueError('Development/test cluster overlap')
    for f in range(5):
        a, b = dev & meta.cv_fold.ne(f).to_numpy(), dev & meta.cv_fold.eq(f).to_numpy()
        if not b.any() or set(meta.loc[a, 'cluster_id']) & set(meta.loc[b, 'cluster_id']):
            raise ValueError('Missing CV fold or fold cluster overlap')
        if set(meta.loc[a, 'label']) != {0, 1} or set(meta.loc[b, 'label']) != {0, 1}:
            raise ValueError('Training/validation fold lacks a class')
    return meta, x, mask, manifest


def score_model(model, x, mask, rows, scale, batch_size):
    model.eval()
    output = []
    with torch.inference_mode():
        for start in range(0, len(rows), batch_size):
            idx = torch.as_tensor(rows[start:start + batch_size], device=x.device)
            output.append(model(x[idx] / scale, mask[idx]).cpu().numpy())
    return np.concatenate(output)


def fit_model(x, mask, y, rows, scale, config, seed, validation=None, epochs=None):
    seed_all(seed)
    model = Classifier(**config['model']).to(x.device)
    targets = torch.as_tensor(y, device=x.device, dtype=torch.float32)
    pos = int(y[rows].sum())
    if pos == 0 or pos == len(rows):
        raise ValueError('Training requires both labels')
    weight = (len(rows) - pos) / pos
    loss_fn = torch.nn.BCEWithLogitsLoss(pos_weight=torch.tensor(weight, device=x.device))
    optimizer = torch.optim.Adam(model.parameters(), lr=config['lr'], weight_decay=config['weight_decay'])
    rng = np.random.default_rng(seed)
    best_auc, best_ap, best_epoch, state, best_scores = -1., -1., 0, None, None
    history = []
    for epoch in range(1, (epochs or config['max_epochs']) + 1):
        model.train()
        ordered = rng.permutation(rows)
        loss_sum = 0.
        for start in range(0, len(rows), config['batch_size']):
            idx = torch.as_tensor(ordered[start:start + config['batch_size']], device=x.device)
            optimizer.zero_grad(set_to_none=True)
            loss = loss_fn(model(x[idx] / scale, mask[idx]), targets[idx])
            if not torch.isfinite(loss):
                raise ValueError('Nonfinite training loss')
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.)
            optimizer.step()
            loss_sum += loss.item() * len(idx)
        row = dict(epoch=epoch, training_loss=loss_sum / len(rows))
        if validation is not None:
            scores = score_model(model, x, mask, validation, scale, config['batch_size'])
            auc, ap = roc_auc_score(y[validation], scores), average_precision_score(y[validation], scores)
            row.update(validation_auc=float(auc), validation_ap=float(ap))
            if auc > best_auc + 1e-6 or (abs(auc - best_auc) <= 1e-6 and ap > best_ap):
                best_auc, best_ap, best_epoch = auc, ap, epoch
                state = {k: v.detach().cpu().clone() for k,v in model.state_dict().items()}
                best_scores = scores.copy()
        history.append(row)
        print(f'epoch={epoch} {row}', flush=True)
        if validation is not None and epoch - best_epoch >= config['patience']:
            break
    if validation is not None:
        model.load_state_dict(state)
    else:
        best_epoch = len(history)
    return model, best_epoch, best_scores, history, weight


def train(features, splits, out, device='cpu', seed=2026, lr=.0003, batch_size=64,
          max_epochs=60, patience=8, dropout=.2):
    if min(batch_size, max_epochs, patience) < 1 or lr <= 0 or not 0 <= dropout < 1:
        raise ValueError('Invalid training hyperparameters')
    meta, x, mask, manifest = load_development(features, splits)
    out = new_output(out)
    config = dict(seed=seed, lr=lr, batch_size=batch_size, max_epochs=max_epochs,
                  patience=patience, weight_decay=.01,
                  model=dict(input_dim=x.shape[2], channels=128, kernel=15, heads=8, ff_dim=512, dropout=dropout))
    protocol = dict(config=config, preprocessing=manifest['preprocessing'],
                    features=str(Path(features).resolve()), features_sha256=sha256(features),
                    splits=str(Path(splits).resolve()), split_sha256=sha256(splits),
                    early_stopping='validation ROC-AUC, AP tie-break; same fold used for scoring so CV is optimistic',
                    final_epochs='median of five best checkpoint epoch counts',
                    threshold='training OOF MCC maximized over fixed score quantiles',
                    test_used_for_model_selection=False, torch=str(torch.__version__))
    write_json(out / 'protocol.json', protocol)
    xt, mt = torch.as_tensor(x, device=device), torch.as_tensor(mask, device=device)
    development = np.flatnonzero(meta.split.eq('development').to_numpy())
    y = meta.label.to_numpy(dtype=int).copy()
    y[meta.split.eq('test').to_numpy()] = -999
    oof = np.full(len(meta), np.nan)
    fold_results, best_epochs = [], []
    for fold in range(5):
        a = development[meta.iloc[development].cv_fold.to_numpy() != fold]
        b = development[meta.iloc[development].cv_fold.to_numpy() == fold]
        scale = fit_rms(x, mask, a)
        model, best, scores, history, weight = fit_model(xt, mt, y, a, scale, config, seed + fold, validation=b)
        oof[b] = scores
        best_epochs.append(best)
        fold_results.append(dict(fold=fold, best_epoch=best, scale=scale, positive_weight=weight,
                                 train_samples=len(a), validation_samples=len(b), **metrics(y[b], scores)))
        pd.DataFrame(history).to_csv(out / f'history_fold{fold}.csv', index=False)
        torch.save(dict(state_dict={k: v.cpu() for k,v in model.state_dict().items()},
                        model_config=config['model'], preprocessing=manifest['preprocessing'], scale=scale),
                   out / f'fold{fold}.pt')
        pd.DataFrame(fold_results).to_csv(out / 'cv_metrics.csv', index=False)
    if not np.isfinite(oof[development]).all():
        raise ValueError('Incomplete OOF predictions')
    threshold = select_threshold(y[development], oof[development])
    epochs = max(1, int(np.rint(np.median(best_epochs))))
    scale = fit_rms(x, mask, development)
    selection = dict(cv_auc_mean=float(np.mean([r['roc_auc'] for r in fold_results])),
                     cv_auc_std=float(np.std([r['roc_auc'] for r in fold_results], ddof=1)),
                     final_epochs=epochs, best_epochs=best_epochs, threshold=threshold,
                     oof_metrics=metrics(y[development], oof[development], threshold),
                     test_used_for_model_selection=False)
    write_json(out / 'selection_before_refit.json', selection)
    frame = meta.iloc[development][['sample_id', 'label', 'cluster_id', 'cv_fold']].copy()
    frame['decision_score'] = oof[development]
    frame['predicted_label'] = (oof[development] >= threshold).astype(int)
    frame.to_csv(out / 'oof_predictions.csv', index=False)
    model, _, _, history, weight = fit_model(xt, mt, y, development, scale, config, seed, epochs=epochs)
    pd.DataFrame(history).to_csv(out / 'history_final.csv', index=False)
    torch.save(dict(state_dict={k: v.cpu() for k,v in model.state_dict().items()},
                    model_config=config['model'], scale=scale, threshold=threshold,
                    preprocessing=manifest['preprocessing'], features_sha256=sha256(features),
                    split_sha256=sha256(splits)), out / 'model.pt')
    selection.update(model_sha256=sha256(out / 'model.pt'), protocol_sha256=sha256(out / 'protocol.json'))
    write_json(out / 'selection.json', selection)
    print(selection, flush=True)


def checkpoint_model(path, device):
    seed_all(2026)
    bundle = torch.load(path, map_location='cpu', weights_only=True)
    model = Classifier(**bundle['model_config']).to(device)
    model.load_state_dict(bundle['state_dict'])
    model.eval()
    return model, bundle


def prediction_table(ids, scores, threshold):
    return pd.DataFrame(dict(sample_id=ids, decision_score=scores,
                             positive_score=torch.sigmoid(torch.as_tensor(scores)).numpy(),
                             predicted_label=(scores >= threshold).astype(int)))


def evaluate(features, splits, run, device='cpu'):
    run = Path(run)
    if (run / 'test_metrics.json').exists():
        raise ValueError('Test already evaluated; inspect saved results instead of tuning against test')
    selection = json.loads((run / 'selection.json').read_text())
    if sha256(run / 'model.pt') != selection['model_sha256'] or sha256(run / 'protocol.json') != selection['protocol_sha256']:
        raise ValueError('Selected model/protocol changed')
    meta, x, mask, manifest = load_development(features, splits)
    model, bundle = checkpoint_model(run / 'model.pt', device)
    if bundle['features_sha256'] != sha256(features) or bundle['split_sha256'] != sha256(splits):
        raise ValueError('Model/feature/split provenance mismatch')
    rows = np.flatnonzero(meta.split.eq('test').to_numpy())
    scores = score_model(model, torch.as_tensor(x, device=device), torch.as_tensor(mask, device=device), rows, bundle['scale'], 64)
    pred = prediction_table(meta.iloc[rows].sample_id.to_numpy(), scores, bundle['threshold'])
    pred['label'] = meta.iloc[rows].label.to_numpy()
    pred['cluster_id'] = meta.iloc[rows].cluster_id.to_numpy()
    pred.to_csv(run / 'test_predictions.csv', index=False)
    result = dict(test_samples=len(rows), threshold=bundle['threshold'],
                  metrics=metrics(pred.label.to_numpy(), scores, bundle['threshold']),
                  test_used_for_model_selection=False, train_test_cluster_overlap=0,
                  predictions_sha256=sha256(run / 'test_predictions.csv'))
    write_json(run / 'test_metrics.json', result)
    print(result, flush=True)


def predict(features, checkpoint, out, device='cpu'):
    ids, x, mask, manifest = load_features(features)
    model, bundle = checkpoint_model(checkpoint, device)
    if manifest['preprocessing'] != bundle['preprocessing']:
        raise ValueError('Prediction preprocessing/checkpoint differs from training; use the same ESM-2 weights, layer, radius and context length')
    out = Path(out)
    if out.exists():
        raise ValueError('Prediction output exists; choose a new filename')
    out.parent.mkdir(parents=True, exist_ok=True)
    scores = score_model(model, torch.as_tensor(x, device=device), torch.as_tensor(mask, device=device), np.arange(len(ids)), bundle['scale'], 64)
    prediction_table(ids, scores, bundle['threshold']).to_csv(out, index=False)
    print(f'Saved {out}', flush=True)
