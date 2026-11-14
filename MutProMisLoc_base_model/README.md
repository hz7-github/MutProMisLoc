# MutProMisLoc_base_model

A protein binary classification baseline using ESM-2 WT−MUT difference features with Conv1D + Transformer. Its environment, code, data and outputs are independent of the main project's literature processing workflow.

## Code and data

[Download code and data (repository ZIP)](https://github.com/hz7-github/MutProMisLoc/archive/refs/heads/main.zip), extract the archive and enter `MutProMisLoc_base_model/`.

- [src/mutpromisloc/](src/mutpromisloc/): data validation, homology-based splitting, feature extraction, training and prediction.
- [data/dataset/](data/dataset/): 3,995 samples, with 3,233 development samples and 762 independent test samples. Five-fold training/validation files are in `folds/fold_0/` through `folds/fold_4/`.
- `proteins.fa` / `samples.csv` contain the complete inputs; `splits.csv` defines the fixed partitions. Each partition includes matching FASTA and CSV files. CSV fields are `sample_id,wt_id,mut_id,position,label`; positions are 1-based, only single amino acid substitutions are supported, and labels retain the source dataset's 0/1 definitions.
- Homology criteria: sequence identity ≥40% and coverage of the longer protein ≥80%. Connected homology clusters remain within a partition. Use the supplied splits without rerunning `prepare`; the test set is excluded from training and tuning.
- `data/example/` contains synthetic test data; `tests/` contains code tests. Pretrained weights, trained models and extracted features must be obtained or generated separately.

## Installation and usage

Requires Linux, Conda and Python 3.12. Run all commands from this directory.

```bash
conda env create -f environment.yml
conda activate MutProMisLoc_base_model
python -m pip install -e . --no-deps
python -m mutpromisloc --help
```

Place the ESM-2 checkpoint `esm2_t33_650M_UR50D.pt` and its matching `esm2_t33_650M_UR50D-contact-regression.pt` in the same directory, then set the checkpoint path:

```bash
export ESM2_CHECKPOINT=/path/to/esm2_t33_650M_UR50D.pt
python -m mutpromisloc validate --fasta data/dataset/proteins.fa --samples data/dataset/samples.csv
python -m mutpromisloc extract \
  --fasta data/dataset/proteins.fa --samples data/dataset/samples.csv \
  --splits data/dataset/splits.csv --checkpoint "$ESM2_CHECKPOINT" \
  --out features/base_model --device cuda:0 --batch-size 2
python -m mutpromisloc train \
  --features features/base_model/features.npz --splits data/dataset/splits.csv \
  --out runs/base_model --device cuda:0
python -m mutpromisloc evaluate \
  --features features/base_model/features.npz --splits data/dataset/splits.csv \
  --run runs/base_model --device cuda:0
```

Outputs are saved under `features/` and `runs/`; existing outputs are not overwritten. Without a compatible GPU, replace `--device cuda:0` with `--device cpu`. To predict new samples, use a paired CSV without the `label` column, run `extract` without `--splits`, then run `predict`. See each command's `--help` for arguments.

Check data integrity by running `sha256sum -c SHA256SUMS` in `data/dataset/`. Run code tests from this directory with `python -m unittest discover -s tests -v`.
