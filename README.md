# MutProMisLoc: A Comprehensive Database of Subcellular Mislocalization of Human Mutated Proteins

A biomedical literature workflow for mutation-associated protein localization:
retrieve PubMed metadata, download public full text, convert PDFs with MinerU,
and extract candidate records with DeepSeek.

## Requirements

- Ubuntu with Bash and Conda/Miniconda; Python 3.12.
- Network access for literature retrieval, model downloads and API requests.
- A GPU environment with compatible drivers and inference dependencies for Step 03.
- Your own DeepSeek API key for new requests in Steps 04 or 05.

Official resources:

- [MinerU repository and installation](https://github.com/opendatalab/MinerU)
- [MinerU documentation](https://opendatalab.github.io/MinerU/)
- [DeepSeek API documentation](https://api-docs.deepseek.com/)
- [Conda getting started](https://docs.conda.io/projects/conda/en/stable/user-guide/getting-started.html)

## Installation

Clone or download this repository onto your Ubuntu server. Replace the example
path below. If the environment already exists, activate it instead of recreating it.

```bash
cd /path/to/MutProMisLoc
conda create -n MutProMisLoc python=3.12 pip -y
conda activate MutProMisLoc
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python -m pip install "mineru>=4.0,<5"
python -m pip check
"$CONDA_PREFIX/bin/mineru-kit" --help
```

MinerU is installed directly into Conda; do not additionally use Step 03 `--install`.
For NVIDIA GPUs, optionally install the Torch extra for the same MinerU version:

```bash
nvidia-smi
MINERU_VERSION=$(python -c "from importlib.metadata import version; print(version('mineru'))")
python -m pip install "mineru[torch]==${MINERU_VERSION}"
python -c "import torch; print('CUDA available:', torch.cuda.is_available())"
```

This does not install NVIDIA drivers. Follow the
[MinerU hardware and installation instructions](https://github.com/opendatalab/MinerU).
Confirm acceleration using runtime logs, CUDA availability and `nvidia-smi`.
The first conversion may download models.

## Repository contents

The repository includes six entry points, shared runtime modules, example configuration and Demo metadata.
PDFs, converted full text, model responses, database caches and results are generated when you run the pipeline.

| Step | Program | Main output |
|---|---|---|
| 01 | [01_get_pubmed_metadata.py](scripts/01_get_pubmed_metadata.py) | PubMed metadata XLSX/CSV |
| 02 | [02_download_papers.py](scripts/02_download_papers.py) | Public PDFs and download status |
| 03 | [03_mineru_pdf_to_md.py](scripts/03_mineru_pdf_to_md.py) | Raw and clean Markdown |
| 04 | [04_extract_deepseek.py](scripts/04_extract_deepseek.py) | `extracted_records.json`, per-paper records and model responses |
| 05 | [05_build_records.py](scripts/05_build_records.py) | Source-adjudicated `records.xlsx/.json/.csv` |
| 06 | [06_enrich_records.py](scripts/06_enrich_records.py) | Enriched `completed_records.xlsx/.json/.csv` |

Shared runtime modules live in `scripts/_internal/` and are loaded by the numbered entry points.
Demo metadata consists of `papers_available.xlsx`, `papers_sample100.xlsx` and `papers_sample100.csv`
in `demo/data/01_sample_papers/`.

## API configuration and reference data

Copy `.env.example` to `.env` at the project root and set `DEEPSEEK_API_KEY`.
Steps 04 and 05 use model requests; Step 05 requires `--online` to allow new requests.
Prepare the HGNC name reference JSON before Step 06:

```bash
mkdir -p data/06_completed_records/reference
curl --fail --location https://storage.googleapis.com/public-download-files/hgnc/json/json/hgnc_complete_set.json --output data/06_completed_records/reference/hgnc_complete_set.json
```

Alternatively, supply an existing reference with `--hgnc`. Use `--online` for new database queries;
subsequent runs can reuse locally stored caches.

## Run the main pipeline

The main pipeline defaults to `expanded`; no preset argument is needed at each step. Step 01 writes `papers_expanded.xlsx/.csv`, and subsequent steps use the matching files and directories.

```bash
python scripts/01_get_pubmed_metadata.py
python scripts/02_download_papers.py
python scripts/03_mineru_pdf_to_md.py --mineru-command "$CONDA_PREFIX/bin/mineru-kit" --model-source modelscope
python scripts/04_extract_deepseek.py
python scripts/05_build_records.py --online --output-dir data/05_records/run1
python scripts/06_enrich_records.py --input data/05_records/run1/records.json --output-dir data/06_completed_records/run1 --online
```

Step 01 defaults to a publication cutoff of 2025-08-01; use `--end-date` to change it.
Manually obtained PDFs can be named `<PMID>.pdf` in the matching download directory before Step 03.

## Run the Demo

Demo entry points explicitly use the sample metadata and step directories under `demo/data/`; the main preset default does not change the Demo input scope.

The 100-paper sample metadata and sampling input are included. Skip Step 01 when using the supplied sample;
add `--overwrite` only when deliberately replacing that sample table.

```bash
python demo/scripts/01_demo_sample_papers.py
python demo/scripts/02_demo_download_papers.py
python demo/scripts/03_demo_mineru_pdf_to_md.py --mineru-command "$CONDA_PREFIX/bin/mineru-kit" --model-source modelscope
python demo/scripts/04_demo_extract_deepseek.py
python demo/scripts/05_demo_build_records.py --fresh --online --output-dir demo/data/05_records/run1
python scripts/06_enrich_records.py --input demo/data/05_records/run1/records.json --output-dir demo/data/06_completed_records/run1 --online --include-unreviewed
```

Demo Step 05 must use `--fresh --online` to adjudicate the new extraction.
Saved adjudication snapshots are not bundled; the default snapshot mode requires a complete snapshot supplied separately.

| Step | Output path relative to the project root |
|---|---|
| 01 | `demo/data/01_sample_papers/papers_sample100.xlsx/.csv` |
| 02 | `demo/data/02_download_papers/<PMID>.pdf` and `download_status.csv` |
| 03 | `demo/data/03_mineru_pdf_to_md/raw/`, `clean/` and `conversion_status.csv` |
| 04 | `demo/data/04_extract_deepseek/extracted_records.json`, per-paper results and reports |
| 05 | `demo/data/05_records/run1/records.xlsx/.json/.csv` |
| 06 | `demo/data/06_completed_records/run1/completed_records.xlsx/.json/.csv` |

Final Excel files contain a single `Records` worksheet, one header row, exactly 26 columns and no cell fills.
Step 05 exports reviewed records; accessions, review status and reference sequences are enriched at the final Step 06.
Step 06 exports records with complete target fields and no detected conflicts; database localization is required only where applicable.
Statistics and audits are separate JSON/CSV files. Complete sequences remain in JSON and `wild_sequences.fasta`.
Use JSON between Steps 05 and 06. Counts depend on the actual download, extraction, review and enrichment results.

## Validation and resume

```bash
python scripts/01_get_pubmed_metadata.py --self-test
python scripts/05_build_records.py --help
python scripts/06_enrich_records.py --help
```

Use a fresh output directory for each result set, such as `run2`. Downloads, conversion and database queries retain status or caches.
When inputs or rules change, use a separate directory. See each entry point's `--help` for all available options.
