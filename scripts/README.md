# Command guide

Run commands from the `MutProMisLoc/` project root. The numbered entry points execute the six pipeline steps.
Shared runtime modules are in `scripts/_internal/`; they are not separate command-line steps.

| Step | Entry point | Main output |
|---|---|---|
| 01 | [01_get_pubmed_metadata.py](01_get_pubmed_metadata.py) | PubMed metadata XLSX/CSV |
| 02 | [02_download_papers.py](02_download_papers.py) | Public PDFs and download status |
| 03 | [03_mineru_pdf_to_md.py](03_mineru_pdf_to_md.py) | Raw and clean Markdown with source information |
| 04 | [04_extract_deepseek.py](04_extract_deepseek.py) | Candidate JSON, per-paper results, responses and reports |
| 05 | [05_build_records.py](05_build_records.py) | Reviewed `records.xlsx/.json/.csv` |
| 06 | [06_enrich_records.py](06_enrich_records.py) | Enriched `completed_records.xlsx/.json/.csv` |

The main pipeline defaults to `expanded`. Step 01 writes `papers_expanded.xlsx/.csv`; later steps use the matching paths. Demo entry points keep their explicit `demo/data/` inputs and outputs.

## Step 05: Review and export

```bash
python scripts/05_build_records.py --online --output-dir data/05_records/run1
```

Default input: `data/04_extract_deepseek/expanded/extracted_records.json`, with matching metadata and Markdown.
Override paths with `--input`, `--metadata` and `--md-dir`. Default output is `data/05_records/<preset>/`.

`--online` enables new model review requests and requires `DEEPSEEK_API_KEY`.
`--max-requests` limits newly selected batches across initial review and follow-up rounds; it is not a precise cost cap.
`--repair-rounds` sets the maximum follow-up rounds, defaulting to one.
Without new requests or compatible saved responses, records remain pending and the exported table may be empty.
Use `--response-dir` to reuse a compatible initial response cache.

Saved review rounds can be supplied explicitly:

```bash
python scripts/05_build_records.py --workset INPUT_WORKSET.json --base FIRST_REVIEWS.json --followup FOLLOWUP.json --metadata TITLES.csv --output-dir data/05_records/reviewed_run
```

Repeat `--followup` for multiple rounds. Optional `--event-decisions` and `--context-flags` supply source-supported
event links and review holds. Identical gene/variant labels alone do not establish an identical experimental event.

| Output | Contents |
|---|---|
| `records.xlsx` | One `Records` worksheet, one header row, 26 columns, no cell fills |
| `records.json` | Complete exported rows and source provenance; recommended input for Step 06 |
| `records.csv` | Flat exported data |
| `export_manifest.json` | Counts, input/output hashes and export information |
| `audit/` | Review decisions, pending records and linked duplicates |
| `work/` | Intermediate inputs and results for a full processing run |

## Step 06: Enrich and export

```bash
python scripts/06_enrich_records.py --input data/05_records/run1/records.json --output-dir data/06_completed_records/run1 --online
```

Prepare the HGNC reference as described in the [project README](../README.md).
The default HGNC reference is `data/06_completed_records/reference/hgnc_complete_set.json`.
`--hgnc`, `--uniprot-cache` and `--reference-cache` override reference paths.
`--online` permits public database queries; existing caches can be reused without new queries.
`--max-requests` covers HTTP attempts, including retries, across enrichment and identity recovery.
`--include-unreviewed` permits unreviewed entries under the existing identity-matching rules.

Completion requires a valid `Accession_new`, a `Reviewed` status of reviewed/unreviewed, and a complete
`Wild Sequence_new` without detected reference-residue or position conflicts. `Subcellular_Location_Values`
is required only for applicable, explicit No change records. Other optional experiment fields may remain blank.
Records with detected identity conflicts or conflicting existing database values remain outside the final table.
Complex variant notation retains its not-checked status.

Use Step 05 JSON when possible. XLSX inputs use the `Records` data worksheet; complete values and provenance
are preserved more reliably in JSON. An unlabelled CSV is not a reviewed-input format.
When original sources are available separately, use `--identity-input`, `--original-records`,
or `--source-audit` with `--dedup-map`.

| Output | Contents |
|---|---|
| `completed_records.xlsx` | One `Records` worksheet, 26 columns, no cell fills |
| `completed_records.json/.csv` | Complete exported rows |
| `wild_sequences.fasta` | Full reference sequences; oversized Excel cells point to FASTA/JSON |
| `run_report.json` | Input/export counts, gaps, validation information and hashes |
| `audit/completeness.json/.csv` | Export decisions and missing fields for each input record |
| `audit/all_records.json` | All processed records, including incomplete ones |
| `work/` | Database enrichment, identity recovery and the input-bound workset |

## Demo

Only the sampling input and 100-paper metadata sample are bundled. PDFs, full text, model responses,
saved review snapshots and database caches must be generated or supplied separately.
The Demo entry points 01–05 are in `demo/scripts/`; Step 06 uses the main entry point.

After downloading, converting and extracting the sample:

```bash
python demo/scripts/05_demo_build_records.py --fresh --online --output-dir demo/data/05_records/run1
python scripts/06_enrich_records.py --input demo/data/05_records/run1/records.json --output-dir demo/data/06_completed_records/run1 --online --include-unreviewed
```

The default Demo snapshot mode needs a complete saved review snapshot. Use `--fresh --online` for a new run.
See the [project README](../README.md) for installation and all six Demo commands.

## Checks and reruns

```bash
python scripts/01_get_pubmed_metadata.py --self-test
python scripts/05_build_records.py --help
python scripts/06_enrich_records.py --help
```

Choose a new output directory for each result set. Downloads, conversion and database queries retain status or caches.
Use separate directories when inputs or policies change. Inspect each script's `--help` for all parameters.
