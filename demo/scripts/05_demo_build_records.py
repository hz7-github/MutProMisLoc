"""Export the saved Demo snapshot; --fresh processes the raw extraction."""
from pathlib import Path
import runpy
import sys


def main():
    project = Path(__file__).resolve().parents[2]
    data = project / 'demo/data'
    extra = sys.argv[1:]
    fresh = '--fresh' in extra
    extra = [a for a in extra if a != '--fresh']
    command = [str(project / 'scripts/05_build_records.py'),
        '--metadata', str(data / '01_sample_papers/papers_sample100.csv'),
        '--output-dir', str(data / '05_records')]
    if fresh:
        command += ['--input', str(data / '04_extract_deepseek/extracted_records.json'),
                    '--md-dir', str(data / '03_mineru_pdf_to_md/clean')]
    else:
        base = data / '06_source_adjudication/human254_v1'
        command += ['--audit', str(base / 'export_A_event_dedup_audit.json'),
                    '--pipeline-root', str(base), '--workset', str(base / 'input/workset.json'),
                    '--baseline', str(data / '07_admit_records/final_records.json')]
    if '--help' in extra or '-h' in extra:
        print('Demo: saved snapshot by default; --fresh starts from raw extraction. --fresh --online permits new requests.')
    sys.argv = command + extra
    runpy.run_path(sys.argv[0], run_name='__main__')


if __name__ == '__main__':
    main()
