"""Run the main MinerU converter on demo PDFs; forward extra CLI options."""
from pathlib import Path
import runpy
import sys


def main():
    project = Path(__file__).resolve().parents[2]
    data = project / 'demo/data'
    command = [str(project / 'scripts/03_mineru_pdf_to_md.py'),
               '--preset', 'expanded', '--input-dir', str(data / '02_download_papers'),
               '--output-dir', str(data / '03_mineru_pdf_to_md'), *sys.argv[1:]]
    sys.argv = command
    runpy.run_path(command[0], run_name='__main__')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
