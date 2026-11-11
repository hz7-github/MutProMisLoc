"""Run the main PDF downloader on demo sample metadata; forward extra CLI options."""
from pathlib import Path
import runpy
import sys


def main():
    project = Path(__file__).resolve().parents[2]
    data = project / 'demo/data'
    command = [str(project / 'scripts/02_download_papers.py'),
               '--preset', 'expanded', '--input', str(data / '01_sample_papers/papers_sample100.xlsx'),
               '--output-dir', str(data / '02_download_papers'), *sys.argv[1:]]
    sys.argv = command
    runpy.run_path(command[0], run_name='__main__')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
