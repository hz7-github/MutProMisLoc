"""Sample paper metadata deterministically without downloading or extracting."""
from __future__ import annotations

import argparse
import csv
import random
from pathlib import Path

from openpyxl import Workbook, load_workbook
from openpyxl.styles import Font, PatternFill
from openpyxl.worksheet.table import Table, TableStyleInfo

DATA = Path(__file__).resolve().parents[1] / 'data/01_sample_papers'


def read_papers(path):
    if path.suffix.lower() == '.csv':
        with path.open(encoding='utf-8-sig', newline='') as stream:
            reader = csv.DictReader(stream)
            headers = list(reader.fieldnames or [])
            rows = list(reader)
    elif path.suffix.lower() == '.xlsx':
        book = load_workbook(path, read_only=True, data_only=True)
        try:
            if 'Papers' not in book.sheetnames:
                raise ValueError('Input workbook must contain a Papers sheet.')
            values = book['Papers'].iter_rows(values_only=True)
            headers = list(next(values))
            rows = [dict(zip(headers, row)) for row in values if any(v is not None for v in row)]
        finally:
            book.close()
    else:
        raise ValueError('Input must be XLSX or CSV.')
    if 'PMID' not in headers or len(set(headers)) != len(headers):
        raise ValueError('Unique headers including PMID are required.')
    seen = set()
    for row in rows:
        pmid = str(row.get('PMID') or '').strip()
        if not pmid.isascii() or not pmid.isdigit():
            raise ValueError('Invalid PMID: ' + repr(pmid))
        if pmid in seen:
            raise ValueError('Duplicate PMID: ' + pmid)
        seen.add(pmid)
        row['PMID'] = pmid
    return [h for h in headers if h != 'Sample_Order'], rows


def select_papers(rows, count=100, seed=2026):
    if not 1 <= count <= len(rows):
        raise ValueError('Count must be between 1 and the population size.')
    population = sorted(rows, key=lambda row: int(row['PMID']))
    return [dict(row, Sample_Order=i) for i, row in
            enumerate(random.Random(seed).sample(population, count), 1)]


def write_outputs(base, headers, rows, *, seed, population_size, source):
    headers = ['Sample_Order'] + [h for h in headers if h != 'Sample_Order']
    base.parent.mkdir(parents=True, exist_ok=True)
    book = Workbook()
    sheet = book.active
    sheet.title = 'Papers'
    sheet.append(headers)
    for row in rows:
        sheet.append([row.get(h) for h in headers])
    for row in sheet:
        for cell in row:
            if isinstance(cell.value, str):
                cell.data_type = 's'
    for cell in sheet[1]:
        cell.font = Font(bold=True, color='FFFFFF')
        cell.fill = PatternFill('solid', fgColor='24445F')
    sheet.freeze_panes = 'C2'
    for column in sheet.columns:
        name = column[0].value
        sheet.column_dimensions[column[0].column_letter].width = 65 if name in ('Title', 'Abstract') else 23
    table = Table(displayName='SamplePapers', ref=sheet.dimensions)
    table.tableStyleInfo = TableStyleInfo(name='TableStyleMedium2', showRowStripes=True)
    sheet.add_table(table)
    details = book.create_sheet('Sampling Details')
    for row in [('Seed', seed), ('Count', len(rows)), ('Population', population_size),
                ('Source', str(source)), ('Method', 'Numeric PMID sort, then random.Random(seed).sample'),
                ('Scope', 'Metadata selection only; no PDF or model processing')]:
        details.append(row)
    details.column_dimensions['A'].width = 24
    details.column_dimensions['B'].width = 85
    book.save(base.with_suffix('.xlsx'))
    with base.with_suffix('.csv').open('w', encoding='utf-8-sig', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=headers, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=DATA / 'papers_available.xlsx')
    parser.add_argument('--output-name', default='papers_sample100')
    parser.add_argument('--seed', type=int, default=2026)
    parser.add_argument('--count', type=int, default=100)
    parser.add_argument('--overwrite', action='store_true', help='Replace existing sample metadata only')
    args = parser.parse_args()
    if Path(args.output_name).name != args.output_name or any(c in args.output_name for c in '/\\:.') or not args.output_name.strip():
        parser.error('--output-name must be a plain filename stem.')
    base = DATA / args.output_name
    if args.input.resolve() in {base.with_suffix(s).resolve() for s in ('.xlsx', '.csv')}:
        parser.error('Output must not replace the population input.')
    if not args.overwrite and any(base.with_suffix(s).exists() for s in ('.xlsx', '.csv')):
        parser.error('Sample files already exist. Use them as supplied, or add --overwrite.')
    try:
        headers, population = read_papers(args.input)
        sample = select_papers(population, args.count, args.seed)
        write_outputs(base, headers, sample, seed=args.seed, population_size=len(population), source=args.input)
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    print(f'Sampled {len(sample)} of {len(population)} papers with seed {args.seed}.')
    print(base.with_suffix('.xlsx'))
    print(base.with_suffix('.csv'))


if __name__ == '__main__':
    main()
