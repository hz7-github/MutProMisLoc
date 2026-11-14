"""Command-line entry points; all paths are explicit and portable."""
import argparse
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(prog='mutpromisloc', description='ESM-2 WT/MUT Conv1D + Transformer classification')
    commands = parser.add_subparsers(dest='command', required=True)
    p = commands.add_parser('validate', help='Validate paired FASTA + samples.csv')
    p.add_argument('--fasta', required=True, type=Path)
    p.add_argument('--samples', required=True, type=Path)
    p.add_argument('--unlabeled', action='store_true')
    p = commands.add_parser('prepare', help='MMseqs2 homology graph, 4:1 split and five development folds')
    p.add_argument('--fasta', required=True, type=Path)
    p.add_argument('--samples', required=True, type=Path)
    p.add_argument('--out', required=True, type=Path)
    p.add_argument('--mmseqs', default='mmseqs')
    p.add_argument('--threads', type=int, default=8)
    p.add_argument('--seed', type=int, default=2026)
    p = commands.add_parser('export', help='Export frozen train/validation/test FASTA + CSV files for readers')
    p.add_argument('--fasta', required=True, type=Path)
    p.add_argument('--splits', required=True, type=Path)
    p.add_argument('--out', required=True, type=Path)
    p = commands.add_parser('extract', help='Extract frozen ESM-2 WT-MUT local token features')
    p.add_argument('--fasta', required=True, type=Path)
    p.add_argument('--samples', required=True, type=Path)
    p.add_argument('--out', required=True, type=Path)
    p.add_argument('--checkpoint', required=True, type=Path)
    p.add_argument('--splits', type=Path)
    p.add_argument('--device', default='cpu')
    p.add_argument('--batch-size', type=int, default=2)
    p.add_argument('--radius', type=int, default=30)
    p.add_argument('--max-context', type=int, default=1022)
    p = commands.add_parser('train', help='Five-fold development CV and median-epoch final refit')
    p.add_argument('--features', required=True, type=Path)
    p.add_argument('--splits', required=True, type=Path)
    p.add_argument('--out', required=True, type=Path)
    p.add_argument('--device', default='cpu')
    p.add_argument('--seed', type=int, default=2026)
    p.add_argument('--lr', type=float, default=.0003)
    p.add_argument('--batch-size', type=int, default=64)
    p.add_argument('--max-epochs', type=int, default=60)
    p.add_argument('--patience', type=int, default=8)
    p.add_argument('--dropout', type=float, default=.2)
    p = commands.add_parser('evaluate', help='Score the frozen independent test partition once')
    p.add_argument('--features', required=True, type=Path)
    p.add_argument('--splits', required=True, type=Path)
    p.add_argument('--run', required=True, type=Path)
    p.add_argument('--device', default='cpu')
    p = commands.add_parser('predict', help='Predict new unlabeled WT/MUT pairs from extracted features')
    p.add_argument('--features', required=True, type=Path)
    p.add_argument('--checkpoint', required=True, type=Path)
    p.add_argument('--out', required=True, type=Path)
    p.add_argument('--device', default='cpu')
    args = vars(parser.parse_args())
    command = args.pop('command')
    if command == 'validate':
        from .data import read_fasta, load_samples
        seqs = read_fasta(args['fasta'])
        samples = load_samples(args['samples'], seqs, require_labels=not args['unlabeled'])
        print(f'Valid: {len(seqs)} protein sequences, {len(samples)} WT/MUT samples')
    elif command == 'prepare':
        from .split import prepare
        prepare(**args)
    elif command == 'extract':
        from .features import extract
        extract(**args)
    elif command == 'export':
        from .release import export_existing
        export_existing(**args)
    else:
        from . import training
        getattr(training, command)(**args)


if __name__ == '__main__':
    main()
