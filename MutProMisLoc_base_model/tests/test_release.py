import hashlib
import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from mutpromisloc.data import read_fasta, write_fasta, sha256, write_json
from mutpromisloc.release import export_dataset, export_existing


def fixture():
    sequences, rows = {}, []
    # Each independent cluster has one WT and two variants, one per label.
    for i in range(6):
        wt = 'ACDEFGHIKLMNPQRSTVWY' + 'A' * i
        ids = [f'w{i}', f'm{i}a', f'm{i}b']
        sequences[ids[0]] = wt
        sequences[ids[1]] = wt[:2] + 'V' + wt[3:]
        sequences[ids[2]] = wt[:2] + 'A' + wt[3:]
        for label in [0, 1]:
            rows.append(dict(sample_id=f's{i}_{label}', wt_id=ids[0], mut_id=ids[1 + label],
                             position=3, label=label, cluster_id=f'g{i}',
                             split='test' if i == 5 else 'development',
                             cv_fold=-1 if i == 5 else i))
    return sequences, pd.DataFrame(rows)


class ReleaseTests(unittest.TestCase):
    def test_export_contains_fixed_test_and_all_five_train_validation_pairs(self):
        sequences, splits = fixture()
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / 'dataset'
            export_dataset(sequences, splits, out)
            pd.testing.assert_frame_equal(pd.read_csv(out / 'splits.csv'), splits)
            self.assertEqual(read_fasta(out / 'proteins.fa'), sequences)
            fixed_test = set(splits.loc[splits.split == 'test', 'sample_id'])
            development = set(splits.loc[splits.split == 'development', 'sample_id'])
            self.assertEqual(set(pd.read_csv(out / 'test.csv').sample_id), fixed_test)
            self.assertEqual(len(read_fasta(out / 'test.fa')), 3)  # Shared WT exported once.
            validations = []
            for f in range(5):
                fold = out / 'folds' / f'fold_{f}'
                train, val = pd.read_csv(fold / 'train.csv'), pd.read_csv(fold / 'validation.csv')
                self.assertFalse(set(train.sample_id) & set(val.sample_id))
                self.assertFalse((set(train.sample_id) | set(val.sample_id)) & fixed_test)
                self.assertEqual(set(train.sample_id) | set(val.sample_id), development)
                self.assertFalse(set(train.cluster_id) & set(val.cluster_id))
                for name, table in [('train', train), ('validation', val)]:
                    subset = read_fasta(fold / f'{name}.fa')
                    self.assertEqual(set(subset), set(table.wt_id) | set(table.mut_id))
                    self.assertTrue(all(sequences[key] == seq for key, seq in subset.items()))
                validations += val.sample_id.tolist()
            self.assertEqual(set(validations), development)
            self.assertEqual(len(validations), len(development))
            manifest = json.loads((out / 'dataset_manifest.json').read_text())
            for name, entry in manifest['files'].items():
                self.assertEqual(sha256(out / name), entry['sha256'])
            self.assertEqual(manifest['partitions']['test']['samples'], 2)
            self.assertEqual(manifest['partitions']['development']['samples'], 10)

    def test_export_rejects_cluster_crossing_inner_folds_before_writing(self):
        sequences, splits = fixture()
        splits.loc[1, 'cv_fold'] = 1  # Same cluster now appears on both sides of fold 0.
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / 'dataset'
            with self.assertRaisesRegex(ValueError, 'cluster'):
                export_dataset(sequences, splits, out)
            self.assertFalse(out.exists())

    def test_export_rejects_same_sequence_under_different_ids_crossing_test(self):
        sequences, splits = fixture()
        sequences['w5'] = sequences['w0']
        sequences['m5a'] = sequences['m0a']
        sequences['m5b'] = sequences['m0b']
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(ValueError, 'sequence'):
                export_dataset(sequences, splits, Path(tmp) / 'dataset')

    def test_existing_export_checks_frozen_hashes_and_does_not_rewrite_split(self):
        sequences, splits = fixture()
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fasta, split = root / 'input.fa', root / 'splits.csv'
            write_fasta(fasta, sequences)
            splits.to_csv(split, index=False)
            write_json(root / 'manifest.json', dict(fasta_sha256=sha256(fasta), split_sha256=sha256(split)))
            original = split.read_bytes()
            export_existing(fasta, split, root / 'dataset')
            self.assertEqual(split.read_bytes(), original)
            with self.assertRaisesRegex(ValueError, 'not empty'):
                export_existing(fasta, split, root / 'dataset')
            fasta.write_text(fasta.read_text() + '\n')
            with self.assertRaisesRegex(ValueError, 'hash'):
                export_existing(fasta, split, root / 'dataset2')


if __name__ == '__main__':
    unittest.main()
