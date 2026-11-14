import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
import torch

from mutpromisloc.data import read_fasta, load_samples, context_window, local_tokens
from mutpromisloc.split import is_homolog, component_groups, assign_partitions
from mutpromisloc.model import Classifier, masked_pool
from mutpromisloc.training import fit_rms, select_threshold


class DataTests(unittest.TestCase):
    def test_multiline_fasta_and_duplicate_id_rejection(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / 'x.fa'
            p.write_text('>a description\nACD\nEFG\n>b\nVVVV\n')
            self.assertEqual(read_fasta(p), {'a': 'ACDEFG', 'b': 'VVVV'})
            p.write_text('>a\nACD\n>a\nVVV\n')
            with self.assertRaisesRegex(ValueError, 'Duplicate'):
                read_fasta(p)

    def test_samples_validate_pair_and_accept_unlabeled_prediction(self):
        seqs = {'wt': 'ACDEF', 'mut': 'ACVEF'}
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / 'samples.csv'
            p.write_text('sample_id,wt_id,mut_id,position,label\ns1,wt,mut,3,1\n')
            self.assertEqual(load_samples(p, seqs).label.tolist(), [1])
            p.write_text('sample_id,wt_id,mut_id,position\ns1,wt,mut,3\n')
            self.assertEqual(len(load_samples(p, seqs, require_labels=False)), 1)
            with self.assertRaisesRegex(ValueError, 'label'):
                load_samples(p, seqs)
            p.write_text('sample_id,wt_id,mut_id,position,label\ns1,wt,mut,2,1\n')
            with self.assertRaisesRegex(ValueError, 'position'):
                load_samples(p, seqs)

    def test_context_keeps_mutation_at_sequence_end(self):
        self.assertEqual(context_window('ACDEFGHIKL', 10, 6), ('FGHIKL', 5, 4))

    def test_local_tokens_exclude_special_tokens_and_pad_at_edges(self):
        h = np.array([[999], [1], [2], [3], [888]], dtype=np.float32)
        w, mask = local_tokens(h, 3, 0, 2)
        np.testing.assert_array_equal(w[:, 0], [0, 0, 1, 2, 3])
        np.testing.assert_array_equal(mask, [False, False, True, True, True])


class HomologyTests(unittest.TestCase):
    def test_thresholds_are_inclusive_and_longer_protein_is_used(self):
        h = dict(nident=64, alnlen=160, qlen=100, tlen=200, qstart=1, qend=70, tstart=1, tend=160)
        self.assertTrue(is_homolog(h))  # shorter coverage =70%, longer=80%
        self.assertFalse(is_homolog({**h, 'tend': 159}))
        self.assertFalse(is_homolog({**h, 'nident': 63}))
        swapped = dict(nident=64, alnlen=160, qlen=200, tlen=100, qstart=1, qend=160, tstart=1, tend=70)
        self.assertTrue(is_homolog(swapped))

    def test_equal_lengths_check_both_residue_spans(self):
        h = dict(nident=40, alnlen=100, qlen=100, tlen=100, qstart=1, qend=80, tstart=1, tend=79)
        self.assertFalse(is_homolog(h))

    def test_components_merge_transitive_edges_pairs_and_exact_duplicates(self):
        sequences = {'a': 'AAA', 'b': 'AAC', 'c': 'ACC', 'd': 'ACC', 'x': 'XXX'}
        samples = pd.DataFrame([{'wt_id': 'a', 'mut_id': 'b'}, {'wt_id': 'c', 'mut_id': 'd'}])
        groups = component_groups(sequences, samples, [('b', 'c')])
        self.assertEqual(groups['a'], groups['d'])
        self.assertNotEqual(groups['a'], groups['x'])

    def test_outer_inner_groups_disjoint_deterministic_and_class_balanced(self):
        meta = pd.DataFrame({'sample_id': [f's{i}' for i in range(120)],
                             'label': np.tile([0, 1], 60), 'cluster_id': [f'g{i // 2}' for i in range(120)]})
        s = assign_partitions(meta, seed=2026)
        pd.testing.assert_frame_equal(s, assign_partitions(meta, seed=2026))
        a = s[s.split == 'development']; b = s[s.split == 'test']
        self.assertFalse(set(a.cluster_id) & set(b.cluster_id))
        self.assertTrue(0.15 <= len(b) / len(s) <= 0.25)
        for f in range(5):
            self.assertFalse(set(a.loc[a.cv_fold == f, 'cluster_id']) & set(a.loc[a.cv_fold != f, 'cluster_id']))
        self.assertEqual(set(a.cv_fold), set(range(5)))


class ModelTests(unittest.TestCase):
    def test_masked_pool_and_all_masked_rejection(self):
        x = torch.tensor([[[1., 2.], [3., 4.], [999., 999.]]])
        mean, maximum = masked_pool(x, torch.tensor([[True, True, False]]))
        torch.testing.assert_close(mean, torch.tensor([[2., 3.]]))
        torch.testing.assert_close(maximum, torch.tensor([[3., 4.]]))
        with self.assertRaises(ValueError):
            masked_pool(x, torch.zeros(1, 3, dtype=torch.bool))

    def test_padding_invariance_finite_gradients_and_output_shape(self):
        torch.manual_seed(4)
        model = Classifier(input_dim=4, channels=8, kernel=3, heads=2, ff_dim=16, dropout=0).eval()
        x = torch.randn(2, 5, 4)
        mask = torch.ones(2, 5, dtype=torch.bool)
        first = model(x, mask)
        padded = torch.cat([x, torch.full((2, 3, 4), 1234.)], dim=1)
        padding_mask = torch.cat([mask, torch.zeros(2, 3, dtype=torch.bool)], dim=1)
        torch.testing.assert_close(first, model(padded, padding_mask), atol=1e-5, rtol=1e-5)
        self.assertEqual(tuple(first.shape), (2,))
        first.sum().backward()
        self.assertTrue(all(torch.isfinite(p.grad).all() for p in model.parameters()))

    def test_rms_ignores_padding_and_holdout(self):
        x = np.array([[[3., 4.], [100., 100.]], [[999., 999.], [999., 999.]]])
        mask = np.array([[True, False], [True, True]])
        self.assertAlmostEqual(fit_rms(x, mask, [0]), np.sqrt(12.5))
        self.assertTrue(np.isfinite(select_threshold(np.array([0, 1]), np.array([-1., 1.]))))


if __name__ == '__main__':
    unittest.main()
