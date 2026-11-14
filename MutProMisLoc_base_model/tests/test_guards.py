import tempfile
import unittest
from pathlib import Path
import numpy as np
import pandas as pd

from mutpromisloc.data import load_samples, new_output
from mutpromisloc.split import component_groups, assign_partitions
from mutpromisloc.training import evaluate, load_features


class GuardTests(unittest.TestCase):
    def test_exact_duplicates_merge_even_without_alignment_or_pair_link(self):
        seqs = {'a': 'ACD', 'am': 'AVD', 'b': 'ACD', 'bm': 'AGD'}
        meta = pd.DataFrame([dict(wt_id='a', mut_id='am'), dict(wt_id='b', mut_id='bm')])
        groups = component_groups(seqs, meta, [])
        self.assertEqual(groups['am'], groups['bm'])

    def test_conflicting_duplicate_pair_labels_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / 'samples.csv'
            p.write_text('sample_id,wt_id,mut_id,position,label\na,w,m,2,1\nb,w,m,2,0\n')
            with self.assertRaisesRegex(ValueError, 'conflicting'):
                load_samples(p, {'w': 'ACD', 'm': 'AVD'})

    def test_too_few_clusters_rejected(self):
        meta = pd.DataFrame(dict(sample_id=['a', 'b'], label=[0, 1], cluster_id=['x', 'y']))
        with self.assertRaisesRegex(ValueError, 'clusters'):
            assign_partitions(meta)

    def test_nonempty_outputs_and_repeated_test_evaluation_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)
            (p / 'test_metrics.json').write_text('{}')
            with self.assertRaisesRegex(ValueError, 'not empty'):
                new_output(p)
            with self.assertRaisesRegex(ValueError, 'already evaluated'):
                evaluate(p / 'absent.npz', p / 'absent.csv', p)

    def test_modified_feature_file_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)
            (p / 'manifest.json').write_text('{"features_sha256":"incorrect"}')
            np.savez(p / 'features.npz', sample_id=['a'], difference=np.zeros((1,3,2)), mask=np.ones((1,3),dtype=bool))
            with self.assertRaisesRegex(ValueError, 'hash changed'):
                load_features(p / 'features.npz')


if __name__ == '__main__':
    unittest.main()
