import unittest
import pandas as pd
from mutpromisloc.features import build_requests


class FeatureTests(unittest.TestCase):
    def test_paired_windows_share_bounds_and_deduplicate_repeated_wt(self):
        seqs = {'w': 'ACDEFGHIKL', 'm': 'ACDEFGHIKA', 'n': 'ACDEFGHIKV'}
        samples = pd.DataFrame([dict(sample_id='a', wt_id='w', mut_id='m', position=10),
                                dict(sample_id='b', wt_id='w', mut_id='n', position=10)])
        requests, refs, windows = build_requests(seqs, samples, max_context=6)
        self.assertEqual(len(requests), 3)
        self.assertEqual(windows.window_start_0based.tolist(), [4, 4])
        self.assertEqual(windows.center_0based.tolist(), [5, 5])
        self.assertTrue(all(len(seq) == 6 for seq in requests.values()))
        shared = [v for v in refs.values() if len(v) == 2]
        self.assertEqual(len(shared), 1)
        self.assertTrue(all(role == 'wt' for _, role, _ in shared[0]))


if __name__ == '__main__':
    unittest.main()
