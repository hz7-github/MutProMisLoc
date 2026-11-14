"""Regression: real MMseqs2 must export counted identities, including WT/MUT hits."""
import os
from pathlib import Path
import shutil
import tempfile
import unittest

import pandas as pd

from mutpromisloc.split import prepare

MMSEQS = os.environ.get('MMSEQS_TEST_BINARY') or shutil.which('mmseqs')


@unittest.skipUnless(MMSEQS, 'Set MMSEQS_TEST_BINARY or put mmseqs on PATH')
class MMseqsIntegrationTest(unittest.TestCase):
    def test_real_alignment_identity_counts_and_known_mutation_edges(self):
        example = Path(__file__).resolve().parents[1] / 'data/example'
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / 'prepared'
            prepare(example/'proteins.fa', example/'samples.csv', out,
                    mmseqs=MMSEQS, threads=1, seed=2026)
            hits = pd.read_csv(out/'alignments.tsv', sep='\t',
                               names='query target nident alnlen qlen tlen qstart qend tstart tend'.split())
            self_hits = hits[hits['query'] == hits['target']]
            self.assertGreater(len(self_hits), 0)
            self.assertTrue(self_hits.nident.eq(self_hits.qlen).all(),
                            'Self alignments must have actual identical-residue counts')
            meta = pd.read_csv(example/'samples.csv')
            groups = pd.read_csv(out/'sequence_groups.csv').set_index('sequence_id')
            observed = set(zip(hits['query'], hits['target']))
            for r in meta.itertuples(index=False):
                a,b = groups.loc[r.wt_id,'tool_id'], groups.loc[r.mut_id,'tool_id']
                self.assertIn((a,b), observed)
                h = hits[(hits['query']==a)&(hits['target']==b)].iloc[0]
                self.assertGreaterEqual(h.nident / h.alnlen, .95)
            edges = pd.read_csv(out/'homology_edges.csv')
            self.assertGreater(len(edges), 0)


if __name__ == '__main__':
    unittest.main()
