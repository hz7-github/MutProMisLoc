import tempfile
import unittest
from pathlib import Path

import numpy as np
import torch
from mutpromisloc.model import Classifier
from mutpromisloc.training import checkpoint_model, score_model


class DeviceConsistencyTests(unittest.TestCase):
    @unittest.skipUnless(torch.cuda.is_available(), 'CUDA optional')
    def test_checkpoint_cpu_gpu_predictions_close_without_tf32_rounding(self):
        torch.manual_seed(12)
        torch.set_num_threads(2)
        model_config = dict(input_dim=1280, channels=128, kernel=15, heads=8, ff_dim=512, dropout=0.)
        model = Classifier(**model_config)
        x = torch.randn(4, 61, 1280)
        mask = torch.ones(4, 61, dtype=torch.bool)
        # Simulate a fresh runtime's default cuDNN TF32 setting.
        torch.backends.cudnn.allow_tf32 = True
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'model.pt'
            torch.save(dict(model_config=model_config, state_dict=model.state_dict()), path)
            cpu, _ = checkpoint_model(path, 'cpu')
            expected = score_model(cpu, x, mask, np.arange(4), 1., 64)
            gpu, _ = checkpoint_model(path, 'cuda:0')
            actual = score_model(gpu, x.cuda(0), mask.cuda(0), np.arange(4), 1., 64)
        np.testing.assert_allclose(actual, expected, atol=2e-5, rtol=2e-5)


if __name__ == '__main__':
    unittest.main()
