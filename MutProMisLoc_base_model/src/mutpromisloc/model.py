"""Conv1D + Transformer + masked mean/max binary classifier."""
import torch
from torch import nn


def masked_pool(x, mask):
    if x.ndim != 3 or mask.shape != x.shape[:2] or not mask.any(dim=1).all():
        raise ValueError('Expected [B,L,C], matching mask, and at least one valid position per row')
    valid = mask.unsqueeze(-1)
    mean = x.masked_fill(~valid, 0).sum(dim=1) / valid.sum(dim=1)
    maximum = x.masked_fill(~valid, torch.finfo(x.dtype).min).amax(dim=1)
    return mean, maximum


class Classifier(nn.Module):
    def __init__(self, input_dim=1280, channels=128, kernel=15, heads=8, ff_dim=512, dropout=0.2):
        super().__init__()
        if kernel % 2 != 1 or channels % heads:
            raise ValueError('Kernel must be odd and channels divisible by heads')
        self.conv = nn.Conv1d(input_dim, channels, kernel, padding=kernel // 2)
        layer = nn.TransformerEncoderLayer(channels, heads, ff_dim, dropout,
                                           activation='relu', batch_first=True)
        self.encoder = nn.TransformerEncoder(layer, num_layers=1, enable_nested_tensor=False)
        self.output = nn.Linear(channels * 2, 1)

    def forward(self, difference, mask):
        if difference.ndim != 3 or mask.shape != difference.shape[:2] or not mask.any(dim=1).all():
            raise ValueError('Input/mask shape mismatch or empty sequence')
        mask = mask.bool()
        x = difference.float().masked_fill(~mask.unsqueeze(-1), 0)
        x = torch.relu(self.conv(x.transpose(1, 2)).transpose(1, 2))
        x = x.masked_fill(~mask.unsqueeze(-1), 0)
        x = self.encoder(x, src_key_padding_mask=~mask)
        mean, maximum = masked_pool(x, mask)
        return self.output(torch.cat([mean, maximum], dim=-1)).squeeze(-1)
