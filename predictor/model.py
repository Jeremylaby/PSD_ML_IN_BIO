"""Dummy structure model that outputs random CA coordinates."""

import torch
from torch import nn

from predictor.vocab import VOCAB_SIZE


class DummyStructureModel(nn.Module):
    """Ignores sequence content and predicts random CA coordinates.

    Input:  x of shape (B, D, N), mask of shape (B, D)
    Output: coordinates of shape (B, D, 3); padded positions are zero.
    """

    def __init__(self, vocab_size: int = VOCAB_SIZE, scale: float = 10.0):
        super().__init__()
        self.vocab_size = vocab_size
        self.scale = scale

    def forward(self, x: torch.Tensor, mask: torch.Tensor) -> torch.Tensor:
        if x.ndim != 3 or x.shape[-1] != self.vocab_size:
            raise ValueError(
                f"expected input of shape (B, D, {self.vocab_size}), got {tuple(x.shape)}"
            )
        if mask.shape != x.shape[:2]:
            raise ValueError(f"mask shape {tuple(mask.shape)} does not match {tuple(x.shape[:2])}")
        batch_size, max_len, _ = x.shape
        # Model "na niby": losuje współrzędne zamiast je przewidywać.
        coords = torch.randn(batch_size, max_len, 3, device=x.device) * self.scale
        return coords * mask.unsqueeze(-1)  # zerujemy pozycje wypełnienia

    @torch.no_grad()
    def predict(self, x: torch.Tensor, mask: torch.Tensor) -> list[torch.Tensor]:
        """Return one (L_i, 3) tensor per sequence, dropping padded positions."""
        coords = self(x, mask)
        return [c[m] for c, m in zip(coords, mask)]  # zostawiamy tylko prawdziwe residua
