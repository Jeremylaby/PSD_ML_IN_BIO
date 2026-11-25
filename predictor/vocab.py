"""Amino-acid vocabulary and one-hot encoding of protein sequences."""

import torch

STANDARD_AA = "ACDEFGHIKLMNPQRSTVWY"
UNK = "X"
PAD = "<PAD>"

# Słownik: 20 aminokwasów + X (nieznany) + PAD (wypełnienie) = 22 "litery".
VOCAB: list[str] = [*STANDARD_AA, UNK, PAD]
TOKEN_TO_INDEX: dict[str, int] = {token: i for i, token in enumerate(VOCAB)}
VOCAB_SIZE = len(VOCAB)  # N = 22
UNK_INDEX = TOKEN_TO_INDEX[UNK]
PAD_INDEX = TOKEN_TO_INDEX[PAD]

THREE_LETTER: dict[str, str] = {
    "A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE",
    "G": "GLY", "H": "HIS", "I": "ILE", "K": "LYS", "L": "LEU",
    "M": "MET", "N": "ASN", "P": "PRO", "Q": "GLN", "R": "ARG",
    "S": "SER", "T": "THR", "V": "VAL", "W": "TRP", "Y": "TYR",
    "X": "UNK",
}


def normalize(sequence: str) -> str:
    """Upper-case the sequence and map every non-standard residue to X."""
    return "".join(c if c in STANDARD_AA else UNK for c in sequence.upper())


def encode_batch(
    sequences: list[str], max_len: int = 1024
) -> tuple[torch.Tensor, torch.Tensor, list[int]]:
    """One-hot encode sequences, padding or truncating them to ``max_len``.

    Returns:
        x: float tensor of shape (B, D, N)
        mask: bool tensor of shape (B, D), True for real residues
        lengths: L_i = min(len(sequence_i), D) for each sequence
    """
    batch_size = len(sequences)
    # Zaczynamy od tablicy wypełnionej samym PAD, potem wpisujemy sekwencje.
    indices = torch.full((batch_size, max_len), PAD_INDEX, dtype=torch.long)
    lengths: list[int] = []
    for i, seq in enumerate(sequences):
        seq = normalize(seq)[:max_len]  # za długie sekwencje są obcinane do max_len
        lengths.append(len(seq))
        if seq:
            indices[i, : len(seq)] = torch.tensor(
                [TOKEN_TO_INDEX[c] for c in seq], dtype=torch.long
            )
    # Maska: True = prawdziwe residuum, False = wypełnienie (PAD).
    mask = indices != PAD_INDEX
    # One-hot: każde residuum to wektor 22 zer z jedną jedynką na pozycji swojej litery.
    x = torch.nn.functional.one_hot(indices, num_classes=VOCAB_SIZE).float()
    return x, mask, lengths
