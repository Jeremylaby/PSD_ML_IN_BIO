"""FASTA parsing with scikit-bio."""

from pathlib import Path

import skbio


def read_fasta(path: str | Path) -> list[tuple[str, str]]:
    """Read a multi-record FASTA file into a list of (id, sequence) pairs."""
    records = []
    for seq in skbio.io.read(str(path), format="fasta", constructor=skbio.Sequence):
        records.append((seq.metadata.get("id", ""), str(seq)))
    return records
