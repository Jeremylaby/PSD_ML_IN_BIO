"""Write CA-only PDB files following the wwPDB v3.3 ATOM record format."""

from pathlib import Path
from typing import TextIO

import torch

from predictor.vocab import THREE_LETTER, normalize


def format_ca_atom(serial: int, res_name: str, chain_id: str, res_seq: int,
                   x: float, y: float, z: float) -> str:
    """Format one 80-column ATOM record for a CA atom.

    Columns (1-based): 1-6 record name, 7-11 serial, 13-16 atom name,
    17 altLoc, 18-20 resName, 22 chainID, 23-26 resSeq, 27 iCode,
    31-38 x, 39-46 y, 47-54 z, 55-60 occupancy, 61-66 tempFactor,
    77-78 element, 79-80 charge.
    """
    # Format PDB ma sztywne kolumny: każde pole musi zaczynać się w ściśle określonym miejscu.
    return (
        f"{'ATOM':<6}{serial % 100000:>5} {' CA ':<4} {res_name:>3} {chain_id:1}"
        f"{res_seq % 10000:>4}    {x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{0.0:6.2f}"
        f"{'':10}{'C':>2}{'':2}"
    )


def write_model(handle: TextIO, model_id: int, sequence: str, coords: torch.Tensor,
                chain_id: str = "A") -> None:
    """Write one MODEL block with one CA record per processed residue."""
    sequence = normalize(sequence)[: len(coords)]
    # Każde białko z pliku FASTA trafia do osobnego bloku MODEL ... ENDMDL.
    handle.write(f"MODEL     {model_id:>4}\n")
    for i, (aa, (x, y, z)) in enumerate(zip(sequence, coords.tolist()), start=1):
        handle.write(format_ca_atom(i, THREE_LETTER[aa], chain_id, i, x, y, z) + "\n")
    if sequence:
        last = THREE_LETTER[sequence[-1]]
        handle.write(f"{'TER':<6}{len(sequence) + 1:>5}      {last:>3} {chain_id:1}{len(sequence):>4}\n")
    handle.write("ENDMDL\n")


def write_pdb(path: str | Path, records: list[tuple[str, str, torch.Tensor]]) -> None:
    """Write (id, sequence, coords) records, one MODEL per protein."""
    with open(path, "w", encoding="ascii") as handle:
        for model_id, (seq_id, sequence, coords) in enumerate(records, start=1):
            handle.write(f"REMARK   1 MODEL {model_id} {seq_id}"[:80] + "\n")
        for model_id, (_, sequence, coords) in enumerate(records, start=1):
            write_model(handle, model_id, sequence, coords)
        handle.write("END\n")
