"""Build the info.json metadata saved next to every prediction."""

import csv
import json
import platform
from datetime import datetime
from functools import lru_cache
from pathlib import Path

import torch

ROOT = Path(__file__).resolve().parents[1]
SPLITS_TSV = ROOT / "splited_protein_ids" / "splits.tsv"
STRUCTURES_DIR = ROOT / "data" / "protein_data"


@lru_cache(maxsize=1)
def load_splits(path: Path = SPLITS_TSV) -> dict[str, dict[str, str]]:
    """Map protein ID (e.g. 7ZX8_2) to its split and sequence cluster."""
    if not path.is_file():
        return {}
    with path.open(newline="", encoding="utf-8") as f:
        return {row["protein_id"]: row for row in csv.DictReader(f, delimiter="\t")}


def relative(path: Path) -> str:
    """Show paths inside the repository relative to it, others as given."""
    path = Path(path).resolve()
    return str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)


def protein_source(seq_id: str) -> dict:
    """Where a protein comes from: split, cluster, PDB entry, entity, CIF file.

    Proteins that are not in splits.tsv (e.g. your own sequences) only get
    their FASTA ID.
    """
    row = load_splits().get(seq_id)
    if row is None:
        return {"in_dataset": False}
    pdb_id, entity_id = seq_id.rsplit("_", 1)
    cif = STRUCTURES_DIR / row["split"] / f"{pdb_id}.cif"
    return {
        "in_dataset": True,
        "split": row["split"],
        "cluster": int(row["cluster"]),
        "pdb_id": pdb_id,
        "entity_id": entity_id,
        "cif_file": relative(cif) if cif.is_file() else None,
        "rcsb_url": f"https://www.rcsb.org/structure/{pdb_id}",
    }


def build_info(*, input_fasta: str, output_pdb: str, max_len: int, batch_size: int,
               seed: int, proteins: list[dict]) -> dict:
    return {
        "created": datetime.now().isoformat(timespec="seconds"),
        "input_fasta": input_fasta,
        "output_pdb": relative(Path(output_pdb)),
        "model": "DummyStructureModel (random coordinates)",
        "parameters": {"max_len": max_len, "batch_size": batch_size, "seed": seed},
        "environment": {"python": platform.python_version(), "torch": torch.__version__},
        "n_proteins": len(proteins),
        "n_ca_atoms": sum(p["processed_length"] for p in proteins),
        "proteins": proteins,
    }


def write_info(path: str | Path, info: dict) -> None:
    Path(path).write_text(json.dumps(info, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
