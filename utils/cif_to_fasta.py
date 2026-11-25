"""Extract entity sequences from downloaded mmCIF files into FASTA files.

Usage:
    uv run python utils/cif_to_fasta.py            # all splits
    uv run python utils/cif_to_fasta.py test       # one split

Reads IDs (e.g. 10ZN_1 = PDB entry 10ZN, entity 1) from splited_protein_ids/{split}.txt,
takes the full canonical sequence of that entity from data/protein_data/{split}/{PDB}.cif
and writes data/fasta/{split}.fasta.
"""

import argparse
import io
import sys
from pathlib import Path

from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from tqdm import tqdm

ROOT = Path(__file__).resolve().parents[1]
SPLITS = ("train", "eval", "test")


def as_list(value: str | list[str]) -> list[str]:
    return value if isinstance(value, list) else [value]


def read_header(cif_path: str | Path) -> str:
    """Return the file text up to the _atom_site block.

    Sequences live in _entity_poly near the top of the file; the coordinates
    below make up most of the file, so skipping them makes parsing much faster.
    """
    lines = []
    with open(cif_path, encoding="utf-8", errors="replace") as f:
        for line in f:
            # Dalej są już tylko współrzędne atomów - nie są nam tu potrzebne.
            if line.startswith("_atom_site."):
                if lines and lines[-1].strip() == "loop_":
                    lines.pop()
                break
            lines.append(line)
    return "".join(lines)


def entity_sequences(cif_path: str | Path) -> dict[str, str]:
    """Return {entity_id: canonical one-letter sequence} for protein entities.

    Uses _entity_poly.pdbx_seq_one_letter_code_can: the full sequence of the
    entity (also residues not observed in the experiment), with modified
    residues mapped to their standard parent (e.g. MSE -> M).
    """
    d = MMCIF2Dict(io.StringIO(read_header(cif_path)))
    ids = as_list(d.get("_entity_poly.entity_id", []))
    types = as_list(d.get("_entity_poly.type", []))
    seqs = as_list(d.get("_entity_poly.pdbx_seq_one_letter_code_can", []))
    result = {}
    for entity_id, poly_type, seq in zip(ids, types, seqs):
        if poly_type != "polypeptide(L)":
            continue  # pomijamy DNA/RNA i inne polimery
        # W pliku sekwencja jest łamana na wiele linii - sklejamy ją w jeden ciąg.
        result[entity_id] = "".join(seq.split()).upper()
    return result


def write_fasta(path: Path, records: list[tuple[str, str]], width: int = 60) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="ascii") as f:
        for seq_id, seq in records:
            f.write(f">{seq_id}\n")
            for i in range(0, len(seq), width):
                f.write(seq[i : i + width] + "\n")


def convert_split(split: str, structures_dir: Path, output_dir: Path) -> int:
    """Write {output_dir}/{split}.fasta. Returns the number of IDs that failed."""
    ids = [line.strip() for line in (ROOT / "splited_protein_ids" / f"{split}.txt").read_text().splitlines() if line.strip()]
    cache: dict[str, dict[str, str]] = {}  # jeden plik CIF czytamy tylko raz
    records, failed = [], 0
    for protein_id in tqdm(ids, desc=split, unit="id"):
        pdb_id, entity_id = protein_id.rsplit("_", 1)
        cif = structures_dir / split / f"{pdb_id}.cif"
        try:
            if pdb_id not in cache:
                cache[pdb_id] = entity_sequences(cif)
            seq = cache[pdb_id].get(entity_id)
            if not seq:
                raise ValueError(f"no protein entity {entity_id} in {cif.name}")
            records.append((protein_id, seq))
        except (OSError, ValueError) as error:
            failed += 1
            tqdm.write(f"{split}/{protein_id}: {error}")
    out = output_dir / f"{split}.fasta"
    write_fasta(out, records)
    tqdm.write(f"{split}: wrote {len(records)} sequences to {out.relative_to(ROOT) if out.is_relative_to(ROOT) else out}, {failed} failed")
    return failed


def main() -> int:
    parser = argparse.ArgumentParser(description="Extract entity sequences from mmCIF files to FASTA.")
    parser.add_argument("splits", nargs="*", choices=SPLITS, help="Default: train eval test")
    parser.add_argument("--structures-dir", type=Path, default=ROOT / "data" / "protein_data")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "data" / "fasta")
    args = parser.parse_args()
    failures = sum(convert_split(s, args.structures_dir, args.output_dir) for s in (args.splits or SPLITS))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
