"""CLI: predict CA coordinates for every sequence in a FASTA file.

Each run is saved to its own folder: outputs/run_<timestamp>/predicted.pdb + info.json

Usage:
    uv run python -m predictor.predict examples/example.fasta
    uv run python -m predictor.predict examples/example.fasta -o my.pdb
"""

import argparse
from datetime import datetime
from pathlib import Path

import torch

from predictor.fasta import read_fasta
from predictor.model import DummyStructureModel
from predictor.pdb_writer import write_pdb
from predictor.run_info import build_info, protein_source, relative, write_info
from predictor.vocab import encode_batch


def make_run_dir(base: str | Path = "outputs") -> Path:
    """Create a fresh outputs/run_<timestamp>/ folder for one prediction run."""
    stamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    run_dir = Path(base) / f"run_{stamp}"
    # Dwa uruchomienia w tej samej sekundzie nie nadpiszą się nawzajem.
    suffix = 1
    while run_dir.exists():
        run_dir = Path(base) / f"run_{stamp}_{suffix}"
        suffix += 1
    run_dir.mkdir(parents=True)
    return run_dir


def predict(fasta_path: str, output_path: str, max_len: int = 1024,
            batch_size: int = 16, seed: int | None = None,
            info_path: str | None = None, input_label: str | None = None) -> int:
    """Run the dummy model over a FASTA file and write a PDB. Returns #proteins.

    If ``info_path`` is given, also writes run metadata (see run_info.py).
    ``input_label`` is the input name stored there (defaults to ``fasta_path``).
    """
    if seed is None:
        # Losujemy seed i zapisujemy go w info.json, żeby każdy run dało się powtórzyć.
        seed = int(torch.randint(0, 2**31 - 1, (1,)))
    torch.manual_seed(seed)  # ten sam seed = te same "losowe" współrzędne
    records = read_fasta(fasta_path)
    model = DummyStructureModel().eval()

    results, proteins = [], []
    # Białka przetwarzamy paczkami (batchami) po batch_size sekwencji naraz.
    for batch_index, start in enumerate(range(0, len(records), batch_size)):
        batch = records[start : start + batch_size]
        x, mask, lengths = encode_batch([seq for _, seq in batch], max_len=max_len)
        for pos, ((seq_id, seq), coords) in enumerate(zip(batch, model.predict(x, mask))):
            results.append((seq_id, seq, coords))
            proteins.append({
                "id": seq_id,
                "fasta_index": start + pos,          # kolejność w pliku FASTA (od 0)
                "pdb_model": len(results),           # numer bloku MODEL w pliku PDB (od 1)
                "batch_index": batch_index,
                "position_in_batch": pos,
                "sequence_length": len(seq),
                "processed_length": lengths[pos],
                "truncated": len(seq) > max_len,
                "source": protein_source(seq_id),
            })

    write_pdb(output_path, results)
    if info_path is not None:
        write_info(info_path, build_info(
            input_fasta=input_label or relative(Path(fasta_path)), output_pdb=output_path,
            max_len=max_len, batch_size=batch_size, seed=seed, proteins=proteins,
        ))
    return len(results)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("fasta", help="input multi-protein FASTA file")
    parser.add_argument("-o", "--output", default=None,
                        help="output PDB file (default: outputs/run_<timestamp>/predicted.pdb)")
    parser.add_argument("--outputs-dir", default="outputs", help="base folder for run folders")
    parser.add_argument("--max-len", type=int, default=1024, help="maximum sequence length D")
    parser.add_argument("--batch-size", type=int, default=16, help="batch size B")
    parser.add_argument("--seed", type=int, default=None, help="random seed (default: random, saved in info.json)")
    parser.add_argument("--input-label", default=None, help=argparse.SUPPRESS)  # used by predict_single.sh
    args = parser.parse_args()

    if args.output:
        output = Path(args.output)
        info = output.with_suffix(".info.json")
    else:
        run_dir = make_run_dir(args.outputs_dir)
        output, info = run_dir / "predicted.pdb", run_dir / "info.json"
    n = predict(args.fasta, str(output), args.max_len, args.batch_size, args.seed,
                info_path=str(info), input_label=args.input_label)
    print(f"Wrote {n} protein(s) to {output} (metadata: {info})")


if __name__ == "__main__":
    main()
