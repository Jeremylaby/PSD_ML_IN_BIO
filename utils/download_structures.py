import argparse
import shutil
import sys
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import urlopen

from tqdm import tqdm

ROOT = Path(__file__).resolve().parents[1]


def download_split(split: str) -> int:
    input_path = ROOT / "splited_protein_ids" / f"{split}.txt"
    output_dir = ROOT / "data" / "protein_data" / split
    pdb_ids = sorted({line.strip().split("_")[0] for line in input_path.read_text().splitlines() if line.strip()})
    output_dir.mkdir(parents=True, exist_ok=True)
    downloaded = skipped = failed = 0

    with tqdm(pdb_ids, desc=split, unit="file") as progress:
        for pdb_id in progress:
            destination = output_dir / f"{pdb_id}.cif"
            temporary = destination.with_suffix(".cif.part")
            if destination.is_file() and destination.stat().st_size > 0:
                skipped += 1
            else:
                try:
                    url = f"https://files.rcsb.org/download/{pdb_id}.cif"
                    with urlopen(url, timeout=60) as response, temporary.open("wb") as output:
                        shutil.copyfileobj(response, output)
                        expected_size = response.headers.get("Content-Length")
                    size = temporary.stat().st_size
                    if size == 0 or (expected_size is not None and size != int(expected_size)):
                        raise OSError("Incomplete or empty file")
                    temporary.replace(destination)
                    downloaded += 1
                except (OSError, HTTPError, URLError) as error:
                    failed += 1
                    tqdm.write(f"{split}/{pdb_id}: {error}")
                finally:
                    temporary.unlink(missing_ok=True)
            progress.set_postfix(downloaded=downloaded, skipped=skipped, errors=failed)

    return failed


def main() -> int:
    parser = argparse.ArgumentParser(description="Download PDB structures; skip completed files.")
    parser.add_argument("splits", nargs="*", choices=("train", "eval", "test"), help="Default: train eval test")
    args = parser.parse_args()
    try:
        failures = sum(download_split(split) for split in (args.splits or ("train", "eval", "test")))
    except KeyboardInterrupt:
        tqdm.write("Interrupted. Run again to download missing files.")
        return 130
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
