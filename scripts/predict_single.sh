#!/usr/bin/env bash
# Predict the structure of a single protein from data/fasta/ by its ID.
#
# Usage:
#   scripts/predict_single.sh 7ZX8_2                # search train, eval and test
#   scripts/predict_single.sh 7ZX8_2 test           # search only one split
#   scripts/predict_single.sh 7ZX8_2 -- --seed 0    # extra options for the predictor
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
FASTA_DIR="$ROOT/data/fasta"

if [[ $# -lt 1 || "$1" == "-h" || "$1" == "--help" ]]; then
    sed -n '2,7p' "$0" | sed 's/^# \{0,1\}//'
    exit 1
fi

ID="$1"; shift
SPLITS=(train eval test)
if [[ $# -gt 0 && "$1" != "--" ]]; then
    SPLITS=("$1"); shift
fi
[[ $# -gt 0 && "$1" == "--" ]] && shift   # wszystko po "--" trafia do predyktora

if [[ ! -d "$FASTA_DIR" ]]; then
    echo "Brak $FASTA_DIR - najpierw uruchom: uv run python utils/cif_to_fasta.py" >&2
    exit 1
fi

# Plik tymczasowy z jednym białkiem; usuwamy go automatycznie na końcu.
TMP_DIR="$(mktemp -d)"
trap 'rm -rf "$TMP_DIR"' EXIT
TMP_FASTA="$TMP_DIR/$ID.fasta"

for split in "${SPLITS[@]}"; do
    file="$FASTA_DIR/$split.fasta"
    [[ -f "$file" ]] || continue
    # Wypisujemy tylko rekord, którego nagłówek to dokładnie ">ID".
    awk -v id="$ID" '/^>/{keep = ($1 == ">" id)} keep' "$file" > "$TMP_FASTA"
    if [[ -s "$TMP_FASTA" ]]; then
        echo "Found $ID in $split"
        cd "$ROOT"
        uv run python -m predictor.predict "$TMP_FASTA" --input-label "data/fasta/$split.fasta#$ID" "$@"
        exit $?
    fi
done

echo "Nie znaleziono $ID w: ${SPLITS[*]}" >&2
exit 1
