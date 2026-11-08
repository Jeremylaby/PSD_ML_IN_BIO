# PSD_ML_IN_BIO

## Setup and shared data (DVC + Azure)

After cloning the repository, run the following commands from its root directory.
Install [uv](https://docs.astral.sh/uv/getting-started/installation/) first if needed.

```sh
uv sync
```

This installs the project dependencies, including DVC with Azure support and tqdm.
Python 3.13 or newer is required.

### Configure access to Azure

DVC is already initialized in this repository. The shared `.dvc/config` defines:

- Default remote: `storage`
- Storage account: `psdmlinbio2026`
- Blob container: `dvc` (`azure://dvc`)

Do not run `dvc init` or add the remote again after cloning. Get the Azure storage
account access key from the project owner privately, then configure it locally:

```sh
uv run dvc remote modify --local storage account_key 'YOUR_AZURE_STORAGE_ACCOUNT_KEY'
```

Replace the placeholder with the access key, not a connection string or SAS token.
The `--local` flag stores the key in `.dvc/config.local`, which is ignored by Git.
Never commit the key or paste it into the README or `.dvc/config`.

See the [DVC Azure authentication documentation](https://doc.dvc.org/user-guide/data-management/remote-storage/azure-blob-storage).

### Download the shared dataset

Once the project owner has uploaded the data, run:

```sh
uv run dvc pull
```

This restores the version of `data/` referenced by `data.dvc`. To get an updated
dataset later, run `git pull` followed by `uv run dvc pull`.

### Download structures from PDB

To download the structures directly from PDB instead, run:

```sh
uv run python utils/download_structures.py
```

To download only one split:

```sh
uv run python utils/download_structures.py train
```

The script reads IDs from `splited_protein_ids/` and saves whole PDB structures to
`data/protein_data/{train,eval,test}/`. Multiple entity IDs pointing to the same
PDB entry are downloaded once per split. The tqdm bar shows downloaded, skipped,
and failed files.

Existing nonempty `.cif` files are skipped. Press `Ctrl+C` to interrupt and rerun
the same command to continue. An interrupted file is downloaded again from the
beginning; temporary `.part` files are not treated as completed downloads.
Failed files are reported and can be retried by rerunning the command.

### Upload a new dataset version

After downloading or changing data, publish it to Azure and record its version:

```sh
uv run dvc add data
uv run dvc push
git add data.dvc .gitignore .dvc/config
git commit -m "Update protein dataset"
git push
```

Commit the project code and dependency files as usual. Git stores code, ID lists,
and DVC metadata; Azure stores the contents of `data/`.

## Original assignment and reference snippets

Protein structure prediction
Download protein clusters from https://cdn.rcsb.org/resources/sequence/clusters/clusters-by-entity-50.txt (one cluster per line; note that in a realistic scenario we would use 30% sequence identity clustering)
Consider clusters coming from specific source:
Team 1 (PDB): exclude protein IDs prepended with AF_ and MA_
Team 2 (AFDB): consider only protein IDs prepended with AF_
Select clusters with at least 10 members
Shuffle clusters (rows)
For each cluster select 10 random members
Assign whole clusters to splits: this should result in 4000/500/500 proteins in train/eval/test
Save IDs (one protein ID per line) to train_pdb.txt (train_af.txt) etc.
Download structure files using this script (download_structures.sh; note that both PDB and AFDB IDs need preprocessing):
```cmd
#!/usr/bin/env bash
# Usage: ./download_structures.sh train pdb
#        ./download_structures.sh eval  af

SPLIT=$1   # train / eval / test
DB=$2      # pdb / af
OUT="structures/${SPLIT}/${DB}"
mkdir -p "$OUT"

if [[ "$DB" == "pdb" ]]; then
    awk -F'_' '{print "https://files.rcsb.org/download/"$1".cif"}' "${SPLIT}_pdb.txt" \
        | sort -u | wget -i - -P "$OUT"

elif [[ "$DB" == "af" ]]; then
    python3 -c "
import re
for line in open('${SPLIT}_af.txt'):
    m = re.match(r'AF_AF(.+?)(F\d+)_\d+$', line.strip())
    if m:
        print(f'https://alphafold.ebi.ac.uk/files/AF-{m.group(1)}-{m.group(2)}-model_v6.cif')
" | sort -u | wget -i - -P "$OUT"
fi
```

   9. Team 1: split PDB files into separate chains using (extract_chain.py; this can be easily parallelized e.g. with xargs or GNU parallel):

```cmd #!/usr/bin/env python3
# Usage:
#   cat train_pdb.txt | python extract_chain.py structures/train/pdb chains/train/

import sys
from pathlib import Path
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB import MMCIFParser, MMCIFIO, Select

entry_id, in_dir, out_dir = sys.argv[1:4]
pdb_id, entity_id = entry_id.rsplit("_", 1)

d = MMCIF2Dict(f"{in_dir}/{pdb_id}.cif")
chains = {c for c, e in zip(d["_struct_asym.id"], d["_struct_asym.entity_id"]) if e == entity_id}

class Keep(Select):
    def accept_chain(self, chain): return chain.id in chains

st = MMCIFParser(QUIET=True).get_structure(pdb_id, f"{in_dir}/{pdb_id}.cif")
io = MMCIFIO()
io.set_structure(st)
Path(out_dir).mkdir(parents=True, exist_ok=True)
io.save(f"{out_dir}/{entry_id}.cif", Keep())
```
