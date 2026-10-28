import random
from pathlib import Path

MEMBERS_PER_CLUSTER = 10
SPLIT_PROTEINS = {"train": 4000, "eval": 500, "test": 500}
SEED = 42


def split_clusters(input_path: str = "filtered_ids.txt", output_dir: str = "splits") -> None:
	rng = random.Random(SEED)

	# 1. Select clusters with at least MEMBERS_PER_CLUSTER members
	with Path(input_path).open("r", encoding="utf-8") as f:
		clusters = [line.split() for line in f]
	clusters = [(idx, members) for idx, members in enumerate(clusters) if len(members) >= MEMBERS_PER_CLUSTER]

	# 2. Shuffle clusters
	rng.shuffle(clusters)

	needed = sum(SPLIT_PROTEINS.values()) // MEMBERS_PER_CLUSTER
	if len(clusters) < needed:
		raise ValueError(f"Need {needed} clusters, got {len(clusters)}")

	out = Path(output_dir)
	out.mkdir(exist_ok=True)
	start = 0
	with (out / "splits.tsv").open("w", encoding="utf-8") as f_tsv:
		f_tsv.write("protein_id\tcluster\tsplit\n")
		# 4. Assign whole clusters to splits
		for split, n_proteins in SPLIT_PROTEINS.items():
			n_clusters = n_proteins // MEMBERS_PER_CLUSTER
			with (out / f"{split}.txt").open("w", encoding="utf-8") as f_split:
				for cluster_idx, members in clusters[start:start + n_clusters]:
					# 3. Select MEMBERS_PER_CLUSTER random members of the cluster
					for protein_id in rng.sample(members, MEMBERS_PER_CLUSTER):
						f_split.write(f"{protein_id}\n")
						f_tsv.write(f"{protein_id}\t{cluster_idx}\t{split}\n")
			start += n_clusters


if __name__ == "__main__":
	split_clusters()
