from pathlib import Path

EXCLUDED_PREFIXES = ("AF_", "MA_")
MIN_CLUSTER_SIZE = 10


def filter_ids(input_path: str = "ids.txt", output_path: str = "filtered_ids.txt") -> None:
	input_file = Path(input_path)
	output_file = Path(output_path)

	# Each line is one cluster of space-separated IDs, so filter every ID
	# individually and keep the line structure; drop clusters that have
	# fewer than MIN_CLUSTER_SIZE IDs left after filtering.
	with input_file.open("r", encoding="utf-8") as f_in, output_file.open("w", encoding="utf-8") as f_out:
		for line in f_in:
			kept = [protein_id for protein_id in line.split() if not protein_id.startswith(EXCLUDED_PREFIXES)]
			if len(kept) >= MIN_CLUSTER_SIZE:
				f_out.write(" ".join(kept) + "\n")


if __name__ == "__main__":
	filter_ids()
