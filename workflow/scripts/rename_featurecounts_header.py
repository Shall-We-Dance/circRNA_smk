from __future__ import annotations

import argparse
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Rename featureCounts sample columns based on configured sample names."
    )
    parser.add_argument("--counts", required=True, help="Path to featureCounts output file")
    parser.add_argument(
        "--sample-names",
        required=True,
        help="Comma-separated sample names in the desired order",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    counts_path = Path(args.counts)
    sample_names = args.sample_names.split(",") if args.sample_names else []

    lines = counts_path.read_text().splitlines()
    if len(lines) < 2:
        raise ValueError("featureCounts output is malformed: %s" % counts_path)

    header = lines[1].split("\t")
    fixed_columns = 6
    if len(header) - fixed_columns != len(sample_names):
        raise ValueError(
            "featureCounts sample columns count does not match configured samples: "
            "%s vs %s" % (len(header) - fixed_columns, len(sample_names))
        )

    header[fixed_columns:] = sample_names
    lines[1] = "\t".join(header)
    counts_path.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
