from __future__ import annotations

import argparse
import sys
from pathlib import Path


FIXED_COLUMNS = 6


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


def _find_header_index(lines: list[str]) -> int:
    for idx, line in enumerate(lines):
        if line.startswith("Geneid\t"):
            return idx
    raise ValueError("Could not locate featureCounts header line starting with 'Geneid'.")


def main() -> None:
    args = parse_args()
    counts_path = Path(args.counts)
    sample_names = [name for name in args.sample_names.split(",") if name]

    lines = counts_path.read_text().splitlines()
    if len(lines) < 2:
        raise ValueError(f"featureCounts output is malformed: {counts_path}")

    header_index = _find_header_index(lines)
    header = lines[header_index].split("\t")

    if len(header) <= FIXED_COLUMNS:
        raise ValueError(
            f"featureCounts header has too few columns ({len(header)}): {counts_path}"
        )

    observed_sample_cols = len(header) - FIXED_COLUMNS
    if observed_sample_cols != len(sample_names):
        print(
            "[rename_featurecounts_header] Warning: sample column count mismatch "
            f"(featureCounts={observed_sample_cols}, configured={len(sample_names)}). "
            "Renaming by position for overlapping columns only.",
            file=sys.stderr,
        )

    rename_count = min(observed_sample_cols, len(sample_names))
    for i in range(rename_count):
        header[FIXED_COLUMNS + i] = sample_names[i]

    lines[header_index] = "\t".join(header)
    counts_path.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
