#!/usr/bin/env python3
"""Select the top-N circRNAs from a DCC CircRNACount table to feed circtools
primex. Selection criterion: total BSJ count across all samples (descending);
ties broken by genomic position for reproducibility.

Why this default: primex requires either a gene-name list or a circRNA-ID
list. The most universally useful "validate these" list is the most reliably
detected circRNAs (highest BSJ support across the dataset), and that signal
is independent of organism, comparisons, or whether circtest ran.

Output format: one circRNA ID per line in primex's expected shape
'GeneName_Chr_Start_End_Strand'. When gene or strand are unknown we use a
single dot, matching the placeholder DCC uses in CircCoordinates.
"""
import argparse
import csv
import sys
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--circ-count", required=True,
                   help="Path to a DCC CircRNACount table.")
    p.add_argument("--circ-coord",
                   help="Optional path to the matching CircCoordinates file; "
                        "used to pull gene name and strand for the primex ID.")
    p.add_argument("--top-n", type=int, default=25,
                   help="Number of candidates to select [default 25].")
    p.add_argument("--out", required=True, help="Output candidate-ID file.")
    return p.parse_args()


def read_coords(path):
    """Return {(chrom, start, end) -> (gene, strand)} from CircCoordinates.

    DCC writes the columns: Chr, Start, End, Gene, JunctionType, Strand.
    Missing fields are filled with '.'.
    """
    if not path:
        return {}
    coord_file = Path(path)
    if not coord_file.exists():
        return {}
    out = {}
    with coord_file.open() as fh:
        reader = csv.reader(fh, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            return {}
        # Header tolerance: allow either an actual header or an unlabeled
        # numeric first row (some DCC builds omit the header). If header[0]
        # looks like a chromosome (no leading 'Chr' label), pretend we already
        # consumed the row.
        if header and header[0].lower() not in {"chr", "chrom", "chromosome"}:
            try:
                start = int(header[1])
                end = int(header[2])
                gene = header[3] if len(header) > 3 else "."
                strand = header[5] if len(header) > 5 else "."
                out[(header[0], start, end)] = (gene or ".", strand or ".")
            except (ValueError, IndexError):
                pass  # first row was likely a real header after all
        for row in reader:
            if len(row) < 3:
                continue
            try:
                start = int(row[1])
                end = int(row[2])
            except ValueError:
                continue
            gene = row[3] if len(row) > 3 else "."
            strand = row[5] if len(row) > 5 else "."
            out[(row[0], start, end)] = (gene or ".", strand or ".")
    return out


def main():
    args = parse_args()
    circ_count_path = Path(args.circ_count)
    if not circ_count_path.exists() or circ_count_path.stat().st_size == 0:
        raise SystemExit(f"ERROR: CircRNACount not found or empty: {circ_count_path}")

    # Pair the coordinates table (if any) automatically with CircRNACount when
    # the caller didn't pass one explicitly. They live in the same directory.
    coord_path = args.circ_coord
    if not coord_path:
        candidate = circ_count_path.parent / "CircCoordinates"
        if candidate.exists():
            coord_path = str(candidate)
    coords = read_coords(coord_path)

    # CircRNACount: Chr, Start, End, sample1, sample2, ... (or with a Strand
    # column when present). We sum every column from index 3 onwards that
    # contains an integer, skipping any string-valued columns (Strand).
    entries = []
    with circ_count_path.open() as fh:
        reader = csv.reader(fh, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            raise SystemExit(f"ERROR: CircRNACount is empty: {circ_count_path}")

        # Identify numeric columns from the header by sniffing the first data
        # row. Anything that fails int-cast on row 1 is treated as non-count.
        first_row = next(reader, None)
        if first_row is None:
            raise SystemExit(f"ERROR: CircRNACount has no data rows: {circ_count_path}")

        count_cols = []
        for i in range(3, len(first_row)):
            try:
                int(first_row[i])
                count_cols.append(i)
            except ValueError:
                continue
        if not count_cols:
            raise SystemExit(f"ERROR: no numeric sample columns found in {circ_count_path}")

        def process(row):
            try:
                chrom = row[0]
                start = int(row[1])
                end = int(row[2])
                total = sum(int(row[i]) for i in count_cols if i < len(row))
            except (ValueError, IndexError):
                return None
            return (chrom, start, end, total)

        first_processed = process(first_row)
        if first_processed:
            entries.append(first_processed)
        for row in reader:
            processed = process(row)
            if processed:
                entries.append(processed)

    # Sort by total BSJ desc, then by chrom/start/end asc for reproducibility.
    entries.sort(key=lambda x: (-x[3], x[0], x[1], x[2]))
    top = entries[: max(0, int(args.top_n))]
    if not top:
        raise SystemExit("ERROR: no circRNA entries to select from (top-N empty).")

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as out_fh:
        for chrom, start, end, total in top:
            gene, strand = coords.get((chrom, start, end), (".", "."))
            # primex expects IDs like 'GeneName_Chr_Start_End_Strand'.
            out_fh.write(f"{gene}_{chrom}_{start}_{end}_{strand}\n")
    print(f"[select_primex_candidates] wrote {len(top)} candidate IDs to {out_path}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
