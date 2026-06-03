#!/usr/bin/env python3
"""Subset a DCC-format directory (CircRNACount, CircCoordinates, LinearCount)
to a chosen ordered list of samples, producing a new DCC directory containing
only those samples.

This is the common path that prepares per-comparison inputs for circtools'
circtest and quickcheck rules, regardless of whether the canonical DCC came
from `circtools detect` (real DCC) or `ciri3_to_dcc_format.py` (synthesized
from CIRI3 matrices). Both layouts are supported because we look up sample
columns by header NAME, not by position.

Output layout invariants we MUST produce (for circtools 2.0.5 compatibility):

  CircRNACount header: Chr<TAB>Start<TAB>End<TAB>sample1<TAB>sample2<...>
                       (NO Strand column -- see "Strand normalization" below)
  LinearCount header:  Chr<TAB>Start<TAB>End<TAB>sample1<TAB>sample2<...>
                       (no Strand column, same as input format)
  CircCoordinates:     copied verbatim (no per-sample columns to trim)

Strand normalization (why we drop CircRNACount.Strand):
   Real `circtools detect` writes CircRNACount with a Strand column (between
   End and the sample columns) but LinearCount without it. circtools' own
   circtest then breaks on that asymmetry because its column slicing uses
       circle_description = list(range(3))   # hardcoded
       all_cols = circle_description + cond_cols
       circ_slice  = CircRNACount.iloc[:, all_cols]
       linear_slice = LinearCount.iloc[:, all_cols]
   The same indices are applied to both files, so the `-c` flag misaligns
   between them when CircRNACount has Strand. We strip Strand from
   CircRNACount on output so BOTH count files have identical 3-column
   metadata. Strand information remains preserved in CircCoordinates, which
   is where primex reads it.

   This also keeps quickcheck consistent: its `data_start` detection
   (`4 if any(col == 'strand' for col in base_cols) else 3`) will now read 3
   from both files instead of seeing the asymmetric 4/3.
"""

import argparse
import csv
import sys
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(
        description="Subset a DCC dir to a specific ordered set of samples.",
    )
    p.add_argument("--in-dir", required=True,
                   help="Source DCC directory (must contain CircRNACount, "
                        "CircCoordinates, LinearCount).")
    p.add_argument("--out-dir", required=True,
                   help="Destination DCC directory; will be created.")
    p.add_argument("--samples", required=True,
                   help="Comma-separated sample names, in the desired column order. "
                        "All listed names must exist as headers in CircRNACount "
                        "AND LinearCount.")
    return p.parse_args()


# The four reserved-header tokens that appear before sample columns in DCC
# count tables. Anything in the header line that isn't one of these is treated
# as a sample column name.
_META_TOKENS = {"Chr", "Start", "End", "Strand"}

# DCC-derived sample column headers carry the basename of whatever file was fed
# into `circtools detect` for that sample. CircRNACount columns come from the
# *_Chimeric.out.junction sample sheet, while LinearCount columns may instead
# carry the *_Aligned.sortedByCoord.out.bam basename from the -B sheet. Neither
# matches the bare sample name the rest of the pipeline uses, so we strip these
# known suffixes to recover the sample name. Listed longest-first so the most
# specific suffix wins. The CIRI3-synthesized DCC path (ciri3_to_dcc_format.py)
# already emits bare sample names, which pass through unchanged.
_DCC_SAMPLE_SUFFIXES = (
    "_Chimeric.out.junction",
    ".Chimeric.out.junction",
    "_Aligned.sortedByCoord.out.bam",
    ".Aligned.sortedByCoord.out.bam",
    "Chimeric.out.junction",
    "Aligned.sortedByCoord.out.bam",
)


def _normalize_sample_name(col):
    """Recover the bare sample name from a DCC sample-column header.

    Strips any directory component and a single known DCC filename suffix, then
    trims a trailing separator. A header that carries none of these suffixes
    (e.g. a CIRI3-synthesized bare sample name) is returned unchanged.
    """
    name = col.rsplit("/", 1)[-1].rsplit("\\", 1)[-1]
    for suffix in _DCC_SAMPLE_SUFFIXES:
        if name.endswith(suffix) and len(name) > len(suffix):
            return name[: -len(suffix)].rstrip("_.")
    return name


def _split_header(header_row):
    """Return (meta_cols_in_order, sample_cols_in_order) from a count-file header."""
    meta = []
    samples = []
    seen_sample_section = False
    for col in header_row:
        if not seen_sample_section and col in _META_TOKENS:
            meta.append(col)
        else:
            # Once we hit a non-meta column, everything after is a sample.
            # This guards against a sample whose name happens to be e.g. "End"
            # by requiring meta columns to appear contiguously at the start.
            seen_sample_section = True
            samples.append(col)
    return meta, samples


def subset_count_file(in_path: Path, out_path: Path, want_samples, drop_strand: bool):
    """Subset and reorder sample columns of CircRNACount / LinearCount.

    Preserves the row order of the input file exactly. If ``drop_strand`` is
    True, any Strand metadata column is removed from the output so both count
    files end with the canonical 3-column metadata layout (Chr/Start/End).
    Raises SystemExit on missing samples.
    """
    with in_path.open() as fh:
        reader = csv.reader(fh, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            raise SystemExit(f"ERROR: {in_path} is empty")

        meta_cols, sample_cols = _split_header(header)

        # Map each requested bare sample name to its column index. DCC headers
        # may be suffixed (e.g. "Mock_R1_Chimeric.out.junction") so we register
        # BOTH the raw header and its normalized bare name as lookup keys. Raw
        # names are registered first and win ties; normalized names are only
        # added when they don't shadow a different column, and any genuinely
        # ambiguous normalized key is recorded so we can fail loudly if (and
        # only if) a requested sample actually hits it.
        sample_to_col_idx = {}
        ambiguous_keys = set()

        def _register(key, idx):
            if key in sample_to_col_idx:
                if sample_to_col_idx[key] != idx:
                    ambiguous_keys.add(key)
                return
            sample_to_col_idx[key] = idx

        for i, raw_name in enumerate(sample_cols):
            col_idx = len(meta_cols) + i
            _register(raw_name, col_idx)
        for i, raw_name in enumerate(sample_cols):
            col_idx = len(meta_cols) + i
            norm_name = _normalize_sample_name(raw_name)
            if norm_name != raw_name:
                _register(norm_name, col_idx)

        ambiguous_wanted = [s for s in want_samples if s in ambiguous_keys]
        if ambiguous_wanted:
            raise SystemExit(
                f"ERROR: {in_path} has ambiguous sample column(s) for: "
                f"{', '.join(ambiguous_wanted)}. Multiple headers normalize to "
                f"the same sample name. Headers seen: {', '.join(sample_cols)}"
            )

        missing = [s for s in want_samples if s not in sample_to_col_idx]
        if missing:
            raise SystemExit(
                f"ERROR: {in_path} is missing required sample column(s): "
                f"{', '.join(missing)}. Available sample columns: "
                f"{', '.join(sample_cols)}"
            )

        want_idx = [sample_to_col_idx[s] for s in want_samples]
        meta_len = len(meta_cols)

        # Build the indices of metadata columns we keep on output. If
        # ``drop_strand`` is set, omit any column whose name is "Strand"
        # (case-insensitive). The remaining metadata columns retain their
        # original relative order.
        if drop_strand:
            keep_meta_idx = [
                i for i, name in enumerate(meta_cols)
                if name.lower() != "strand"
            ]
        else:
            keep_meta_idx = list(range(meta_len))
        kept_meta_names = [meta_cols[i] for i in keep_meta_idx]

        out_path.parent.mkdir(parents=True, exist_ok=True)
        with out_path.open("w") as out_fh:
            writer = csv.writer(out_fh, delimiter="\t", lineterminator="\n")
            writer.writerow(kept_meta_names + list(want_samples))
            for row in reader:
                if not row:
                    continue
                # Pad short rows defensively (DCC sometimes truncates trailing
                # empty fields). Use "0" for sample columns and "" for meta.
                while len(row) < meta_len:
                    row.append("")
                while len(row) < meta_len + len(sample_cols):
                    row.append("0")
                out_row = [row[i] for i in keep_meta_idx] + [row[i] for i in want_idx]
                writer.writerow(out_row)


def copy_coordinates(in_path: Path, out_path: Path):
    """CircCoordinates has no per-sample columns; copy verbatim."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_bytes(in_path.read_bytes())


def main():
    args = parse_args()
    samples = [s for s in args.samples.split(",") if s]
    if len(samples) < 2:
        raise SystemExit("ERROR: at least two samples are required.")
    if len(set(samples)) != len(samples):
        raise SystemExit("ERROR: --samples contains duplicate names.")

    in_dir = Path(args.in_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    crc_in = in_dir / "CircRNACount"
    coord_in = in_dir / "CircCoordinates"
    lin_in = in_dir / "LinearCount"
    for p in (crc_in, coord_in, lin_in):
        if not p.exists():
            raise SystemExit(f"ERROR: required DCC file is missing: {p}")

    # Drop Strand from both count files so circtest's hardcoded
    # circle_description=[0,1,2] indexes both files correctly. CircCoordinates
    # retains Strand because it's not column-sliced the same way.
    subset_count_file(crc_in, out_dir / "CircRNACount", samples, drop_strand=True)
    subset_count_file(lin_in, out_dir / "LinearCount", samples, drop_strand=True)
    copy_coordinates(coord_in, out_dir / "CircCoordinates")

    print(f"[subset_dcc] subset to {len(samples)} samples: {', '.join(samples)}",
          file=sys.stderr)
    print(f"[subset_dcc] wrote {out_dir / 'CircRNACount'}", file=sys.stderr)
    print(f"[subset_dcc] wrote {out_dir / 'CircCoordinates'}", file=sys.stderr)
    print(f"[subset_dcc] wrote {out_dir / 'LinearCount'}", file=sys.stderr)


if __name__ == "__main__":
    main()
