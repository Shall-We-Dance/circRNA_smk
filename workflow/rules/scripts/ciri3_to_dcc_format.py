#!/usr/bin/env python3
"""Synthesize DCC-format CircRNACount/CircCoordinates/LinearCount tables from
CIRI3 BSJ_Matrix + FSJ_Matrix outputs.

DCC layout we write here (3 metadata columns in BOTH count tables -- see
"layout caveat" below):

  CircRNACount:    Chr<TAB>Start<TAB>End<TAB>sample1<TAB>sample2<...>
                   (one row per circRNA; integer BSJ counts per sample)
  CircCoordinates: Chr<TAB>Start<TAB>End<TAB>Gene<TAB>JunctionType<TAB>Strand
                   (one row per circRNA, same row order as CircRNACount)
  LinearCount:     Chr<TAB>Start<TAB>End<TAB>sample1<TAB>sample2<...>
                   (one row per circRNA, same row order, integer FSJ/host counts)

Layout caveat (why we DROP Strand from CircRNACount):
   circtools detect writes CircRNACount with a Strand column between End and
   the sample columns, but LinearCount without it (see
   circtools/detect/CombineCounts.py::writeouput vs ::writeouput_linear).
   However, circtools circtest's wrapper (circtools/circtest/
   circtools_circtest_wrapper.py) hardcodes
       circle_description = list(range(3))
       all_cols = circle_description + cond_cols
       circ_slice  = CircRNACount.iloc[:, all_cols]
       linear_slice = LinearCount.iloc[:, all_cols]
   so it applies the SAME column indices to both files. If CircRNACount has a
   Strand column and LinearCount doesn't, those indices misalign and circtest
   raises 'positional indexers are out-of-bounds'.
   The workaround is to keep BOTH count files at 3 metadata columns. We carry
   strand information in CircCoordinates instead, which is where downstream
   circtools rules (primex) look for it anyway.

CIRI3's BSJ/FSJ matrices have rows keyed by "chr:start|end" with a leading
circRNA_ID column followed by per-sample integer columns; we parse the ID into
chr/start/end and pair the two matrices on identical IDs. Strand and gene come
from all_samples.ciri3 (merged CIRI3 detection table) when present.
"""

import argparse
import csv
import re
import sys
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(
        description="Convert CIRI3 BSJ/FSJ matrices to DCC-format CircRNACount + "
        "CircCoordinates + LinearCount.",
    )
    p.add_argument("--bsj-matrix", required=True,
                   help="all_samples.ciri3.BSJ_Matrix (CIRI3 merged BSJ matrix).")
    p.add_argument("--fsj-matrix", required=True,
                   help="all_samples.ciri3.FSJ_Matrix (CIRI3 merged FSJ matrix).")
    p.add_argument("--ciri3-annot",
                   help="Optional all_samples.ciri3 to source gene/strand annotation.")
    p.add_argument("--samples", required=True,
                   help="Comma-separated sample names, in the desired column order.")
    p.add_argument("--out-dir", required=True,
                   help="Output directory; will contain CircRNACount, CircCoordinates, LinearCount.")
    return p.parse_args()


_ID_RE = re.compile(r"^(?P<chr>[^:]+):(?P<start>\d+)[|\-](?P<end>\d+)$")


def parse_circ_id(circ_id):
    """Parse a CIRI3 circRNA_ID like 'chr1:12345|67890' into (chr, start, end).

    Supports both '|' and '-' as the separator between start and end, since
    different CIRI3 builds have used both. Returns None on parse failure so the
    caller can skip / warn cleanly rather than crash mid-file.
    """
    m = _ID_RE.match(circ_id.strip())
    if not m:
        return None
    return m.group("chr"), int(m.group("start")), int(m.group("end"))


def read_matrix(path, samples):
    """Read a CIRI3 BSJ/FSJ matrix into {circ_id: {sample: int_count}}.

    The matrix header must include 'circRNA_ID' as col 0 plus one column per
    sample we care about. We tolerate extra columns and out-of-order samples;
    columns not in `samples` are silently ignored.
    """
    path = Path(path)
    with path.open() as fh:
        reader = csv.reader(fh, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            raise SystemExit(f"ERROR: {path} is empty")

        sample_to_idx = {col: i for i, col in enumerate(header)}
        missing = [s for s in samples if s not in sample_to_idx]
        if missing:
            raise SystemExit(
                f"ERROR: matrix {path} is missing required sample column(s): "
                f"{', '.join(missing)}. Available columns: {', '.join(header[1:])}"
            )

        out = {}
        for row in reader:
            if not row:
                continue
            circ_id = row[0]
            counts = {}
            for s in samples:
                idx = sample_to_idx[s]
                try:
                    counts[s] = int(float(row[idx])) if idx < len(row) and row[idx] != "" else 0
                except ValueError:
                    counts[s] = 0
            out[circ_id] = counts
        return out


def read_ciri3_annotation(path):
    """Extract {circ_id -> (gene, strand, junction_type)} from all_samples.ciri3.

    Returns an empty dict if path is None or unreadable; the caller falls back
    to '.' for gene and strand in that case. Junction type defaults to 'NA'
    since CIRI3's annotation is gene-level, not splice-class-level.
    """
    if not path:
        return {}
    path = Path(path)
    if not path.exists():
        return {}
    out = {}
    with path.open() as fh:
        reader = csv.reader(fh, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            return {}

        def col(name):
            return header.index(name) if name in header else -1

        c_id = col("circRNA_ID")
        c_strand = col("strand")
        c_gene = col("gene_id")
        c_type = col("circRNA_type")
        if c_id < 0:
            return {}

        for row in reader:
            if not row or len(row) <= c_id:
                continue
            cid = row[c_id]
            if cid in out:
                continue  # first-row wins (merged table has one row per sample)
            gene = row[c_gene] if c_gene >= 0 and c_gene < len(row) else "."
            strand = row[c_strand] if c_strand >= 0 and c_strand < len(row) else "."
            jtype = row[c_type] if c_type >= 0 and c_type < len(row) else "NA"
            out[cid] = (gene or ".", strand or ".", jtype or "NA")
        return out


def write_dcc_tables(circs, bsj_counts, fsj_counts, annot, samples, out_dir):
    """Write the three DCC files matching CombineCounts.writeouput layout."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    crc_path = out_dir / "CircRNACount"
    coord_path = out_dir / "CircCoordinates"
    lin_path = out_dir / "LinearCount"

    with crc_path.open("w") as crc, coord_path.open("w") as coord, lin_path.open("w") as lin:
        # CircRNACount and LinearCount both use 3 metadata columns
        # (Chr, Start, End) -- see "layout caveat" in module docstring for why
        # we drop Strand from CircRNACount.
        crc.write("Chr\tStart\tEnd\t" + "\t".join(samples) + "\n")
        coord.write("Chr\tStart\tEnd\tGene\tJunctionType\tStrand\n")
        lin.write("Chr\tStart\tEnd\t" + "\t".join(samples) + "\n")

        n_written = 0
        for circ_id in circs:
            parsed = parse_circ_id(circ_id)
            if parsed is None:
                continue
            chrom, start, end = parsed
            gene, strand, jtype = annot.get(circ_id, (".", ".", "NA"))

            bsj = bsj_counts.get(circ_id, {})
            fsj = fsj_counts.get(circ_id, {})
            bsj_row = [str(int(bsj.get(s, 0))) for s in samples]
            # Some circRNAs may have a BSJ row but no FSJ row (FSJ matrix only
            # carries loci with detected linear reads). Default to 0 so the row
            # order between CircRNACount and LinearCount stays aligned.
            fsj_row = [str(int(fsj.get(s, 0))) for s in samples]

            crc.write(f"{chrom}\t{start}\t{end}\t" + "\t".join(bsj_row) + "\n")
            coord.write(f"{chrom}\t{start}\t{end}\t{gene}\t{jtype}\t{strand}\n")
            lin.write(f"{chrom}\t{start}\t{end}\t" + "\t".join(fsj_row) + "\n")
            n_written += 1

    return crc_path, coord_path, lin_path, n_written


def main():
    args = parse_args()
    samples = [s for s in args.samples.split(",") if s]
    if len(samples) < 2:
        raise SystemExit("ERROR: at least two samples are required.")

    bsj = read_matrix(args.bsj_matrix, samples)
    fsj = read_matrix(args.fsj_matrix, samples)
    annot = read_ciri3_annotation(args.ciri3_annot) if args.ciri3_annot else {}

    # Union of circRNA IDs across the two matrices, sorted for reproducibility.
    # We use the BSJ-matrix order primarily (the canonical circRNA set), then
    # append any FSJ-only loci so the linear-only signal isn't silently dropped.
    bsj_ids = list(bsj.keys())
    fsj_extra = [cid for cid in fsj.keys() if cid not in bsj]
    circs = sorted(bsj_ids) + sorted(fsj_extra)

    crc, coord, lin, n_kept = write_dcc_tables(circs, bsj, fsj, annot, samples, args.out_dir)
    n_skipped = len(circs) - n_kept
    print(f"[ciri3_to_dcc_format] BSJ ids={len(bsj_ids)} FSJ-only extras={len(fsj_extra)} "
          f"rows_written={n_kept} unparseable_skipped={n_skipped}", file=sys.stderr)
    if n_kept == 0:
        raise SystemExit("ERROR: no parseable circRNA_IDs found in the BSJ/FSJ matrices.")
    print(f"[ciri3_to_dcc_format] wrote {crc}", file=sys.stderr)
    print(f"[ciri3_to_dcc_format] wrote {coord}", file=sys.stderr)
    print(f"[ciri3_to_dcc_format] wrote {lin}", file=sys.stderr)


if __name__ == "__main__":
    main()
