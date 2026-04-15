from __future__ import annotations

import csv
import re
from collections import Counter, defaultdict
from pathlib import Path

import pysam


bsj_matrix_path = Path(snakemake.input.bsj)
ciri3_merged_path = Path(snakemake.input.ciri3)
fasta_path = Path(snakemake.input.fasta)

site_out = Path(snakemake.output.site_table)
summary_out = Path(snakemake.output.motif_summary)

flank = int(snakemake.params.get("flank", 30))
kmer_size = int(snakemake.params.get("kmer", 4))
weight_mode = str(snakemake.params.get("weight_mode", "sum")).strip().lower()


if kmer_size <= 0:
    raise ValueError("motif.kmer must be a positive integer")
if flank < kmer_size:
    raise ValueError("motif.flank must be >= motif.kmer")
if weight_mode not in {"sum", "mean", "none"}:
    raise ValueError("motif.weight_mode must be one of: sum, mean, none")


_comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq: str) -> str:
    return seq.translate(_comp)[::-1]


def parse_bsj_id(circ_id: str):
    """
    Parse common CIRI-style BSJ ids, e.g.
      chr1:100|200
      chr1:100-200
      chr1:100..200
    """
    m = re.match(r"^([^:]+):(\d+)(?:\||-|\.\.)(\d+)$", circ_id)
    if not m:
        return None
    chrom, start_s, end_s = m.groups()
    start, end = int(start_s), int(end_s)
    if end < start:
        start, end = end, start
    return chrom, start, end


def get_ciri3_strand_map(ciri3_file: Path):
    strand_map = {}
    with ciri3_file.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None:
            return strand_map

        id_col = None
        strand_col = None
        for candidate in ("circRNA_ID", "circRNA", "circ_id", "junction"):
            if candidate in reader.fieldnames:
                id_col = candidate
                break
        for candidate in ("strand", "Strand", "bsj_strand"):
            if candidate in reader.fieldnames:
                strand_col = candidate
                break

        if id_col is None or strand_col is None:
            return strand_map

        for row in reader:
            cid = (row.get(id_col) or "").strip()
            strand = (row.get(strand_col) or "").strip()
            if cid and strand in {"+", "-"} and cid not in strand_map:
                strand_map[cid] = strand

    return strand_map


def read_bsj_counts(bsj_file: Path):
    bsj_weight = {}
    with bsj_file.open() as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, None)
        if header is None:
            raise ValueError("BSJ matrix is empty")
        if len(header) < 2:
            raise ValueError("BSJ matrix must include circRNA id + >=1 sample columns")

        for row in reader:
            if not row:
                continue
            circ_id = row[0].strip()
            if not circ_id:
                continue
            values = []
            for val in row[1:]:
                try:
                    values.append(float(val))
                except (TypeError, ValueError):
                    values.append(0.0)
            if weight_mode == "none":
                bsj_weight[circ_id] = 1.0
            elif weight_mode == "mean":
                bsj_weight[circ_id] = sum(values) / max(len(values), 1)
            else:
                bsj_weight[circ_id] = sum(values)
    return bsj_weight


def extract_site_sequences(fasta: pysam.FastaFile, chrom: str, pos_1based: int, flank_size: int):
    start0 = max(0, pos_1based - flank_size - 1)
    end0 = pos_1based + flank_size
    seq = fasta.fetch(chrom, start0, end0).upper()
    return seq


strand_map = get_ciri3_strand_map(ciri3_merged_path)
bsj_weights = read_bsj_counts(bsj_matrix_path)

site_rows = []
motif_total = Counter()
motif_unique = Counter()
motif_weighted = defaultdict(float)

with pysam.FastaFile(str(fasta_path)) as fa:
    for circ_id, weight in bsj_weights.items():
        parsed = parse_bsj_id(circ_id)
        if parsed is None:
            continue

        chrom, bsj_start, bsj_end = parsed
        if chrom not in fa.references:
            continue

        strand = strand_map.get(circ_id, "+")

        donor_seq = extract_site_sequences(fa, chrom, bsj_start, flank)
        acceptor_seq = extract_site_sequences(fa, chrom, bsj_end, flank)

        if strand == "-":
            donor_seq = revcomp(donor_seq)
            acceptor_seq = revcomp(acceptor_seq)

        site_rows.append(
            {
                "circRNA_ID": circ_id,
                "chrom": chrom,
                "bsj_start": bsj_start,
                "bsj_end": bsj_end,
                "strand": strand,
                "weight": f"{weight:.6g}",
                "donor_window_seq": donor_seq,
                "acceptor_window_seq": acceptor_seq,
            }
        )

        seen_once = set()
        for side, seq in (("donor", donor_seq), ("acceptor", acceptor_seq)):
            for i in range(0, len(seq) - kmer_size + 1):
                motif = seq[i : i + kmer_size]
                if "N" in motif:
                    continue
                key = (side, motif)
                motif_total[key] += 1
                if key not in seen_once:
                    motif_unique[key] += 1
                    motif_weighted[key] += weight
                    seen_once.add(key)

site_out.parent.mkdir(parents=True, exist_ok=True)
with site_out.open("w", newline="") as fh:
    fieldnames = [
        "circRNA_ID",
        "chrom",
        "bsj_start",
        "bsj_end",
        "strand",
        "weight",
        "donor_window_seq",
        "acceptor_window_seq",
    ]
    writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    writer.writerows(site_rows)

summary_out.parent.mkdir(parents=True, exist_ok=True)
with summary_out.open("w", newline="") as fh:
    writer = csv.writer(fh, delimiter="\t")
    writer.writerow(
        [
            "side",
            "motif",
            "total_occurrence",
            "unique_bsj_occurrence",
            "weighted_bsj_occurrence",
            "weight_mode",
            "kmer",
            "flank",
        ]
    )

    def sort_key(item):
        key, _ = item
        return (
            -motif_weighted.get(key, 0.0),
            -motif_unique.get(key, 0),
            -motif_total.get(key, 0),
            key[0],
            key[1],
        )

    for (side, motif), total in sorted(motif_total.items(), key=sort_key):
        writer.writerow(
            [
                side,
                motif,
                total,
                motif_unique[(side, motif)],
                f"{motif_weighted[(side, motif)]:.6g}",
                weight_mode,
                kmer_size,
                flank,
            ]
        )
