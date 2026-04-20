import csv
import re
from pathlib import Path

import pysam


bsj_matrix_path = Path(snakemake.input.bsj)
ciri3_path = Path(snakemake.input.ciri3)
fasta_path = Path(snakemake.input.fasta)

site_out = Path(snakemake.output.site_table)
site_fasta_out = Path(snakemake.output.site_fasta)

flank = int(snakemake.params.flank)
weight_mode = str(snakemake.params.weight_mode).strip().lower()

if flank < 1:
    raise ValueError("motif.flank must be >= 1")
if weight_mode not in {"sum", "mean", "none"}:
    raise ValueError("motif.weight_mode must be one of: sum, mean, none")

_comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq: str) -> str:
    return seq.translate(_comp)[::-1]


def parse_bsj_id(circ_id: str):
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
    gene_map = {}

    with ciri3_file.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None:
            return strand_map, gene_map

        id_col = None
        strand_col = None
        gene_col = None

        for candidate in ("circRNA_ID", "circRNA", "circ_id", "junction"):
            if candidate in reader.fieldnames:
                id_col = candidate
                break

        for candidate in ("strand", "Strand", "bsj_strand"):
            if candidate in reader.fieldnames:
                strand_col = candidate
                break

        for candidate in ("gene_id", "gene", "gene_name"):
            if candidate in reader.fieldnames:
                gene_col = candidate
                break

        if id_col is None:
            return strand_map, gene_map

        for row in reader:
            cid = (row.get(id_col) or "").strip()
            if not cid:
                continue

            if strand_col is not None:
                strand = (row.get(strand_col) or "").strip()
                if strand in {"+", "-"} and cid not in strand_map:
                    strand_map[cid] = strand

            if gene_col is not None:
                gene_id = (row.get(gene_col) or "").strip()
                if gene_id and cid not in gene_map:
                    gene_map[cid] = gene_id

    return strand_map, gene_map


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
    return fasta.fetch(chrom, start0, end0).upper()


strand_map, gene_map = get_ciri3_strand_map(ciri3_path)
bsj_weights = read_bsj_counts(bsj_matrix_path)

site_rows = []

with pysam.FastaFile(str(fasta_path)) as fa:
    for circ_id, weight in bsj_weights.items():
        parsed = parse_bsj_id(circ_id)
        if parsed is None:
            continue

        chrom, bsj_start, bsj_end = parsed
        if chrom not in fa.references:
            continue

        strand = strand_map.get(circ_id, "+")
        gene_id = gene_map.get(circ_id, "NA")

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
                "gene_id": gene_id,
                "weight": f"{weight:.6g}",
                "donor_window_seq": donor_seq,
                "acceptor_window_seq": acceptor_seq,
            }
        )

site_out.parent.mkdir(parents=True, exist_ok=True)
with site_out.open("w", newline="") as fh:
    fieldnames = [
        "circRNA_ID",
        "chrom",
        "bsj_start",
        "bsj_end",
        "strand",
        "gene_id",
        "weight",
        "donor_window_seq",
        "acceptor_window_seq",
    ]
    writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    writer.writerows(site_rows)

site_fasta_out.parent.mkdir(parents=True, exist_ok=True)
with site_fasta_out.open("w") as fh:
    for row in site_rows:
        prefix = (
            f"{row['circRNA_ID']}|{row['chrom']}:{row['bsj_start']}-{row['bsj_end']}"
            f"|strand={row['strand']}|gene_id={row['gene_id']}|weight={row['weight']}"
        )
        fh.write(f">{prefix}|site=donor\n{row['donor_window_seq']}\n")
        fh.write(f">{prefix}|site=acceptor\n{row['acceptor_window_seq']}\n")