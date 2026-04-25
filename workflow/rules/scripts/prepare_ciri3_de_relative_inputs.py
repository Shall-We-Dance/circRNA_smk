import csv
import re
from pathlib import Path


comparison = snakemake.wildcards.comparison
selected_samples = list(snakemake.params.selected_samples)
sample_classes = dict(snakemake.params.sample_classes)

input_ciri3_files = [Path(p) for p in snakemake.input.ciri3]
sample_to_ciri3 = dict(zip(selected_samples, input_ciri3_files))

info_out = Path(snakemake.output.info)
bsj_matrix_out = Path(snakemake.output.bsj_matrix)
circ_gene_out = Path(snakemake.output.circ_gene)
log_path = Path(snakemake.log[0]) if snakemake.log else None


def log_message(message):
    if log_path is None:
        return
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a") as fh:
        fh.write(message + "\n")


def read_tsv(path):
    with Path(path).open() as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    if not rows:
        raise ValueError(f"Input TSV is empty: {path}")
    return rows[0], rows[1:]


missing_ciri3 = [sample for sample in selected_samples if sample not in sample_to_ciri3]
if missing_ciri3:
    raise ValueError(
        "Selected samples missing from per-sample CIRI3 inputs for comparison "
        f"{comparison}: " + ", ".join(missing_ciri3)
    )

info_out.parent.mkdir(parents=True, exist_ok=True)
with info_out.open("w", newline="") as fh:
    writer = csv.writer(fh, delimiter="\t")
    writer.writerow(["Sample", "Path", "Class"])
    for sample in selected_samples:
        writer.writerow([
            sample,
            str(sample_to_ciri3[sample].resolve()),
            sample_classes[sample],
        ])

bsj_header, bsj_rows = read_tsv(snakemake.input.bsj)
missing_bsj = [sample for sample in selected_samples if sample not in bsj_header]
if missing_bsj:
    raise ValueError(
        "Selected samples missing from merged BSJ matrix for comparison "
        f"{comparison}: " + ", ".join(missing_bsj)
    )

keep_idx = [0] + [bsj_header.index(sample) for sample in selected_samples]
bsj_matrix_out.parent.mkdir(parents=True, exist_ok=True)
with bsj_matrix_out.open("w", newline="") as fh:
    writer = csv.writer(fh, delimiter="\t")
    writer.writerow([bsj_header[i] for i in keep_idx])
    for row in bsj_rows:
        padded = row + [""] * (len(bsj_header) - len(row))
        writer.writerow([padded[i] if padded[i] != "" else "0" for i in keep_idx])

pairs = set()
with Path(snakemake.input.ciri3_all).open() as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    if reader.fieldnames is None:
        raise ValueError(f"Empty or malformed CIRI3 table: {snakemake.input.ciri3_all}")

    if "circRNA_ID" not in reader.fieldnames or "gene_id" not in reader.fieldnames:
        raise ValueError("all_samples.ciri3 must contain 'circRNA_ID' and 'gene_id' columns.")

    for row in reader:
        circ_id = row["circRNA_ID"].strip()
        gene_field = row["gene_id"].strip()

        if not circ_id or not gene_field or gene_field == "NA":
            continue

        gene_parts = [
            gene.strip()
            for gene in re.split(r"[;,]", gene_field)
            if gene.strip() and gene.strip() != "NA"
        ]
        for gene_id in gene_parts:
            pairs.add((circ_id, gene_id))

circ_gene_out.parent.mkdir(parents=True, exist_ok=True)
with circ_gene_out.open("w", newline="") as fh:
    writer = csv.writer(fh, delimiter="\t")
    writer.writerow(["circRNA_ID", "gene_id"])
    for circ_id, gene_id in sorted(pairs):
        writer.writerow([circ_id, gene_id])

log_message(
    f"Prepared CIRI3 DE_Relative inputs for {comparison}: "
    f"{len(selected_samples)} samples, {len(bsj_rows)} BSJ rows, {len(pairs)} circ-gene pairs."
)
