import csv
from pathlib import Path

pair = snakemake.params.pair
group_a = snakemake.params.group_a
group_b = snakemake.params.group_b
groups = snakemake.params.groups
sample_names = list(snakemake.params.sample_names)

selected_samples = list(groups[group_a]) + list(groups[group_b])
sample_to_class = {sample: group_a for sample in groups[group_a]}
sample_to_class.update({sample: group_b for sample in groups[group_b]})

input_ciri3_files = [Path(p) for p in snakemake.input.ciri3]
sample_to_ciri3 = dict(zip(sample_names, input_ciri3_files))

for sample in selected_samples:
    if sample not in sample_to_ciri3:
        raise ValueError(f"Sample '{sample}' required by pair '{pair}' not found in per-sample CIRI3 outputs")


def read_tsv(path: Path):
    with path.open() as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    if not rows:
        raise ValueError(f"Input TSV is empty: {path}")
    return rows[0], rows[1:]


# 1) infor.tsv
infor_out = Path(snakemake.output.infor)
infor_out.parent.mkdir(parents=True, exist_ok=True)

class_counts = {group_a: 0, group_b: 0}
with infor_out.open("w", newline="") as fh:
    writer = csv.writer(fh, delimiter="\t")
    writer.writerow(["Sample", "Path", "Class", "Num"])
    for sample in selected_samples:
        cls = sample_to_class[sample]
        class_counts[cls] += 1
        writer.writerow([sample, str(sample_to_ciri3[sample].resolve()), cls, class_counts[cls]])


# 2) BSJ matrix subset
bsj_header, bsj_rows = read_tsv(Path(snakemake.input.bsj_matrix))
bsj_col_idx = {name: idx for idx, name in enumerate(bsj_header)}
missing_bsj = [s for s in selected_samples if s not in bsj_col_idx]
if missing_bsj:
    raise ValueError(f"Missing selected samples in BSJ matrix header: {missing_bsj}")

bsj_out = Path(snakemake.output.bsj_matrix)
bsj_out.parent.mkdir(parents=True, exist_ok=True)
with bsj_out.open("w", newline="") as fh:
    writer = csv.writer(fh, delimiter="\t")
    writer.writerow([bsj_header[0]] + selected_samples)
    for row in bsj_rows:
        if not row:
            continue
        writer.writerow([row[0]] + [row[bsj_col_idx[s]] if bsj_col_idx[s] < len(row) else "0" for s in selected_samples])


# 3) Gene expression matrix subset (featureCounts -> Geneid + selected samples)
gene_counts = Path(snakemake.input.gene_counts)
gene_expr_out = Path(snakemake.output.gene_expression)
gene_expr_out.parent.mkdir(parents=True, exist_ok=True)

def normalize_featurecounts_colname(name: str) -> str:
    """
    featureCounts sample columns are often full BAM paths:
      results/star/<sample>/<sample>.Aligned.sortedByCoord.out.bam
    Normalize those to <sample> so they can match config sample names.
    """
    p = Path(name)
    base = p.name
    suffix = ".Aligned.sortedByCoord.out.bam"
    if base.endswith(suffix):
        return base[: -len(suffix)]
    return base


with gene_counts.open() as in_fh, gene_expr_out.open("w", newline="") as out_fh:
    reader = csv.reader((line for line in in_fh if not line.startswith("#")), delimiter="\t")
    writer = csv.writer(out_fh, delimiter="\t")

    header = next(reader, None)
    if not header:
        raise ValueError(f"FeatureCounts matrix is empty: {gene_counts}")

    col_idx = {name: idx for idx, name in enumerate(header)}
    normalized_idx = {}
    for idx, col_name in enumerate(header):
        norm = normalize_featurecounts_colname(col_name)
        # keep first match if duplicates appear after normalization
        if norm not in normalized_idx:
            normalized_idx[norm] = idx

    missing_gene = [s for s in selected_samples if s not in col_idx and s not in normalized_idx]
    if missing_gene:
        raise ValueError(f"Missing selected samples in featureCounts header: {missing_gene}")

    writer.writerow(["Geneid"] + selected_samples)

    geneid_idx = col_idx.get("Geneid", 0)
    for row in reader:
        if not row:
            continue
        geneid = row[geneid_idx]
        values = []
        for sample in selected_samples:
            sample_idx = col_idx.get(sample, normalized_idx.get(sample))
            values.append(row[sample_idx] if sample_idx is not None and sample_idx < len(row) else "0")
        writer.writerow([geneid] + values)
