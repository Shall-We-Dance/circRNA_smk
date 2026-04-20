from pathlib import Path
import csv

pair = snakemake.params.pair
group_a = snakemake.params.group_a
group_b = snakemake.params.group_b
groups = snakemake.params.groups

selected_samples = list(groups[group_a]) + list(groups[group_b])
sample_to_class = {sample: group_a for sample in groups[group_a]}
sample_to_class.update({sample: group_b for sample in groups[group_b]})


def read_tsv(path: Path):
    with path.open() as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    if not rows:
        raise ValueError(f"Input TSV is empty: {path}")
    return rows[0], rows[1:]


# 1) infor.tsv
infor_out = Path(snakemake.output.infor)
infor_out.parent.mkdir(parents=True, exist_ok=True)
with infor_out.open("w", newline="") as fh:
    writer = csv.writer(fh, delimiter="\t")
    writer.writerow(["Sample", "Class"])
    for sample in selected_samples:
        writer.writerow([sample, sample_to_class[sample]])


# 2) subset matrices
for in_matrix, out_matrix in [
    (Path(snakemake.input.bsj_matrix), Path(snakemake.output.bsj_matrix)),
    (Path(snakemake.input.fsj_matrix), Path(snakemake.output.fsj_matrix)),
]:
    header, rows = read_tsv(in_matrix)
    col_idx = {name: idx for idx, name in enumerate(header)}

    missing = [s for s in selected_samples if s not in col_idx]
    if missing:
        raise ValueError(
            f"Missing selected samples in matrix header for pair '{pair}' ({in_matrix}): {missing}"
        )

    out_matrix.parent.mkdir(parents=True, exist_ok=True)
    with out_matrix.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow([header[0]] + selected_samples)
        for row in rows:
            if not row:
                continue
            writer.writerow([row[0]] + [row[col_idx[s]] if col_idx[s] < len(row) else "0" for s in selected_samples])
