from pathlib import Path
import csv

samples = list(snakemake.params.samples)
ciri3_files = list(map(Path, snakemake.input.ciri3))
bsj_files = list(map(Path, snakemake.input.bsj))
fsj_files = list(map(Path, snakemake.input.fsj))

all_ciri3_out = Path(snakemake.output.ciri3)
all_bsj_out = Path(snakemake.output.bsj)
all_fsj_out = Path(snakemake.output.fsj)


def read_tsv(path: Path):
    with path.open() as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    if not rows:
        raise ValueError(f"{path} is empty")
    return rows[0], rows[1:]


def pick_value_column(header, sample_name):
    if sample_name in header:
        return header.index(sample_name)
    if len(header) >= 2:
        return 1
    raise ValueError(f"No value column found in matrix header: {header}")


def merge_matrix(matrix_files, sample_names, out_path):
    index_name = None
    merged = {}

    for sample_name, matrix_file in zip(sample_names, matrix_files):
        header, rows = read_tsv(matrix_file)
        if index_name is None:
            index_name = header[0]
        value_col = pick_value_column(header, sample_name)

        for row in rows:
            if not row or len(row) <= value_col:
                continue
            circ_id = row[0]
            value = row[value_col]
            if circ_id not in merged:
                merged[circ_id] = {}
            merged[circ_id][sample_name] = value

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow([index_name] + sample_names)
        for circ_id in sorted(merged):
            writer.writerow([circ_id] + [merged[circ_id].get(sample, "0") for sample in sample_names])


# Merge per-sample .ciri3 result tables by appending a sample column.
all_ciri3_out.parent.mkdir(parents=True, exist_ok=True)
with all_ciri3_out.open("w", newline="") as out_fh:
    writer = csv.writer(out_fh, delimiter="\t")
    wrote_header = False

    for sample_name, ciri3_file in zip(samples, ciri3_files):
        header, rows = read_tsv(ciri3_file)
        if not wrote_header:
            writer.writerow(header + ["sample"])
            wrote_header = True
        for row in rows:
            writer.writerow(row + [sample_name])

merge_matrix(bsj_files, samples, all_bsj_out)
merge_matrix(fsj_files, samples, all_fsj_out)
