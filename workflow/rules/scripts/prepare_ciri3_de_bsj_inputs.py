import csv
import keyword
import re
from pathlib import Path


comparison = snakemake.wildcards.comparison
selected_samples = list(snakemake.params.selected_samples)
sample_classes = dict(snakemake.params.sample_classes)

input_ciri3_files = [Path(p) for p in snakemake.input.ciri3]
sample_to_ciri3 = dict(zip(selected_samples, input_ciri3_files))

info_out = Path(snakemake.output.info)
gene_expression_out = Path(snakemake.output.gene_expression)
bsj_matrix_out = Path(snakemake.output.bsj_matrix)
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


def normalize_featurecounts_colname(name):
    base = Path(name).name
    suffix = ".Aligned.sortedByCoord.out.bam"
    if base.endswith(suffix):
        return base[:-len(suffix)]
    return base


R_RESERVED_WORDS = {
    "if", "else", "repeat", "while", "function", "for", "in", "next", "break",
    "TRUE", "FALSE", "NULL", "Inf", "NaN", "NA", "NA_integer_", "NA_real_",
    "NA_complex_", "NA_character_",
}


def r_safe_name(name):
    safe = re.sub(r"[^A-Za-z0-9._]", ".", str(name))
    if not safe:
        safe = "X"
    if safe[0].isdigit() or (safe[0] == "." and len(safe) > 1 and safe[1].isdigit()):
        safe = "X" + safe
    if safe in R_RESERVED_WORDS or keyword.iskeyword(safe):
        safe += "."
    return safe


def unique_names(names):
    used = set()
    counts = {}
    out = []
    for name in names:
        base = r_safe_name(name)
        candidate = base
        while candidate in used:
            counts[base] = counts.get(base, 0) + 1
            candidate = f"{base}.{counts[base]}"
        used.add(candidate)
        counts.setdefault(base, 0)
        out.append(candidate)
    return out


def sample_output_name_map(samples):
    aliases = unique_names(samples)
    return dict(zip(samples, aliases))


sample_output_names = sample_output_name_map(selected_samples)


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
            sample_output_names[sample],
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
    writer.writerow([bsj_header[0]] + [sample_output_names[sample] for sample in selected_samples])
    for row in bsj_rows:
        padded = row + [""] * (len(bsj_header) - len(row))
        writer.writerow([padded[i] if padded[i] != "" else "0" for i in keep_idx])

with Path(snakemake.input.featurecounts).open() as fh:
    raw_lines = [line.rstrip("\n") for line in fh if line.strip()]

fc_lines = [line for line in raw_lines if not line.startswith("#")]
if len(fc_lines) < 2:
    raise ValueError(f"featureCounts file appears malformed: {snakemake.input.featurecounts}")

fc_rows = list(csv.reader(fc_lines, delimiter="\t"))
fc_header = fc_rows[0]

if "Geneid" not in fc_header:
    raise ValueError("featureCounts header does not contain 'Geneid'.")

geneid_idx = fc_header.index("Geneid")
sample_to_fc_idx = {}
count_start_idx = 6
for idx, col_name in enumerate(fc_header[count_start_idx:], start=count_start_idx):
    sample_name = normalize_featurecounts_colname(col_name)
    if sample_name not in sample_to_fc_idx:
        sample_to_fc_idx[sample_name] = idx

missing_fc = [sample for sample in selected_samples if sample not in sample_to_fc_idx]
if missing_fc:
    raise ValueError(
        "Selected samples missing from featureCounts header for comparison "
        f"{comparison}: " + ", ".join(missing_fc)
    )

gene_expression_out.parent.mkdir(parents=True, exist_ok=True)
with gene_expression_out.open("w", newline="") as fh:
    writer = csv.writer(fh, delimiter="\t")
    writer.writerow(["Geneid"] + [sample_output_names[sample] for sample in selected_samples])
    for row in fc_rows[1:]:
        padded = row + ["0"] * (len(fc_header) - len(row))
        out_row = [padded[geneid_idx]]
        for sample in selected_samples:
            value = padded[sample_to_fc_idx[sample]]
            out_row.append(value if value != "" else "0")
        writer.writerow(out_row)

log_message(
    f"Prepared CIRI3 DE_BSJ inputs for {comparison}: "
    f"{len(selected_samples)} samples, {len(bsj_rows)} BSJ rows, {len(fc_rows) - 1} genes."
)
renamed_samples = [
    f"{sample}->{sample_output_names[sample]}"
    for sample in selected_samples
    if sample != sample_output_names[sample]
]
if renamed_samples:
    log_message(
        "Normalized CIRI3 DE sample names for R compatibility: "
        + ", ".join(renamed_samples)
    )
