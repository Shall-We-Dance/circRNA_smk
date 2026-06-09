import csv
import keyword
import re
from pathlib import Path


comparison = snakemake.wildcards.comparison
selected_samples = list(snakemake.params.selected_samples)
sample_classes = dict(snakemake.params.sample_classes)
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


info_out = Path(snakemake.output.info)
info_out.parent.mkdir(parents=True, exist_ok=True)
with info_out.open("w", newline="") as fh:
    writer = csv.writer(fh, delimiter="\t")
    writer.writerow(["Sample", "Class"])
    for sample in selected_samples:
        writer.writerow([sample_output_names[sample], sample_classes[sample]])

for in_matrix, out_matrix in [
    (snakemake.input.bsj, snakemake.output.bsj_matrix),
    (snakemake.input.fsj, snakemake.output.fsj_matrix),
]:
    header, rows = read_tsv(in_matrix)
    missing = [sample for sample in selected_samples if sample not in header]
    if missing:
        raise ValueError(
            f"Selected samples missing from matrix {in_matrix} for comparison "
            f"{comparison}: " + ", ".join(missing)
        )

    keep_idx = [0] + [header.index(sample) for sample in selected_samples]
    out_path = Path(out_matrix)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow([header[0]] + [sample_output_names[sample] for sample in selected_samples])
        for row in rows:
            padded = row + [""] * (len(header) - len(row))
            writer.writerow([padded[i] if padded[i] != "" else "0" for i in keep_idx])

log_message(f"Prepared CIRI3 DE_Ratio inputs for {comparison}: {len(selected_samples)} samples.")
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
