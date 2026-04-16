import csv
import math
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ciri3_path = Path(snakemake.input.ciri3)
bsj_path = Path(snakemake.input.bsj)
fsj_path = Path(snakemake.input.fsj)
sample = snakemake.wildcards.sample

circ_table_out = Path(snakemake.output.circ_table)
summary_out = Path(snakemake.output.summary)
dist_out = Path(snakemake.output.dist)
plot_out = Path(snakemake.output.plot)


def parse_circ_id(circ_id: str):
    match = re.match(r"^([^:]+):(\d+)(?:\||-|\.\.)(\d+)$", circ_id)
    if not match:
        return None
    chrom, start_s, end_s = match.groups()
    start = int(start_s)
    end = int(end_s)
    if end < start:
        start, end = end, start
    return chrom, start, end


def pick_value_column(header, sample_name):
    if sample_name in header:
        return header.index(sample_name)
    if len(header) >= 2:
        return 1
    raise ValueError(f"No sample value column detected in matrix header: {header}")


def read_matrix_values(path: Path, sample_name: str):
    values = {}
    with path.open() as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, None)
        if header is None:
            raise ValueError(f"{path} is empty")
        value_col = pick_value_column(header, sample_name)
        for row in reader:
            if not row:
                continue
            circ_id = row[0].strip()
            if not circ_id:
                continue
            raw = row[value_col].strip() if len(row) > value_col else "0"
            try:
                values[circ_id] = float(raw)
            except ValueError:
                values[circ_id] = 0.0
    return values


def load_ciri3_meta(path: Path):
    strand_map = {}
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None:
            return strand_map

        circ_col = None
        strand_col = None
        for candidate in ("circRNA_ID", "circRNA", "circ_id", "junction"):
            if candidate in reader.fieldnames:
                circ_col = candidate
                break
        for candidate in ("strand", "Strand", "bsj_strand"):
            if candidate in reader.fieldnames:
                strand_col = candidate
                break

        if circ_col is None:
            return strand_map

        for row in reader:
            circ_id = (row.get(circ_col) or "").strip()
            if not circ_id:
                continue
            strand = (row.get(strand_col) or ".").strip() if strand_col else "."
            if strand not in {"+", "-"}:
                strand = "."
            if circ_id not in strand_map:
                strand_map[circ_id] = strand
    return strand_map


def value_bin(v, edges):
    for i in range(len(edges) - 1):
        left = edges[i]
        right = edges[i + 1]
        if left <= v < right:
            return f"[{left},{right})"
    return f">={edges[-1]}"


bsj_edges = [0, 1, 2, 5, 10, 20, 50]
fsj_edges = [0, 1, 2, 5, 10, 20, 50]
ratio_edges = [0, 0.25, 0.5, 1, 2, 5, 10]
span_edges = [0, 200, 500, 1000, 2000, 5000, 10000]

bsj_values = read_matrix_values(bsj_path, sample)
fsj_values = read_matrix_values(fsj_path, sample)
strand_map = load_ciri3_meta(ciri3_path)

all_ids = sorted(set(bsj_values) | set(fsj_values))

rows = []
for circ_id in all_ids:
    parsed = parse_circ_id(circ_id)
    if parsed is None:
        continue
    chrom, start, end = parsed
    span = end - start + 1
    bsj = bsj_values.get(circ_id, 0.0)
    fsj = fsj_values.get(circ_id, 0.0)
    ratio = (bsj + 1.0) / (fsj + 1.0)
    rows.append(
        {
            "sample": sample,
            "circRNA_ID": circ_id,
            "chrom": chrom,
            "bsj_start": start,
            "bsj_end": end,
            "circ_span": span,
            "strand": strand_map.get(circ_id, "."),
            "bsj_count": f"{bsj:.6g}",
            "fsj_count": f"{fsj:.6g}",
            "bsj_fsj_ratio": f"{ratio:.6g}",
            "log2_bsj_fsj_ratio": f"{math.log2(ratio):.6g}",
            "bsj_bin": value_bin(bsj, bsj_edges),
            "fsj_bin": value_bin(fsj, fsj_edges),
            "ratio_bin": value_bin(ratio, ratio_edges),
            "span_bin": value_bin(span, span_edges),
        }
    )

circ_table_out.parent.mkdir(parents=True, exist_ok=True)
with circ_table_out.open("w", newline="") as fh:
    fieldnames = [
        "sample",
        "circRNA_ID",
        "chrom",
        "bsj_start",
        "bsj_end",
        "circ_span",
        "strand",
        "bsj_count",
        "fsj_count",
        "bsj_fsj_ratio",
        "log2_bsj_fsj_ratio",
        "bsj_bin",
        "fsj_bin",
        "ratio_bin",
        "span_bin",
    ]
    writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)

metrics = []
bsj_list = [float(r["bsj_count"]) for r in rows]
fsj_list = [float(r["fsj_count"]) for r in rows]
ratio_list = [float(r["bsj_fsj_ratio"]) for r in rows]
span_list = [int(r["circ_span"]) for r in rows]

if rows:
    strand_plus = sum(1 for r in rows if r["strand"] == "+")
    strand_minus = sum(1 for r in rows if r["strand"] == "-")
    metrics.extend(
        [
            ("sample", sample),
            ("n_circRNAs", len(rows)),
            ("n_bsj_gt_0", sum(v > 0 for v in bsj_list)),
            ("n_fsj_gt_0", sum(v > 0 for v in fsj_list)),
            ("median_bsj", sorted(bsj_list)[len(bsj_list) // 2]),
            ("median_fsj", sorted(fsj_list)[len(fsj_list) // 2]),
            ("median_ratio", sorted(ratio_list)[len(ratio_list) // 2]),
            ("median_span", sorted(span_list)[len(span_list) // 2]),
            ("strand_plus_fraction", strand_plus / len(rows)),
            ("strand_minus_fraction", strand_minus / len(rows)),
        ]
    )
else:
    metrics.extend(
        [
            ("sample", sample),
            ("n_circRNAs", 0),
            ("n_bsj_gt_0", 0),
            ("n_fsj_gt_0", 0),
            ("median_bsj", 0),
            ("median_fsj", 0),
            ("median_ratio", 0),
            ("median_span", 0),
            ("strand_plus_fraction", 0),
            ("strand_minus_fraction", 0),
        ]
    )

with summary_out.open("w", newline="") as fh:
    writer = csv.writer(fh, delimiter="\t")
    writer.writerow(["sample", "metric", "value"])
    for metric, value in metrics:
        writer.writerow([sample, metric, value])


def count_bins(column, ordered_bins):
    c = {b: 0 for b in ordered_bins}
    for row in rows:
        b = row[column]
        c[b] = c.get(b, 0) + 1
    return c


bsj_bins = [f"[{bsj_edges[i]},{bsj_edges[i+1]})" for i in range(len(bsj_edges) - 1)] + [f">={bsj_edges[-1]}"]
fsj_bins = [f"[{fsj_edges[i]},{fsj_edges[i+1]})" for i in range(len(fsj_edges) - 1)] + [f">={fsj_edges[-1]}"]
ratio_bins = [f"[{ratio_edges[i]},{ratio_edges[i+1]})" for i in range(len(ratio_edges) - 1)] + [f">={ratio_edges[-1]}"]
span_bins = [f"[{span_edges[i]},{span_edges[i+1]})" for i in range(len(span_edges) - 1)] + [f">={span_edges[-1]}"]

bin_specs = [
    ("bsj_count", "bsj_bin", bsj_bins),
    ("fsj_count", "fsj_bin", fsj_bins),
    ("bsj_fsj_ratio", "ratio_bin", ratio_bins),
    ("circ_span", "span_bin", span_bins),
]

with dist_out.open("w", newline="") as fh:
    writer = csv.writer(fh, delimiter="\t")
    writer.writerow(["sample", "feature", "bin", "count", "fraction"])
    denom = max(len(rows), 1)
    for feature, col, ordered in bin_specs:
        counts = count_bins(col, ordered)
        for b in ordered:
            count = counts.get(b, 0)
            writer.writerow([sample, feature, b, count, count / denom])

fig, axes = plt.subplots(2, 2, figsize=(12, 8))
for ax, (feature, col, ordered) in zip(axes.flatten(), bin_specs):
    counts = count_bins(col, ordered)
    xs = list(range(len(ordered)))
    ys = [counts.get(b, 0) for b in ordered]
    ax.bar(xs, ys, color="#4C72B0")
    ax.set_xticks(xs)
    ax.set_xticklabels(ordered, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("circRNA count")
    ax.set_title(f"{sample} | {feature}")

fig.tight_layout()
plot_out.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(plot_out, dpi=150)
plt.close(fig)
