import csv
import math
import re
from collections import defaultdict
from pathlib import Path

import matplotlib
import pysam

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ciri3_path = Path(snakemake.input.ciri3)
bsj_path = Path(snakemake.input.bsj)
fsj_path = Path(snakemake.input.fsj)
fasta_path = Path(snakemake.input.fasta)
sample = snakemake.wildcards.sample

circ_table_out = Path(snakemake.output.circ_table)
summary_out = Path(snakemake.output.summary)
dist_out = Path(snakemake.output.dist)
abs_out = Path(snakemake.output.abs)
plot_out = Path(snakemake.output.plot)


_comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq: str):
    return seq.translate(_comp)[::-1]


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


def safe_float(value, default=0.0):
    try:
        if value is None:
            return default
        parsed = float(str(value).strip())
        if math.isnan(parsed):
            return default
        return parsed
    except (TypeError, ValueError):
        return default


def format_float(value):
    return f"{safe_float(value):.6g}"


def format_optional_float(value):
    try:
        parsed = float(str(value).strip())
    except (TypeError, ValueError):
        return "NA"
    if math.isnan(parsed):
        return "NA"
    return f"{parsed:.6g}"


def median(values, default=0):
    if not values:
        return default
    values = sorted(values)
    mid = len(values) // 2
    if len(values) % 2:
        return values[mid]
    return (values[mid - 1] + values[mid]) / 2


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
            values[circ_id] = safe_float(raw)
    return values


def first_existing(fieldnames, candidates):
    for candidate in candidates:
        if candidate in fieldnames:
            return candidate
    return None


def load_ciri3_meta(path: Path):
    meta = {}
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None:
            return meta

        fieldnames = reader.fieldnames
        circ_col = first_existing(
            fieldnames,
            ("circRNA_ID", "circRNA", "circ_id", "junction"),
        )
        if circ_col is None:
            return meta

        strand_col = first_existing(fieldnames, ("strand", "Strand", "bsj_strand"))
        gene_col = first_existing(fieldnames, ("gene_id", "geneID", "gene", "Gene"))
        type_col = first_existing(
            fieldnames,
            ("circRNA_type", "circ_type", "type", "Type"),
        )
        sm_ms_sms_col = first_existing(fieldnames, ("SM_MS_SMS", "sm_ms_sms"))
        reported_ratio_col = first_existing(
            fieldnames,
            (
                "junction_reads_ratio",
                "junction_ratio",
                "junc_ratio",
                "circRNA_junction_ratio",
            ),
        )
        high_quality_col = first_existing(
            fieldnames,
            (
                "high_quality_bsj_reads",
                "high_quality_BSJ_reads",
                "#high_quality_BSJ_reads",
                "high_quality_junction_reads",
                "#high_quality_junction_reads",
                "high_confidence_bsj_reads",
                "high_confidence_BSJ_reads",
                "#high_confidence_BSJ_reads",
            ),
        )
        if high_quality_col is None and len(fieldnames) >= 13:
            high_quality_col = fieldnames[12]

        for row in reader:
            circ_id = (row.get(circ_col) or "").strip()
            if not circ_id or circ_id in meta:
                continue

            strand = (row.get(strand_col) or ".").strip() if strand_col else "."
            if strand not in {"+", "-"}:
                strand = "."

            gene_id = (row.get(gene_col) or "NA").strip() if gene_col else "NA"
            circ_type = (row.get(type_col) or "NA").strip() if type_col else "NA"
            sm_ms_sms = (row.get(sm_ms_sms_col) or "NA").strip() if sm_ms_sms_col else "NA"
            high_quality = (
                (row.get(high_quality_col) or "NA").strip() if high_quality_col else "NA"
            )
            reported_ratio = (
                format_optional_float(row.get(reported_ratio_col))
                if reported_ratio_col
                else "NA"
            )

            meta[circ_id] = {
                "strand": strand,
                "gene_id": gene_id or "NA",
                "circRNA_type": circ_type or "NA",
                "sm_ms_sms": sm_ms_sms or "NA",
                "high_quality_bsj_reads": high_quality or "NA",
                "ciri3_reported_junction_ratio": reported_ratio,
            }
    return meta


def value_bin(v, edges):
    for i in range(len(edges) - 1):
        left = edges[i]
        right = edges[i + 1]
        if left <= v < right:
            return f"[{left},{right})"
    return f">={edges[-1]}"


def fetch_2nt(fa: pysam.FastaFile, chrom: str, pos_1based: int):
    start0 = max(0, pos_1based - 1)
    end0 = start0 + 2
    return fa.fetch(chrom, start0, end0).upper()


def splice_dinuc(fa: pysam.FastaFile, chrom: str, start: int, end: int, strand: str):
    try:
        if strand == "-":
            donor = revcomp(fetch_2nt(fa, chrom, end - 1))
            acceptor = revcomp(fetch_2nt(fa, chrom, start))
        else:
            donor = fetch_2nt(fa, chrom, start)
            acceptor = fetch_2nt(fa, chrom, end - 1)
    except Exception:
        return "NN", "NN"

    if len(donor) != 2:
        donor = "NN"
    if len(acceptor) != 2:
        acceptor = "NN"
    return donor, acceptor


def splice_class(donor: str, acceptor: str):
    pair = f"{donor}-{acceptor}"
    if pair == "GT-AG":
        return "canonical_GU-AG"
    if pair == "GC-AG":
        return "semi_canonical_GC-AG"
    if pair == "AT-AC":
        return "minor_canonical_AU-AC"
    if "N" in donor or "N" in acceptor:
        return "unknown"
    return "non_canonical"


def dna_to_rna(seq: str):
    return seq.replace("T", "U")


def junction_ratio(bsj, fsj):
    denom = 2.0 * bsj + fsj
    if denom <= 0:
        return 0.0
    return (2.0 * bsj) / denom


def assign_abs_events(rows):
    for row in rows:
        row.update(
            {
                "abs_class": "no_abs",
                "is_a5bs": "no",
                "is_a3bs": "no",
                "a5bs_event_id": "NA",
                "a5bs_site_count": "0",
                "a5bs_event_bsj_total": "0",
                "a5bs_pcu": "0",
                "a5bs_rank_by_bsj": "NA",
                "a3bs_event_id": "NA",
                "a3bs_site_count": "0",
                "a3bs_event_bsj_total": "0",
                "a3bs_pcu": "0",
                "a3bs_rank_by_bsj": "NA",
            }
        )

    def build_groups(shared_col, variable_col):
        groups = defaultdict(list)
        for row in rows:
            if safe_float(row["bsj_count"]) <= 0:
                continue
            key = (row["chrom"], row["strand"], row[shared_col])
            groups[key].append(row)
        return {
            key: members
            for key, members in groups.items()
            if len({member[variable_col] for member in members}) >= 2
        }

    event_specs = [
        ("A5BS", "a5bs", "bsj_3p_site", "bsj_5p_site", "shared_3p"),
        ("A3BS", "a3bs", "bsj_5p_site", "bsj_3p_site", "shared_5p"),
    ]

    for event_type, prefix, shared_col, variable_col, shared_label in event_specs:
        groups = build_groups(shared_col, variable_col)
        for (chrom, strand, shared_site), members in groups.items():
            site_count = len({member[variable_col] for member in members})
            total_bsj = sum(safe_float(member["bsj_count"]) for member in members)
            event_id = f"{event_type}|{chrom}|{strand}|{shared_label}:{shared_site}"
            ranked = sorted(
                members,
                key=lambda member: (-safe_float(member["bsj_count"]), member["circRNA_ID"]),
            )
            for rank, row in enumerate(ranked, start=1):
                bsj = safe_float(row["bsj_count"])
                row[f"is_{prefix}"] = "yes"
                row[f"{prefix}_event_id"] = event_id
                row[f"{prefix}_site_count"] = str(site_count)
                row[f"{prefix}_event_bsj_total"] = format_float(total_bsj)
                row[f"{prefix}_pcu"] = format_float(bsj / total_bsj if total_bsj else 0)
                row[f"{prefix}_rank_by_bsj"] = str(rank)

    for row in rows:
        is_a5bs = row["is_a5bs"] == "yes"
        is_a3bs = row["is_a3bs"] == "yes"
        if is_a5bs and is_a3bs:
            row["abs_class"] = "A5BS_and_A3BS"
        elif is_a5bs:
            row["abs_class"] = "A5BS_only"
        elif is_a3bs:
            row["abs_class"] = "A3BS_only"


def make_abs_rows(rows):
    abs_rows = []
    for row in rows:
        for event_type, prefix, shared_col, variable_col in (
            ("A5BS", "a5bs", "bsj_3p_site", "bsj_5p_site"),
            ("A3BS", "a3bs", "bsj_5p_site", "bsj_3p_site"),
        ):
            if row[f"is_{prefix}"] != "yes":
                continue
            abs_rows.append(
                {
                    "sample": row["sample"],
                    "event_type": event_type,
                    "event_id": row[f"{prefix}_event_id"],
                    "circRNA_ID": row["circRNA_ID"],
                    "chrom": row["chrom"],
                    "strand": row["strand"],
                    "gene_id": row["gene_id"],
                    "circRNA_type": row["circRNA_type"],
                    "bsj_start": row["bsj_start"],
                    "bsj_end": row["bsj_end"],
                    "bsj_5p_site": row["bsj_5p_site"],
                    "bsj_3p_site": row["bsj_3p_site"],
                    "shared_back_splice_site": row[shared_col],
                    "alternative_back_splice_site": row[variable_col],
                    "bsj_count": row["bsj_count"],
                    "fsj_count": row["fsj_count"],
                    "junction_ratio": row["junction_ratio"],
                    "event_bsj_total": row[f"{prefix}_event_bsj_total"],
                    "site_count": row[f"{prefix}_site_count"],
                    "pcu": row[f"{prefix}_pcu"],
                    "rank_by_bsj": row[f"{prefix}_rank_by_bsj"],
                    "splice_site_class": row["splice_site_class"],
                }
            )
    return abs_rows


bsj_edges = [0, 1, 2, 5, 10, 20, 50]
fsj_edges = [0, 1, 2, 5, 10, 20, 50]
ratio_edges = [0, 0.25, 0.5, 1, 2, 5, 10]
junction_ratio_edges = [0, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 1]
span_edges = [0, 200, 500, 1000, 2000, 5000, 10000]

bsj_values = read_matrix_values(bsj_path, sample)
fsj_values = read_matrix_values(fsj_path, sample)
ciri3_meta = load_ciri3_meta(ciri3_path)

all_ids = sorted(set(bsj_values) | set(fsj_values))

rows = []
with pysam.FastaFile(str(fasta_path)) as fa:
    for circ_id in all_ids:
        parsed = parse_circ_id(circ_id)
        if parsed is None:
            continue
        chrom, start, end = parsed
        span = end - start + 1
        bsj = bsj_values.get(circ_id, 0.0)
        fsj = fsj_values.get(circ_id, 0.0)
        ratio = (bsj + 1.0) / (fsj + 1.0)
        junc_ratio = junction_ratio(bsj, fsj)
        meta = ciri3_meta.get(circ_id, {})
        strand = meta.get("strand", ".")
        bsj_5p_site = start if strand != "-" else end
        bsj_3p_site = end if strand != "-" else start
        donor, acceptor = splice_dinuc(fa, chrom, start, end, strand)
        pair_dna = f"{donor}-{acceptor}"
        pair_rna = f"{dna_to_rna(donor)}-{dna_to_rna(acceptor)}"
        rows.append(
            {
                "sample": sample,
                "circRNA_ID": circ_id,
                "chrom": chrom,
                "bsj_start": start,
                "bsj_end": end,
                "bsj_5p_site": bsj_5p_site,
                "bsj_3p_site": bsj_3p_site,
                "circ_span": span,
                "strand": strand,
                "circRNA_type": meta.get("circRNA_type", "NA"),
                "gene_id": meta.get("gene_id", "NA"),
                "sm_ms_sms": meta.get("sm_ms_sms", "NA"),
                "high_quality_bsj_reads": meta.get("high_quality_bsj_reads", "NA"),
                "bsj_count": format_float(bsj),
                "fsj_count": format_float(fsj),
                "bsj_fsj_ratio": format_float(ratio),
                "log2_bsj_fsj_ratio": format_float(math.log2(ratio)),
                "junction_ratio": format_float(junc_ratio),
                "ciri3_reported_junction_ratio": meta.get(
                    "ciri3_reported_junction_ratio", "NA"
                ),
                "donor_dinuc_dna": donor,
                "acceptor_dinuc_dna": acceptor,
                "splice_site_pair_dna": pair_dna,
                "splice_site_pair_rna": pair_rna,
                "splice_site_class": splice_class(donor, acceptor),
                "bsj_bin": value_bin(bsj, bsj_edges),
                "fsj_bin": value_bin(fsj, fsj_edges),
                "ratio_bin": value_bin(ratio, ratio_edges),
                "junction_ratio_bin": value_bin(junc_ratio, junction_ratio_edges),
                "span_bin": value_bin(span, span_edges),
            }
        )

assign_abs_events(rows)
abs_rows = make_abs_rows(rows)

circ_table_out.parent.mkdir(parents=True, exist_ok=True)
with circ_table_out.open("w", newline="") as fh:
    fieldnames = [
        "sample",
        "circRNA_ID",
        "chrom",
        "bsj_start",
        "bsj_end",
        "bsj_5p_site",
        "bsj_3p_site",
        "circ_span",
        "strand",
        "circRNA_type",
        "gene_id",
        "sm_ms_sms",
        "high_quality_bsj_reads",
        "bsj_count",
        "fsj_count",
        "bsj_fsj_ratio",
        "log2_bsj_fsj_ratio",
        "junction_ratio",
        "ciri3_reported_junction_ratio",
        "donor_dinuc_dna",
        "acceptor_dinuc_dna",
        "splice_site_pair_dna",
        "splice_site_pair_rna",
        "splice_site_class",
        "abs_class",
        "is_a5bs",
        "is_a3bs",
        "a5bs_event_id",
        "a5bs_site_count",
        "a5bs_event_bsj_total",
        "a5bs_pcu",
        "a5bs_rank_by_bsj",
        "a3bs_event_id",
        "a3bs_site_count",
        "a3bs_event_bsj_total",
        "a3bs_pcu",
        "a3bs_rank_by_bsj",
        "bsj_bin",
        "fsj_bin",
        "ratio_bin",
        "junction_ratio_bin",
        "span_bin",
    ]
    writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)

abs_out.parent.mkdir(parents=True, exist_ok=True)
with abs_out.open("w", newline="") as fh:
    fieldnames = [
        "sample",
        "event_type",
        "event_id",
        "circRNA_ID",
        "chrom",
        "strand",
        "gene_id",
        "circRNA_type",
        "bsj_start",
        "bsj_end",
        "bsj_5p_site",
        "bsj_3p_site",
        "shared_back_splice_site",
        "alternative_back_splice_site",
        "bsj_count",
        "fsj_count",
        "junction_ratio",
        "event_bsj_total",
        "site_count",
        "pcu",
        "rank_by_bsj",
        "splice_site_class",
    ]
    writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    writer.writerows(abs_rows)

bsj_list = [safe_float(r["bsj_count"]) for r in rows]
fsj_list = [safe_float(r["fsj_count"]) for r in rows]
ratio_list = [safe_float(r["bsj_fsj_ratio"]) for r in rows]
junction_ratio_list = [safe_float(r["junction_ratio"]) for r in rows]
span_list = [int(r["circ_span"]) for r in rows]

strand_plus = sum(1 for r in rows if r["strand"] == "+")
strand_minus = sum(1 for r in rows if r["strand"] == "-")
canonical = sum(1 for r in rows if r["splice_site_class"] == "canonical_GU-AG")
non_canonical = sum(1 for r in rows if r["splice_site_class"] == "non_canonical")
a5bs_events = {r["event_id"] for r in abs_rows if r["event_type"] == "A5BS"}
a3bs_events = {r["event_id"] for r in abs_rows if r["event_type"] == "A3BS"}
a5bs_circs = {r["circRNA_ID"] for r in rows if r["is_a5bs"] == "yes"}
a3bs_circs = {r["circRNA_ID"] for r in rows if r["is_a3bs"] == "yes"}
abs_circs = a5bs_circs | a3bs_circs
a5bs_pcus = [safe_float(r["pcu"]) for r in abs_rows if r["event_type"] == "A5BS"]
a3bs_pcus = [safe_float(r["pcu"]) for r in abs_rows if r["event_type"] == "A3BS"]
denom = max(len(rows), 1)

metrics = [
    ("sample", sample),
    ("n_circRNAs", len(rows)),
    ("n_bsj_gt_0", sum(v > 0 for v in bsj_list)),
    ("n_fsj_gt_0", sum(v > 0 for v in fsj_list)),
    ("total_bsj_reads", sum(bsj_list)),
    ("total_fsj_reads", sum(fsj_list)),
    ("median_bsj", median(bsj_list)),
    ("median_fsj", median(fsj_list)),
    ("median_ratio", median(ratio_list)),
    ("median_junction_ratio", median(junction_ratio_list)),
    ("n_junction_ratio_ge_0_5", sum(v >= 0.5 for v in junction_ratio_list)),
    ("median_span", median(span_list)),
    ("strand_plus_fraction", strand_plus / denom),
    ("strand_minus_fraction", strand_minus / denom),
    ("canonical_GU_AG_fraction", canonical / denom),
    ("non_canonical_fraction", non_canonical / denom),
    ("n_a5bs_events", len(a5bs_events)),
    ("n_a3bs_events", len(a3bs_events)),
    ("n_abs_events", len(a5bs_events) + len(a3bs_events)),
    ("n_a5bs_circRNAs", len(a5bs_circs)),
    ("n_a3bs_circRNAs", len(a3bs_circs)),
    ("n_abs_circRNAs", len(abs_circs)),
    ("abs_circRNA_fraction", len(abs_circs) / denom),
    ("median_a5bs_pcu", median(a5bs_pcus)),
    ("median_a3bs_pcu", median(a3bs_pcus)),
]

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


def unique_bins(column, preferred=()):
    seen = []
    for value in preferred:
        if value not in seen:
            seen.append(value)
    for row in rows:
        value = row[column]
        if value not in seen:
            seen.append(value)
    return seen


bsj_bins = [f"[{bsj_edges[i]},{bsj_edges[i+1]})" for i in range(len(bsj_edges) - 1)] + [f">={bsj_edges[-1]}"]
fsj_bins = [f"[{fsj_edges[i]},{fsj_edges[i+1]})" for i in range(len(fsj_edges) - 1)] + [f">={fsj_edges[-1]}"]
ratio_bins = [f"[{ratio_edges[i]},{ratio_edges[i+1]})" for i in range(len(ratio_edges) - 1)] + [f">={ratio_edges[-1]}"]
junction_ratio_bins = [
    f"[{junction_ratio_edges[i]},{junction_ratio_edges[i+1]})"
    for i in range(len(junction_ratio_edges) - 1)
] + [f">={junction_ratio_edges[-1]}"]
span_bins = [f"[{span_edges[i]},{span_edges[i+1]})" for i in range(len(span_edges) - 1)] + [f">={span_edges[-1]}"]
splice_class_bins = [
    "canonical_GU-AG",
    "semi_canonical_GC-AG",
    "minor_canonical_AU-AC",
    "non_canonical",
    "unknown",
]
abs_class_bins = ["no_abs", "A5BS_only", "A3BS_only", "A5BS_and_A3BS"]
circ_type_bins = unique_bins(
    "circRNA_type",
    preferred=("exon", "intron", "intergenic_region", "intergenic", "NA"),
)
pair_bins = sorted({r["splice_site_pair_rna"] for r in rows})

bin_specs = [
    ("bsj_count", "bsj_bin", bsj_bins),
    ("fsj_count", "fsj_bin", fsj_bins),
    ("bsj_fsj_ratio", "ratio_bin", ratio_bins),
    ("junction_ratio", "junction_ratio_bin", junction_ratio_bins),
    ("circ_span", "span_bin", span_bins),
    ("splice_site_class", "splice_site_class", splice_class_bins),
    ("abs_class", "abs_class", abs_class_bins),
    ("circRNA_type", "circRNA_type", circ_type_bins),
    ("splice_site_pair_rna", "splice_site_pair_rna", pair_bins),
]

with dist_out.open("w", newline="") as fh:
    writer = csv.writer(fh, delimiter="\t")
    writer.writerow(["sample", "feature", "bin", "count", "fraction"])
    for feature, col, ordered in bin_specs:
        counts = count_bins(col, ordered)
        for b in ordered:
            count = counts.get(b, 0)
            writer.writerow([sample, feature, b, count, count / denom])

def clean_axis(ax, grid_axis="y"):
    ax.grid(axis=grid_axis, color="#E6E6E6", linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def empty_panel(ax, message):
    ax.axis("off")
    ax.text(
        0.5,
        0.5,
        message,
        ha="center",
        va="center",
        color="#6B6B6B",
        fontsize=10,
        transform=ax.transAxes,
    )


def short_event_label(event_id):
    parts = event_id.split("|")
    if len(parts) >= 4:
        return f"{parts[0]} {parts[3]}"
    return event_id


fig, axes = plt.subplots(2, 3, figsize=(15.5, 8.4), squeeze=False)
axes = axes.flatten()
fig.suptitle(f"{sample} back-splicing summary", fontsize=14, fontweight="bold")

summary_lines = [
    ("circRNAs", len(rows)),
    ("BSJ > 0", sum(v > 0 for v in bsj_list)),
    ("total BSJ reads", int(sum(bsj_list))),
    ("median junction ratio", f"{median(junction_ratio_list):.3g}"),
    ("A5BS events", len(a5bs_events)),
    ("A3BS events", len(a3bs_events)),
    ("ABS circRNAs", len(abs_circs)),
]
axes[0].axis("off")
for i, (label, value) in enumerate(summary_lines):
    y = 0.92 - i * 0.13
    axes[0].text(0.05, y, label, ha="left", va="center", color="#555555", fontsize=10)
    axes[0].text(
        0.95,
        y,
        str(value),
        ha="right",
        va="center",
        color="#1F1F1F",
        fontsize=12,
        fontweight="bold",
    )
axes[0].set_title("Key metrics", loc="left")

bsj_labels = ["0", "1", "2-4", "5-9", "10-19", "20-49", "50+"]
bsj_counts = count_bins("bsj_bin", bsj_bins)
axes[1].bar(
    range(len(bsj_bins)),
    [bsj_counts.get(b, 0) for b in bsj_bins],
    color="#4E79A7",
)
axes[1].set_xticks(range(len(bsj_bins)))
axes[1].set_xticklabels(bsj_labels, fontsize=8)
axes[1].set_ylabel("circRNA count")
axes[1].set_title("BSJ read support")
clean_axis(axes[1])

junction_ratio_labels = ["0-.05", ".05-.1", ".1-.25", ".25-.5", ".5-.75", ".75-.9", ".9-.95", ".95-1", "1"]
junction_ratio_counts = count_bins("junction_ratio_bin", junction_ratio_bins)
axes[2].bar(
    range(len(junction_ratio_bins)),
    [junction_ratio_counts.get(b, 0) for b in junction_ratio_bins],
    color="#59A14F",
)
axes[2].set_xticks(range(len(junction_ratio_bins)))
axes[2].set_xticklabels(junction_ratio_labels, fontsize=8)
axes[2].set_ylim(bottom=0)
axes[2].set_ylabel("circRNA count")
axes[2].set_title("Junction ratio")
clean_axis(axes[2])

class_labels = ["GU-AG", "GC-AG", "AU-AC", "non-canonical", "unknown"]
class_counts = count_bins("splice_site_class", splice_class_bins)
class_values = [class_counts.get(b, 0) / denom for b in splice_class_bins]
axes[3].barh(class_labels, class_values, color=["#4E79A7", "#76B7B2", "#F28E2B", "#E15759", "#BAB0AC"])
axes[3].set_xlim(0, 1)
axes[3].set_xlabel("fraction")
axes[3].set_title("Splice-site class")
clean_axis(axes[3], grid_axis="x")

event_labels = ["A5BS events", "A3BS events", "ABS circRNAs"]
event_values = [len(a5bs_events), len(a3bs_events), len(abs_circs)]
axes[4].barh(event_labels, event_values, color=["#F28E2B", "#E15759", "#76B7B2"])
axes[4].set_xlabel("count")
axes[4].set_title("Alternative back-splicing")
clean_axis(axes[4], grid_axis="x")

event_totals = {}
for row in abs_rows:
    event_id = row["event_id"]
    event_totals[event_id] = (
        row["event_type"],
        safe_float(row["event_bsj_total"]),
    )
top_events = sorted(event_totals.items(), key=lambda item: (-item[1][1], item[0]))[:8]
if top_events:
    labels = [short_event_label(event_id) for event_id, _ in reversed(top_events)]
    values = [value for _, (_, value) in reversed(top_events)]
    colors = [
        "#F28E2B" if event_type == "A5BS" else "#E15759"
        for _, (event_type, _) in reversed(top_events)
    ]
    axes[5].barh(labels, values, color=colors)
    axes[5].set_xlabel("event BSJ reads")
    axes[5].set_title("Top ABS events")
    clean_axis(axes[5], grid_axis="x")
else:
    empty_panel(axes[5], "No alternative back-splicing events detected")
    axes[5].set_title("Top ABS events", loc="left")

fig.tight_layout(rect=[0, 0, 1, 0.94])
plot_out.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(plot_out, dpi=150)
plt.close(fig)
