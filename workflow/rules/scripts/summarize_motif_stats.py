import csv
import math
import re
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


site_files = [Path(p) for p in snakemake.input.site_tables]
known_files = [Path(p) for p in snakemake.input.known_results]
samples = list(snakemake.params.samples)

site_stats_out = Path(snakemake.output.site_stats)
known_out = Path(snakemake.output.known_merged)
homer_summary_out = Path(snakemake.output.homer_summary)
top_motifs_out = Path(snakemake.output.top_motifs)
plot_out = Path(snakemake.output.plot)
heatmap_out = Path(snakemake.output.heatmap)


def safe_float(value):
    if value is None:
        return None
    text = str(value).strip().replace(",", "").replace("%", "")
    if not text or text.upper() in {"NA", "N/A", "NAN", "INF", "-INF"}:
        return None
    match = re.search(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?", text)
    if not match:
        return None
    try:
        parsed = float(match.group(0))
    except ValueError:
        return None
    if math.isnan(parsed) or math.isinf(parsed):
        return None
    return parsed


def format_number(value, default="NA"):
    if value is None:
        return default
    return f"{value:.6g}"


def median(values):
    values = sorted(v for v in values if v is not None)
    if not values:
        return None
    mid = len(values) // 2
    if len(values) % 2:
        return values[mid]
    return (values[mid - 1] + values[mid]) / 2


def get_pair(seq):
    donor = (seq.get("donor_window_seq") or "NN")[:2].upper().replace("T", "U")
    acceptor = (seq.get("acceptor_window_seq") or "NN")[-2:].upper().replace("T", "U")
    return f"{donor}-{acceptor}"


def normalize_homer_row(row):
    """
    HOMER sometimes writes dynamic column names like:
      # of Target Sequences with Motif(of 2758)
      # of Background Sequences with Motif(of 12390)

    Normalize those into stable column names so rows from different samples
    can be merged safely.
    """
    new_row = {}

    for key, value in row.items():
        if key is None:
            continue

        key = key.strip()

        m_target = re.match(r"^# of Target Sequences with Motif\(of (\d+)\)$", key)
        if m_target:
            new_row["Target_sequences_with_motif"] = value
            new_row["Target_sequences_total"] = m_target.group(1)
            continue

        m_bg = re.match(r"^# of Background Sequences with Motif\(of (\d+)\)$", key)
        if m_bg:
            new_row["Background_sequences_with_motif"] = value
            new_row["Background_sequences_total"] = m_bg.group(1)
            continue

        if key == "# of Target Sequences with Motif":
            new_row["Target_sequences_with_motif"] = value
            continue

        if key == "# of Background Sequences with Motif":
            new_row["Background_sequences_with_motif"] = value
            continue

        new_row[key] = value

    return new_row


def first_value(row, names):
    for name in names:
        value = row.get(name)
        if value not in (None, ""):
            return value
    return None


def motif_name(row):
    return (
        first_value(row, ("Motif Name", "Name", "Motif", "Known Motif"))
        or "unknown_motif"
    )


def motif_consensus(row):
    return first_value(row, ("Consensus", "Motif Consensus", "consensus")) or "NA"


def motif_key(row):
    return motif_name(row), motif_consensus(row)


def short_label(text, limit=38):
    text = str(text)
    if len(text) <= limit:
        return text
    return text[: limit - 1] + "..."


def enrich_homer_row(row):
    p_value = safe_float(first_value(row, ("P-value", "P Value", "Pvalue")))
    q_value = safe_float(
        first_value(
            row,
            (
                "q-value (Benjamini)",
                "q-value",
                "FDR",
                "Benjamini",
                "qvalue",
            ),
        )
    )
    log_p = safe_float(first_value(row, ("Log P-value", "log P-value", "LogP", "logP")))

    if p_value is not None and p_value > 0:
        neg_log10_p = -math.log10(p_value)
    elif log_p is not None:
        neg_log10_p = abs(log_p)
    else:
        neg_log10_p = 0.0

    target_pct = safe_float(
        first_value(row, ("% of Target Sequences", "% Target", "Target %"))
    )
    background_pct = safe_float(
        first_value(
            row,
            ("% of Background Sequences", "% Background", "Background %"),
        )
    )
    if target_pct is None:
        target_count = safe_float(row.get("Target_sequences_with_motif"))
        target_total = safe_float(row.get("Target_sequences_total"))
        if target_count is not None and target_total is not None and target_total > 0:
            target_pct = 100.0 * target_count / target_total
    if background_pct is None:
        background_count = safe_float(row.get("Background_sequences_with_motif"))
        background_total = safe_float(row.get("Background_sequences_total"))
        if (
            background_count is not None
            and background_total is not None
            and background_total > 0
        ):
            background_pct = 100.0 * background_count / background_total
    enrichment = None
    if target_pct is not None and background_pct is not None and background_pct > 0:
        enrichment = target_pct / background_pct

    row["P_value_numeric"] = format_number(p_value)
    row["q_value_numeric"] = format_number(q_value)
    row["neg_log10_p"] = format_number(neg_log10_p)
    row["target_pct"] = format_number(target_pct)
    row["background_pct"] = format_number(background_pct)
    row["target_vs_background_enrichment"] = format_number(enrichment)
    return row


sample_rows = {}
for sample, path in zip(samples, site_files):
    rows = []
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            rows.append(row)
    sample_rows[sample] = rows


metrics = []
pair_fraction = defaultdict(dict)
all_pairs = set()

for sample in samples:
    rows = sample_rows.get(sample, [])
    n_sites = len(rows)
    total_weight = 0.0
    guag_w = 0.0
    pair_weight = defaultdict(float)

    for row in rows:
        try:
            w = float(row.get("weight", 1.0))
        except (ValueError, TypeError):
            w = 1.0

        total_weight += w
        pair = get_pair(row)
        pair_weight[pair] += w

        if pair == "GU-AG":
            guag_w += w

    denom = total_weight if total_weight > 0 else 1.0

    metrics.extend(
        [
            {"sample": sample, "metric": "n_bsj_sites", "value": n_sites},
            {"sample": sample, "metric": "sum_weight", "value": f"{total_weight:.6g}"},
            {
                "sample": sample,
                "metric": "canonical_GU_AG_weighted_fraction",
                "value": f"{(guag_w / denom):.6g}",
            },
        ]
    )

    for pair, w in pair_weight.items():
        pair_fraction[sample][pair] = w / denom
        all_pairs.add(pair)


site_stats_out.parent.mkdir(parents=True, exist_ok=True)
with site_stats_out.open("w", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=["sample", "metric", "value"], delimiter="\t")
    writer.writeheader()
    writer.writerows(metrics)


known_rows = []
known_by_sample = defaultdict(list)
all_known_fields = {"sample"}

for sample, path in zip(samples, known_files):
    with path.open() as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = None

        for row in reader:
            if not row:
                continue
            if row[0].startswith("#"):
                continue

            if header is None:
                header = row
                continue

            rec = {k: v for k, v in zip(header, row)}
            rec = normalize_homer_row(rec)
            rec["sample"] = sample
            rec = enrich_homer_row(rec)

            known_rows.append(rec)
            known_by_sample[sample].append(rec)
            all_known_fields.update(rec.keys())


preferred_known_fields = [
    "sample",
    "Motif Name",
    "Consensus",
    "P-value",
    "P_value_numeric",
    "neg_log10_p",
    "Log P-value",
    "log P-value",
    "q-value (Benjamini)",
    "q_value_numeric",
    "# Target Sequences",
    "% of Target Sequences",
    "target_pct",
    "Target_sequences_with_motif",
    "Target_sequences_total",
    "# Background Sequences",
    "% of Background Sequences",
    "background_pct",
    "Background_sequences_with_motif",
    "Background_sequences_total",
    "target_vs_background_enrichment",
]

known_fields = [f for f in preferred_known_fields if f in all_known_fields]
known_fields += sorted(all_known_fields - set(known_fields))

known_out.parent.mkdir(parents=True, exist_ok=True)
with known_out.open("w", newline="") as fh:
    writer = csv.DictWriter(
        fh,
        fieldnames=known_fields,
        delimiter="\t",
        extrasaction="ignore",
    )
    writer.writeheader()
    for row in known_rows:
        writer.writerow(row)


metric_map = defaultdict(dict)
for row in metrics:
    metric_map[row["metric"]][row["sample"]] = float(row["value"])


def row_p(row):
    return safe_float(row.get("P_value_numeric"))


def row_q(row):
    return safe_float(row.get("q_value_numeric"))


def row_neg_log10_p(row):
    return safe_float(row.get("neg_log10_p")) or 0.0


def row_enrichment(row):
    return safe_float(row.get("target_vs_background_enrichment"))


homer_summary_rows = []
for sample in samples:
    rows = known_by_sample.get(sample, [])
    sorted_rows = sorted(rows, key=lambda r: (-row_neg_log10_p(r), motif_name(r)))
    top = sorted_rows[0] if sorted_rows else {}
    top_p = row_p(top)
    top_q = row_q(top)
    top_enrichment = row_enrichment(top)

    homer_summary_rows.append(
        {
            "sample": sample,
            "n_bsj_sites": format_number(metric_map["n_bsj_sites"].get(sample, 0.0)),
            "sum_weight": format_number(metric_map["sum_weight"].get(sample, 0.0)),
            "canonical_GU_AG_weighted_fraction": format_number(
                metric_map["canonical_GU_AG_weighted_fraction"].get(sample, 0.0)
            ),
            "n_known_motifs": len(rows),
            "n_known_motifs_q_lt_0_05": sum(
                1 for r in rows if (row_q(r) is not None and row_q(r) < 0.05)
            ),
            "n_known_motifs_p_lt_0_001": sum(
                1 for r in rows if (row_p(r) is not None and row_p(r) < 0.001)
            ),
            "top_motif": motif_name(top) if top else "NA",
            "top_consensus": motif_consensus(top) if top else "NA",
            "top_motif_p_value": format_number(top_p),
            "top_motif_q_value": format_number(top_q),
            "top_motif_neg_log10_p": format_number(row_neg_log10_p(top) if top else None),
            "top_motif_target_pct": top.get("target_pct", "NA") if top else "NA",
            "top_motif_background_pct": top.get("background_pct", "NA") if top else "NA",
            "top_motif_enrichment": format_number(top_enrichment),
        }
    )

homer_summary_fields = [
    "sample",
    "n_bsj_sites",
    "sum_weight",
    "canonical_GU_AG_weighted_fraction",
    "n_known_motifs",
    "n_known_motifs_q_lt_0_05",
    "n_known_motifs_p_lt_0_001",
    "top_motif",
    "top_consensus",
    "top_motif_p_value",
    "top_motif_q_value",
    "top_motif_neg_log10_p",
    "top_motif_target_pct",
    "top_motif_background_pct",
    "top_motif_enrichment",
]

homer_summary_out.parent.mkdir(parents=True, exist_ok=True)
with homer_summary_out.open("w", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=homer_summary_fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(homer_summary_rows)


motif_groups = defaultdict(list)
for row in known_rows:
    motif_groups[motif_key(row)].append(row)

top_motif_rows = []
for (name, consensus), rows in motif_groups.items():
    best = max(rows, key=lambda r: (row_neg_log10_p(r), motif_name(r)))
    q_samples = {
        r["sample"]
        for r in rows
        if row_q(r) is not None and row_q(r) < 0.05
    }
    p_samples = {
        r["sample"]
        for r in rows
        if row_p(r) is not None and row_p(r) < 0.001
    }
    top_motif_rows.append(
        {
            "motif_name": name,
            "consensus": consensus,
            "samples_detected": len({r["sample"] for r in rows}),
            "samples_q_lt_0_05": len(q_samples),
            "samples_p_lt_0_001": len(p_samples),
            "best_sample": best["sample"],
            "best_p_value": format_number(row_p(best)),
            "best_q_value": format_number(row_q(best)),
            "max_neg_log10_p": format_number(row_neg_log10_p(best)),
            "median_neg_log10_p": format_number(
                median([row_neg_log10_p(r) for r in rows])
            ),
            "median_target_pct": format_number(
                median([safe_float(r.get("target_pct")) for r in rows])
            ),
            "median_background_pct": format_number(
                median([safe_float(r.get("background_pct")) for r in rows])
            ),
            "median_enrichment": format_number(
                median([row_enrichment(r) for r in rows])
            ),
        }
    )

top_motif_rows.sort(
    key=lambda r: (
        -int(r["samples_q_lt_0_05"]),
        -int(r["samples_p_lt_0_001"]),
        -safe_float(r["max_neg_log10_p"] or 0),
        r["motif_name"],
    )
)

top_motif_fields = [
    "motif_name",
    "consensus",
    "samples_detected",
    "samples_q_lt_0_05",
    "samples_p_lt_0_001",
    "best_sample",
    "best_p_value",
    "best_q_value",
    "max_neg_log10_p",
    "median_neg_log10_p",
    "median_target_pct",
    "median_background_pct",
    "median_enrichment",
]

top_motifs_out.parent.mkdir(parents=True, exist_ok=True)
with top_motifs_out.open("w", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=top_motif_fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(top_motif_rows)


def clean_axis(ax, grid_axis="x"):
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


def horizontal_metric(ax, values, title, xlabel, color, xmax=None):
    y = list(range(len(samples)))
    ax.barh(y, values, color=color)
    ax.set_yticks(y)
    ax.set_yticklabels(samples, fontsize=8)
    ax.invert_yaxis()
    if xmax is not None:
        ax.set_xlim(0, xmax)
    ax.set_xlabel(xlabel)
    ax.set_title(title)
    clean_axis(ax)


def stacked_horizontal(ax, fraction_map, bins, colors, title):
    y = list(range(len(samples)))
    left = [0.0 for _ in samples]
    for pair, color in zip(bins, colors):
        vals = [fraction_map[s].get(pair, 0.0) for s in samples]
        ax.barh(y, vals, left=left, color=color, label=pair)
        left = [x + v for x, v in zip(left, vals)]

    ax.set_yticks(y)
    ax.set_yticklabels(samples, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlim(0, 1)
    ax.set_xlabel("weighted fraction")
    ax.set_title(title)
    ax.legend(fontsize=7, title="pair", loc="lower center", bbox_to_anchor=(0.5, -0.34), ncol=3)
    clean_axis(ax)


pair_bins = sorted(all_pairs)
if len(pair_bins) > 8:
    ranked = sorted(
        pair_bins,
        key=lambda p: sum(pair_fraction[s].get(p, 0.0) for s in samples),
        reverse=True,
    )
    pair_bins = ranked[:8]

pair_colors = [
    "#4E79A7",
    "#76B7B2",
    "#F28E2B",
    "#E15759",
    "#59A14F",
    "#B07AA1",
    "#EDC948",
    "#BAB0AC",
][: len(pair_bins)]

fig_height = max(8.5, 3.5 + 0.42 * len(samples))
fig, axes = plt.subplots(2, 3, figsize=(18, fig_height))
axes = axes.flatten()
fig.suptitle("HOMER BSJ motif summary", fontsize=14, fontweight="bold")

horizontal_metric(
    axes[0],
    [metric_map["n_bsj_sites"].get(s, 0.0) for s in samples],
    "BSJ sites exported to HOMER",
    "site count",
    "#4E79A7",
)

horizontal_metric(
    axes[1],
    [metric_map["canonical_GU_AG_weighted_fraction"].get(s, 0.0) for s in samples],
    "Canonical GU-AG weighted fraction",
    "fraction",
    "#59A14F",
    xmax=1,
)

y = list(range(len(samples)))
bar_height = 0.34
q_counts = [
    int(r["n_known_motifs_q_lt_0_05"])
    for r in homer_summary_rows
]
p_counts = [
    int(r["n_known_motifs_p_lt_0_001"])
    for r in homer_summary_rows
]
axes[2].barh([i - bar_height / 2 for i in y], q_counts, height=bar_height, color="#F28E2B", label="q < 0.05")
axes[2].barh([i + bar_height / 2 for i in y], p_counts, height=bar_height, color="#E15759", label="p < 0.001")
axes[2].set_yticks(y)
axes[2].set_yticklabels(samples, fontsize=8)
axes[2].invert_yaxis()
axes[2].set_xlabel("known motif count")
axes[2].set_title("Significant HOMER known motifs")
axes[2].legend(fontsize=8)
clean_axis(axes[2])

top_enrichments = [
    safe_float(r["top_motif_enrichment"]) or 0.0
    for r in homer_summary_rows
]
horizontal_metric(
    axes[3],
    top_enrichments,
    "Top motif target/background enrichment",
    "fold enrichment",
    "#76B7B2",
)
for y_pos, row in enumerate(homer_summary_rows):
    label = short_label(row["top_motif"], limit=24)
    if label != "NA":
        axes[3].text(
            top_enrichments[y_pos],
            y_pos,
            f"  {label}",
            va="center",
            ha="left",
            fontsize=7,
            color="#4A4A4A",
        )

if pair_bins:
    stacked_horizontal(
        axes[4],
        pair_fraction,
        pair_bins,
        pair_colors,
        "BSJ splice dinucleotide composition",
    )
else:
    empty_panel(axes[4], "No BSJ site sequences available")

if top_motif_rows:
    top_rows = top_motif_rows[:10]
    labels = [short_label(r["motif_name"], 34) for r in reversed(top_rows)]
    values = [safe_float(r["max_neg_log10_p"]) or 0.0 for r in reversed(top_rows)]
    axes[5].barh(labels, values, color="#B07AA1")
    axes[5].set_xlabel("-log10(best p-value)")
    axes[5].set_title("Top known motifs across samples")
    clean_axis(axes[5])
else:
    empty_panel(axes[5], "No HOMER known motif rows detected")
    axes[5].set_title("Top known motifs across samples", loc="left")

fig.tight_layout(rect=[0, 0.03, 1, 0.95])
plot_out.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(plot_out, dpi=150, bbox_inches="tight")
plt.close(fig)


def plot_heatmap():
    if not top_motif_rows:
        fig, ax = plt.subplots(1, 1, figsize=(9, 3.5))
        empty_panel(ax, "No HOMER known motif rows detected")
        ax.set_title("Top known motif significance")
        heatmap_out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(heatmap_out, dpi=150, bbox_inches="tight")
        plt.close(fig)
        return

    top_rows = top_motif_rows[: min(20, len(top_motif_rows))]
    keys = [(row["motif_name"], row["consensus"]) for row in top_rows]
    row_index = {key: i for i, key in enumerate(keys)}
    sample_index = {sample: i for i, sample in enumerate(samples)}
    matrix = [[0.0 for _ in samples] for _ in keys]

    for row in known_rows:
        key = motif_key(row)
        if key not in row_index:
            continue
        value = min(row_neg_log10_p(row), 50.0)
        matrix[row_index[key]][sample_index[row["sample"]]] = value

    width = max(8, 0.55 * len(samples) + 5)
    height = max(5, 0.32 * len(keys) + 2.5)
    fig, ax = plt.subplots(1, 1, figsize=(width, height))
    image = ax.imshow(matrix, aspect="auto", interpolation="nearest", cmap="magma")
    ax.set_xticks(range(len(samples)))
    ax.set_xticklabels(samples, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(keys)))
    ax.set_yticklabels([short_label(key[0], 42) for key in keys], fontsize=8)
    ax.set_title("Top HOMER known motifs by sample")
    cbar = fig.colorbar(image, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("-log10(P-value), capped at 50")
    fig.tight_layout()
    heatmap_out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(heatmap_out, dpi=150, bbox_inches="tight")
    plt.close(fig)


plot_heatmap()
