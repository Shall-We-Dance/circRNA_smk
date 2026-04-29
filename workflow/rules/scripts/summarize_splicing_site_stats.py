import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


summary_files = [Path(p) for p in snakemake.input.summaries]
dist_files = [Path(p) for p in snakemake.input.distributions]
abs_files = [Path(p) for p in snakemake.input.abs_events]
samples = list(snakemake.params.samples)

summary_out = Path(snakemake.output.summary)
dists_out = Path(snakemake.output.distributions)
abs_out = Path(snakemake.output.abs_events)
plot_out = Path(snakemake.output.plot)


summary_rows = []
for path in summary_files:
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            summary_rows.append(row)

summary_out.parent.mkdir(parents=True, exist_ok=True)
with summary_out.open("w", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=["sample", "metric", "value"], delimiter="\t")
    writer.writeheader()
    for row in summary_rows:
        writer.writerow(row)

dist_rows = []
for path in dist_files:
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            dist_rows.append(row)

with dists_out.open("w", newline="") as fh:
    writer = csv.DictWriter(
        fh,
        fieldnames=["sample", "feature", "bin", "count", "fraction"],
        delimiter="\t",
    )
    writer.writeheader()
    for row in dist_rows:
        writer.writerow(row)

abs_rows = []
abs_fieldnames = []
for path in abs_files:
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames:
            for fieldname in reader.fieldnames:
                if fieldname not in abs_fieldnames:
                    abs_fieldnames.append(fieldname)
        for row in reader:
            abs_rows.append(row)

if not abs_fieldnames:
    abs_fieldnames = [
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

with abs_out.open("w", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=abs_fieldnames, delimiter="\t")
    writer.writeheader()
    for row in abs_rows:
        writer.writerow(row)


def metric_map(metric_name):
    m = {}
    for row in summary_rows:
        if row["metric"] == metric_name:
            try:
                m[row["sample"]] = float(row["value"])
            except ValueError:
                m[row["sample"]] = 0.0
    return m


def dist_fraction_map(feature_name, bins):
    values = {s: {b: 0.0 for b in bins} for s in samples}
    for row in dist_rows:
        if row["feature"] != feature_name:
            continue
        sample = row["sample"]
        b = row["bin"]
        if sample not in values:
            values[sample] = {x: 0.0 for x in bins}
        if b not in values[sample]:
            values[sample][b] = 0.0
        try:
            values[sample][b] = float(row["fraction"])
        except ValueError:
            values[sample][b] = 0.0
    return values


bsj_detected = metric_map("n_bsj_gt_0")
total_bsj = metric_map("total_bsj_reads")
med_junction_ratio = metric_map("median_junction_ratio")
n_a5bs = metric_map("n_a5bs_events")
n_a3bs = metric_map("n_a3bs_events")

class_bins = [
    "canonical_GU-AG",
    "semi_canonical_GC-AG",
    "minor_canonical_AU-AC",
    "non_canonical",
    "unknown",
]
class_fraction = dist_fraction_map("splice_site_class", class_bins)

abs_class_bins = ["no_abs", "A5BS_only", "A3BS_only", "A5BS_and_A3BS"]
abs_class_fraction = dist_fraction_map("abs_class", abs_class_bins)


def clean_axis(ax, grid_axis="x"):
    ax.grid(axis=grid_axis, color="#E6E6E6", linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


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


def stacked_horizontal(ax, fraction_map, bins, labels, colors, title):
    y = list(range(len(samples)))
    left = [0.0 for _ in samples]
    for b, label, color in zip(bins, labels, colors):
        values = [fraction_map.get(s, {}).get(b, 0.0) for s in samples]
        ax.barh(y, values, left=left, label=label, color=color)
        left = [x + v for x, v in zip(left, values)]
    ax.set_yticks(y)
    ax.set_yticklabels(samples, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlim(0, 1)
    ax.set_xlabel("fraction")
    ax.set_title(title)
    ax.legend(fontsize=7, loc="lower center", bbox_to_anchor=(0.5, -0.34), ncol=2)
    clean_axis(ax)


fig_height = max(8.5, 3.5 + 0.42 * len(samples))
fig, axes = plt.subplots(2, 3, figsize=(18, fig_height))
axes = axes.flatten()
fig.suptitle("Back-splicing overview", fontsize=14, fontweight="bold")

horizontal_metric(
    axes[0],
    [bsj_detected.get(s, 0.0) for s in samples],
    "BSJ-supported circRNAs",
    "count",
    "#4E79A7",
)

horizontal_metric(
    axes[1],
    [total_bsj.get(s, 0.0) for s in samples],
    "Total BSJ reads",
    "reads",
    "#76B7B2",
)

horizontal_metric(
    axes[2],
    [med_junction_ratio.get(s, 0.0) for s in samples],
    "Median junction ratio",
    "2*BSJ/(2*BSJ+FSJ)",
    "#59A14F",
    xmax=1,
)

y = list(range(len(samples)))
a5_values = [n_a5bs.get(s, 0.0) for s in samples]
a3_values = [n_a3bs.get(s, 0.0) for s in samples]
axes[3].barh(y, a5_values, color="#F28E2B", label="A5BS")
axes[3].barh(y, a3_values, left=a5_values, color="#E15759", label="A3BS")
axes[3].set_yticks(y)
axes[3].set_yticklabels(samples, fontsize=8)
axes[3].invert_yaxis()
axes[3].set_xlabel("event count")
axes[3].set_title("Alternative back-splicing events")
axes[3].legend(fontsize=8)
clean_axis(axes[3])

stacked_horizontal(
    axes[4],
    class_fraction,
    class_bins,
    ["GU-AG", "GC-AG", "AU-AC", "non-canonical", "unknown"],
    ["#4E79A7", "#76B7B2", "#F28E2B", "#E15759", "#BAB0AC"],
    "Splice-site class composition",
)

stacked_horizontal(
    axes[5],
    abs_class_fraction,
    abs_class_bins,
    ["no ABS", "A5BS only", "A3BS only", "A5BS+A3BS"],
    ["#D0D0D0", "#F28E2B", "#E15759", "#9467BD"],
    "ABS circRNA composition",
)

fig.tight_layout(rect=[0, 0.03, 1, 0.95])
plot_out.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(plot_out, dpi=150, bbox_inches="tight")
plt.close(fig)
