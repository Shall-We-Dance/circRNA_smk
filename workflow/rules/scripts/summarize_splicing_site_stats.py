import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


summary_files = [Path(p) for p in snakemake.input.summaries]
dist_files = [Path(p) for p in snakemake.input.distributions]
samples = list(snakemake.params.samples)

summary_out = Path(snakemake.output.summary)
dists_out = Path(snakemake.output.distributions)
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


def metric_map(metric_name):
    m = {}
    for row in summary_rows:
        if row["metric"] == metric_name:
            try:
                m[row["sample"]] = float(row["value"])
            except ValueError:
                m[row["sample"]] = 0.0
    return m


n_circ = metric_map("n_circRNAs")
med_ratio = metric_map("median_ratio")
canonical_fraction = metric_map("canonical_GU_AG_fraction")

ratio_rows = [r for r in dist_rows if r["feature"] == "bsj_fsj_ratio"]
ratio_bins = sorted({r["bin"] for r in ratio_rows})

ratio_fraction = {s: {b: 0.0 for b in ratio_bins} for s in samples}
for row in ratio_rows:
    sample = row["sample"]
    b = row["bin"]
    if sample not in ratio_fraction:
        ratio_fraction[sample] = {x: 0.0 for x in ratio_bins}
    try:
        ratio_fraction[sample][b] = float(row["fraction"])
    except ValueError:
        ratio_fraction[sample][b] = 0.0

class_rows = [r for r in dist_rows if r["feature"] == "splice_site_class"]
class_bins = [
    "canonical_GU-AG",
    "semi_canonical_GC-AG",
    "minor_canonical_AU-AC",
    "non_canonical",
    "unknown",
]
class_fraction = {s: {b: 0.0 for b in class_bins} for s in samples}
for row in class_rows:
    sample = row["sample"]
    b = row["bin"]
    if b not in class_bins:
        continue
    try:
        class_fraction[sample][b] = float(row["fraction"])
    except ValueError:
        class_fraction[sample][b] = 0.0

fig, axes = plt.subplots(1, 4, figsize=(20, 4.5))

xs = list(range(len(samples)))
axes[0].bar(xs, [n_circ.get(s, 0.0) for s in samples], color="#55A868")
axes[0].set_xticks(xs)
axes[0].set_xticklabels(samples, rotation=45, ha="right")
axes[0].set_ylabel("count")
axes[0].set_title("Detected circRNAs")

axes[1].bar(xs, [med_ratio.get(s, 0.0) for s in samples], color="#C44E52")
axes[1].set_xticks(xs)
axes[1].set_xticklabels(samples, rotation=45, ha="right")
axes[1].set_ylabel("median (BSJ+1)/(FSJ+1)")
axes[1].set_title("Median BSJ/FSJ ratio")

axes[2].bar(xs, [canonical_fraction.get(s, 0.0) for s in samples], color="#4C72B0")
axes[2].set_xticks(xs)
axes[2].set_xticklabels(samples, rotation=45, ha="right")
axes[2].set_ylabel("fraction")
axes[2].set_ylim(0, 1)
axes[2].set_title("Canonical GU-AG fraction")

bottom = [0.0 for _ in samples]
colors = plt.cm.Blues([0.35 + 0.55 * i / max(len(class_bins), 1) for i in range(len(class_bins))])
for color, b in zip(colors, class_bins):
    values = [class_fraction.get(s, {}).get(b, 0.0) for s in samples]
    axes[3].bar(xs, values, bottom=bottom, label=b, color=color)
    bottom = [x + y for x, y in zip(bottom, values)]
axes[3].set_xticks(xs)
axes[3].set_xticklabels(samples, rotation=45, ha="right")
axes[3].set_ylim(0, 1)
axes[3].set_ylabel("fraction")
axes[3].set_title("Splice-site class composition")
axes[3].legend(fontsize=7, title="site class", loc="upper left", bbox_to_anchor=(1.02, 1.0))

fig.tight_layout()
plot_out.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(plot_out, dpi=150, bbox_inches="tight")
plt.close(fig)
