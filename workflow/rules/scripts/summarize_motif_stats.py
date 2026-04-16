import csv
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
plot_out = Path(snakemake.output.plot)


sample_rows = {}
for sample, path in zip(samples, site_files):
    rows = []
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            rows.append(row)
    sample_rows[sample] = rows


def get_pair(seq):
    donor = (seq.get("donor_window_seq") or "NN")[:2].upper().replace("T", "U")
    acceptor = (seq.get("acceptor_window_seq") or "NN")[-2:].upper().replace("T", "U")
    return f"{donor}-{acceptor}"


metrics = []
pair_fraction = defaultdict(dict)
all_pairs = set()
all_donor = defaultdict(float)
all_acceptor = defaultdict(float)

for sample in samples:
    rows = sample_rows.get(sample, [])
    n_sites = len(rows)
    total_weight = 0.0
    guag_w = 0.0
    pair_weight = defaultdict(float)
    for row in rows:
        try:
            w = float(row.get("weight", 1.0))
        except ValueError:
            w = 1.0
        total_weight += w
        pair = get_pair(row)
        pair_weight[pair] += w
        donor, acceptor = pair.split("-")
        all_donor[donor] += w
        all_acceptor[acceptor] += w
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
            rec = {"sample": sample}
            for k, v in zip(header, row):
                rec[k] = v
            known_rows.append(rec)

known_fields = ["sample"]
if known_rows:
    extras = list(known_rows[0].keys())
    known_fields = ["sample"] + [x for x in extras if x != "sample"]

with known_out.open("w", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=known_fields, delimiter="\t")
    writer.writeheader()
    for row in known_rows:
        writer.writerow(row)


metric_map = defaultdict(dict)
for row in metrics:
    metric_map[row["metric"]][row["sample"]] = float(row["value"])

pair_bins = sorted(all_pairs)
if len(pair_bins) > 8:
    ranked = sorted(
        pair_bins,
        key=lambda p: sum(pair_fraction[s].get(p, 0.0) for s in samples),
        reverse=True,
    )
    pair_bins = ranked[:8]

fig, axes = plt.subplots(1, 3, figsize=(18, 4.5))
xs = list(range(len(samples)))

axes[0].bar(xs, [metric_map["n_bsj_sites"].get(s, 0.0) for s in samples], color="#55A868")
axes[0].set_xticks(xs)
axes[0].set_xticklabels(samples, rotation=45, ha="right")
axes[0].set_ylabel("count")
axes[0].set_title("BSJ sites per sample")

axes[1].bar(
    xs,
    [metric_map["canonical_GU_AG_weighted_fraction"].get(s, 0.0) for s in samples],
    color="#4C72B0",
)
axes[1].set_xticks(xs)
axes[1].set_xticklabels(samples, rotation=45, ha="right")
axes[1].set_ylabel("fraction")
axes[1].set_ylim(0, 1)
axes[1].set_title("Canonical GU-AG weighted fraction")

bottom = [0.0 for _ in samples]
colors = plt.cm.PuBu([0.35 + 0.55 * i / max(len(pair_bins), 1) for i in range(len(pair_bins))])
for color, pair in zip(colors, pair_bins):
    vals = [pair_fraction[s].get(pair, 0.0) for s in samples]
    axes[2].bar(xs, vals, bottom=bottom, color=color, label=pair)
    bottom = [x + y for x, y in zip(bottom, vals)]
axes[2].set_xticks(xs)
axes[2].set_xticklabels(samples, rotation=45, ha="right")
axes[2].set_ylabel("fraction")
axes[2].set_ylim(0, 1)
axes[2].set_title("Splice dinucleotide pair composition")
axes[2].legend(fontsize=7, title="pair", loc="upper left", bbox_to_anchor=(1.02, 1.0))

fig.tight_layout()
plot_out.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(plot_out, dpi=150, bbox_inches="tight")
plt.close(fig)
