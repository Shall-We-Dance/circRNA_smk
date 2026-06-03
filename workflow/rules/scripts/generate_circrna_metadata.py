"""Generate Sample / Condition / Rep metadata for circrna_plots from deg.groups.

Replicate numbers are auto-assigned by 1-based position within each group, so
the first sample listed under a group is Rep 1, the second is Rep 2, etc.
The order of conditions follows the order of keys in deg.groups (Python 3.7+
preserves dict insertion order, which the pipeline relies on elsewhere).
"""
import csv
from pathlib import Path

groups = snakemake.params.groups
out_path = Path(snakemake.output.metadata)
log_path = Path(snakemake.log[0]) if snakemake.log else None


def log_message(message):
    if log_path is None:
        return
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a") as fh:
        fh.write(message + "\n")


if not isinstance(groups, dict) or not groups:
    raise ValueError("snakemake.params.groups must be a non-empty mapping of condition -> samples.")

out_path.parent.mkdir(parents=True, exist_ok=True)
n_rows = 0
with out_path.open("w", newline="") as fh:
    writer = csv.writer(fh, delimiter="\t")
    writer.writerow(["Sample", "Condition", "Rep"])
    for condition, samples in groups.items():
        if isinstance(samples, str) or not hasattr(samples, "__iter__"):
            raise ValueError(
                f"deg.groups.{condition} must be a list of sample names; got {samples!r}."
            )
        for rep_idx, sample in enumerate(samples, start=1):
            writer.writerow([str(sample), str(condition), rep_idx])
            n_rows += 1

log_message(
    f"Wrote circrna_plots metadata for {len(groups)} condition(s) "
    f"and {n_rows} sample(s) to {out_path}."
)

