from pathlib import Path
import csv

pair = snakemake.params.pair
group_a = snakemake.params.group_a
group_b = snakemake.params.group_b
groups = snakemake.params.groups
sample_names = list(snakemake.params.sample_names)

selected_samples = list(groups[group_a]) + list(groups[group_b])
sample_to_class = {sample: group_a for sample in groups[group_a]}
sample_to_class.update({sample: group_b for sample in groups[group_b]})

input_ciri3_files = [Path(p) for p in snakemake.input.ciri3]
sample_to_ciri3 = dict(zip(sample_names, input_ciri3_files))

for sample in selected_samples:
    if sample not in sample_to_ciri3:
        raise ValueError(f"Sample '{sample}' required by pair '{pair}' not found in per-sample CIRI3 outputs")

infor_out = Path(snakemake.output.infor)
infor_out.parent.mkdir(parents=True, exist_ok=True)

with infor_out.open("w", newline="") as fh:
    writer = csv.writer(fh, delimiter="\t")
    writer.writerow(["Sample", "Path", "Class"])
    for sample in selected_samples:
        writer.writerow([sample, str(sample_to_ciri3[sample].resolve()), sample_to_class[sample]])
