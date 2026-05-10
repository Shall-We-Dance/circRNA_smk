from pathlib import Path


def abs_path(path):
    return str(Path(path).resolve())


group_names = list(snakemake.params.group_names)
groups = {group: list(samples) for group, samples in dict(snakemake.params.groups).items()}
ordered_samples = [
    sample
    for group_name in group_names
    for sample in groups[group_name]
]
bams = [abs_path(path) for path in snakemake.input.bams]

if len(ordered_samples) != len(bams):
    raise ValueError(
        f"Expected {len(ordered_samples)} BAMs for all-group sashimi plotting, "
        f"got {len(bams)}."
    )
if len(group_names) < 2:
    raise ValueError("All-group sashimi plotting requires at least two groups.")

sample_to_bam = dict(zip(ordered_samples, bams))
sample_to_index = {sample: index for index, sample in enumerate(ordered_samples, start=1)}

first_group = group_names[0]
first_group_samples = groups[first_group]
other_samples = [
    sample
    for group_name in group_names[1:]
    for sample in groups[group_name]
]

for output_path in [
    snakemake.output.b1,
    snakemake.output.b2,
    snakemake.output.group_info,
    snakemake.output.samples,
]:
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)

Path(snakemake.output.b1).write_text(
    ",".join(sample_to_bam[sample] for sample in first_group_samples) + "\n"
)
Path(snakemake.output.b2).write_text(
    ",".join(sample_to_bam[sample] for sample in other_samples) + "\n"
)

with open(snakemake.output.group_info, "w") as handle:
    for group_name in group_names:
        indices = [str(sample_to_index[sample]) for sample in groups[group_name]]
        handle.write(f"{group_name}: {','.join(indices)}\n")

with open(snakemake.output.samples, "w") as handle:
    handle.write("sample\tgroup\trmats2sashimi_index\tbam\n")
    for group_name in group_names:
        for sample in groups[group_name]:
            handle.write(
                f"{sample}\t{group_name}\t{sample_to_index[sample]}\t{sample_to_bam[sample]}\n"
            )

Path(snakemake.log[0]).parent.mkdir(parents=True, exist_ok=True)
Path(snakemake.log[0]).write_text(
    "Prepared all-group rmats2sashimiplot grouping: "
    + ", ".join(f"{group} ({len(groups[group])} BAMs)" for group in group_names)
    + ".\n"
)
