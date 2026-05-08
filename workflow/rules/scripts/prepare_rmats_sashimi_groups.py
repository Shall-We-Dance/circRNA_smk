from pathlib import Path


def abs_path(path):
    return str(Path(path).resolve())


group_names = list(snakemake.params.group_names)
groups = {group: list(samples) for group, samples in dict(snakemake.params.groups).items()}
reference_group = str(snakemake.params.reference_group)
b1_samples = list(snakemake.params.b1_samples)
b2_samples = list(snakemake.params.b2_samples)
ordered_samples = b1_samples + b2_samples
bams = [abs_path(path) for path in snakemake.input.bams]

if len(ordered_samples) != len(bams):
    raise ValueError(
        f"Expected {len(ordered_samples)} BAMs for sashimi all-condition plotting, "
        f"got {len(bams)}."
    )

sample_to_bam = dict(zip(ordered_samples, bams))
sample_to_index = {sample: index for index, sample in enumerate(ordered_samples, start=1)}

for output_path in [
    snakemake.output.b1,
    snakemake.output.b2,
    snakemake.output.group_info,
    snakemake.output.samples,
]:
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)

Path(snakemake.output.b1).write_text(
    ",".join(sample_to_bam[sample] for sample in b1_samples) + "\n"
)
Path(snakemake.output.b2).write_text(
    ",".join(sample_to_bam[sample] for sample in b2_samples) + "\n"
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
    f"Prepared rmats2sashimiplot grouping with {len(group_names)} groups. "
    f"Reference b1 group: {reference_group}.\n"
)
