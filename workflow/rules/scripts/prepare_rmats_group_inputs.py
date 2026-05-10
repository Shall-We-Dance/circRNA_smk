from pathlib import Path


def abs_path(path):
    return str(Path(path).resolve())


group = str(snakemake.params.group)
samples = list(snakemake.params.samples)
bams = [abs_path(path) for path in snakemake.input.bams]

if len(samples) != len(bams):
    raise ValueError(
        f"Expected {len(samples)} BAM files for rMATS group {group}, got {len(bams)}."
    )

for output_path in [snakemake.output.b1, snakemake.output.samples]:
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)

Path(snakemake.output.b1).write_text(",".join(bams) + "\n")

with open(snakemake.output.samples, "w") as handle:
    handle.write("sample\tgroup\tbam\n")
    for sample, bam in zip(samples, bams):
        handle.write(f"{sample}\t{group}\t{bam}\n")

Path(snakemake.log[0]).parent.mkdir(parents=True, exist_ok=True)
Path(snakemake.log[0]).write_text(
    f"Prepared rMATS single-group inputs for {group}: {len(bams)} BAMs.\n"
)
