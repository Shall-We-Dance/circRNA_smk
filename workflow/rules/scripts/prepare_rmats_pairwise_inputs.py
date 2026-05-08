from pathlib import Path


def abs_path(path):
    return str(Path(path).resolve())


case_samples = list(snakemake.params.case_samples)
control_samples = list(snakemake.params.control_samples)
case_group = str(snakemake.params.case_group)
control_group = str(snakemake.params.control_group)
bams = [abs_path(path) for path in snakemake.input.bams]

expected_count = len(case_samples) + len(control_samples)
if len(bams) != expected_count:
    raise ValueError(
        f"Expected {expected_count} BAM files for {snakemake.wildcards.comparison}, "
        f"got {len(bams)}."
    )

case_bams = bams[: len(case_samples)]
control_bams = bams[len(case_samples) :]

for output_path in [snakemake.output.b1, snakemake.output.b2, snakemake.output.samples]:
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)

Path(snakemake.output.b1).write_text(",".join(case_bams) + "\n")
Path(snakemake.output.b2).write_text(",".join(control_bams) + "\n")

with open(snakemake.output.samples, "w") as handle:
    handle.write("sample\tgroup\trole\tbam\n")
    for sample, bam in zip(case_samples, case_bams):
        handle.write(f"{sample}\t{case_group}\tcase\t{bam}\n")
    for sample, bam in zip(control_samples, control_bams):
        handle.write(f"{sample}\t{control_group}\tcontrol\t{bam}\n")

Path(snakemake.log[0]).parent.mkdir(parents=True, exist_ok=True)
Path(snakemake.log[0]).write_text(
    f"Prepared rMATS pairwise inputs for {snakemake.wildcards.comparison}: "
    f"{case_group} ({len(case_bams)} BAMs) vs {control_group} ({len(control_bams)} BAMs).\n"
)
