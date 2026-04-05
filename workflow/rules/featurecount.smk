# workflow/rules/featurecount.smk
OUTDIR = config["output"]["dir"]
SAMPLES = list(config["samples"].keys())


rule featurecounts_totalrna:
    input:
        bams=expand(f"{OUTDIR}/star/{{sample}}/{{sample}}.Aligned.sortedByCoord.out.bam", sample=SAMPLES)
    output:
        counts=f"{OUTDIR}/featurecount/totalRNA.counts.txt"
    log:
        "logs/featureCounts.log"
    threads: int(config["threads"].get("featurecounts", config["threads"]["star"]))
    conda:
        "envs/ciri3.yaml"
    params:
        gtf=config["reference"]["gtf"],
        sample_names=",".join(SAMPLES)
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.counts}) $(dirname {log})
        featureCounts -T {threads} -p -a {params.gtf} -o {output.counts} {input.bams} > {log} 2>&1

        python - <<'PY'
from pathlib import Path

counts_path = Path("{output.counts}")
sample_names = "{params.sample_names}".split(",")

lines = counts_path.read_text().splitlines()
if len(lines) < 2:
    raise ValueError("featureCounts output is malformed: {}".format(counts_path))

header = lines[1].split("\t")
fixed_columns = 6
if len(header) - fixed_columns != len(sample_names):
    raise ValueError(
        "featureCounts sample columns count does not match configured samples: "
        "{} vs {}".format(len(header) - fixed_columns, len(sample_names))
    )

header[fixed_columns:] = sample_names
lines[1] = "\t".join(header)
counts_path.write_text("\n".join(lines) + "\n")
PY
        """
