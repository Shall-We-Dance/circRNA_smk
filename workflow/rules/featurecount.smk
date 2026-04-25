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
        gtf=config["reference"]["gtf"]
    shell:
        r"""
        set -euo pipefail
        output_dir=$(dirname "{output.counts}")
        log_dir=$(dirname "{log}")
        mkdir -p "$output_dir" "$log_dir"
        featureCounts -T {threads} -p -a "{params.gtf}" -o "{output.counts}" {input.bams:q} > "{log}" 2>&1
        """
