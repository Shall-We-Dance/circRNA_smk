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

        python workflow/scripts/rename_featurecounts_header.py \
            --counts {output.counts} \
            --sample-names {params.sample_names}
        """
