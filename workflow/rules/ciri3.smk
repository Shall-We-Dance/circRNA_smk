# workflow/rules/ciri3.smk
import os

OUTDIR = config["output"]["dir"]
SAMPLES = list(config["samples"].keys())


rule star_bam_index_for_ciri3:
    input:
        bam=f"{OUTDIR}/tmp/star/{{sample}}/{{sample}}.Aligned.sortedByCoord.out.bam"
    output:
        bai=f"{OUTDIR}/tmp/star/{{sample}}/{{sample}}.Aligned.sortedByCoord.out.bam.bai"
    threads: int(config["threads"]["samtools"])
    conda:
        "envs/samtools.yaml"
    shell:
        r"""
        set -euo pipefail
        samtools index -@ {threads} {input.bam}
        """


rule featurecounts_totalrna:
    input:
        bams=expand(f"{OUTDIR}/tmp/star/{{sample}}/{{sample}}.Aligned.sortedByCoord.out.bam", sample=SAMPLES)
    output:
        counts=f"{OUTDIR}/ciri_star/totalRNA.counts.txt"
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
        mkdir -p $(dirname {output.counts}) $(dirname {log})
        featureCounts -T {threads} -p -a {params.gtf} -o {output.counts} {input.bams} > {log} 2>&1
        """


rule ciri3_samples_tsv:
    input:
        bams=expand(f"{OUTDIR}/tmp/star/{{sample}}/{{sample}}.Aligned.sortedByCoord.out.bam", sample=SAMPLES)
    output:
        tsv=f"{OUTDIR}/ciri3/samples.tsv"
    run:
        os.makedirs(os.path.dirname(output.tsv), exist_ok=True)
        with open(output.tsv, "w", encoding="utf-8") as fout:
            for bam in input.bams:
                sample = os.path.basename(bam).replace(".Aligned.sortedByCoord.out.bam", "")
                fout.write(f"{sample}\t{bam}\n")


rule ciri3_detect:
    input:
        tsv=f"{OUTDIR}/ciri3/samples.tsv",
        fasta=config["reference"]["fasta"],
        gtf=config["reference"]["gtf"],
        bais=expand(f"{OUTDIR}/tmp/star/{{sample}}/{{sample}}.Aligned.sortedByCoord.out.bam.bai", sample=SAMPLES)
    output:
        result=f"{OUTDIR}/ciri3/all_samples.ciri3",
        bsj=f"{OUTDIR}/ciri3/all_samples.ciri3.BSJ_Matrix",
        fsj=f"{OUTDIR}/ciri3/all_samples.ciri3.FSJ_Matrix"
    log:
        "logs/ciri3.detect.log"
    conda:
        "envs/ciri3.yaml"
    params:
        jar=config["ciri3"]["jar"],
        outprefix=f"{OUTDIR}/ciri3/all_samples.ciri3",
        ma=int(config.get("ciri3", {}).get("ma", 1)),
        w=int(config.get("ciri3", {}).get("w", 1)),
        a=("-A" if bool(config.get("ciri3", {}).get("a", True)) else "")
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.result}) $(dirname {log})
        java -jar {params.jar} \
          -G {input.gtf} \
          -Ma {params.ma} \
          -W {params.w} \
          -I {input.tsv} \
          -O {params.outprefix} \
          -F {input.fasta} \
          {params.a} \
          > {log} 2>&1
        """
