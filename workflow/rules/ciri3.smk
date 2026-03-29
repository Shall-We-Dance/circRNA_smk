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


rule ciri3_samples_tsv_star:
    input:
        bams=expand(f"{OUTDIR}/tmp/star/{{sample}}/{{sample}}.Aligned.sortedByCoord.out.bam", sample=SAMPLES)
    output:
        tsv=f"{OUTDIR}/ciri3/samples.star.tsv"
    run:
        os.makedirs(os.path.dirname(output.tsv), exist_ok=True)
        with open(output.tsv, "w", encoding="utf-8") as fout:
            for bam in input.bams:
                sample = os.path.basename(bam).replace(".Aligned.sortedByCoord.out.bam", "")
                fout.write(f"{sample}\t{bam}\n")


rule ciri3_samples_tsv_unique:
    input:
        bams=expand(f"{OUTDIR}/star/{{sample}}/{{sample}}.unique.mapq11.sorted.bam", sample=SAMPLES)
    output:
        tsv=f"{OUTDIR}/ciri3/samples.unique.tsv"
    run:
        os.makedirs(os.path.dirname(output.tsv), exist_ok=True)
        with open(output.tsv, "w", encoding="utf-8") as fout:
            for bam in input.bams:
                sample = os.path.basename(bam).replace(".unique.mapq11.sorted.bam", "")
                fout.write(f"{sample}\t{bam}\n")


rule ciri3_detect_star:
    input:
        tsv=f"{OUTDIR}/ciri3/samples.star.tsv",
        fasta=config["reference"]["fasta"],
        gtf=config["reference"]["gtf"],
        bais=expand(f"{OUTDIR}/tmp/star/{{sample}}/{{sample}}.Aligned.sortedByCoord.out.bam.bai", sample=SAMPLES)
    output:
        result=f"{OUTDIR}/ciri3/star.all_samples.ciri3",
        bsj=f"{OUTDIR}/ciri3/star.all_samples.ciri3.BSJ_Matrix",
        fsj=f"{OUTDIR}/ciri3/star.all_samples.ciri3.FSJ_Matrix"
    log:
        "logs/ciri3.detect.star.log"
    conda:
        "envs/ciri3.yaml"
    params:
        jar=config["ciri3"]["jar"],
        outprefix=f"{OUTDIR}/ciri3/star.all_samples.ciri3",
        Ma=int(config.get("ciri3", {}).get("Ma", 1)),
        W=int(config.get("ciri3", {}).get("W", 1))
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.result}) $(dirname {log})
        java -jar {params.jar} \
          -A {input.gtf} \
          -Ma {params.Ma} \
          -W {params.W} \
          -I {input.tsv} \
          -O {params.outprefix} \
          -F {input.fasta} \
          > {log} 2>&1

        # CIRI3 naming differs across releases; normalize to expected targets.
        if [ ! -e {output.result} ] && [ -e {params.outprefix}.ciri3 ]; then cp {params.outprefix}.ciri3 {output.result}; fi
        if [ ! -e {output.bsj} ] && [ -e {params.outprefix}.ciri3.BSJ_Matrix ]; then cp {params.outprefix}.ciri3.BSJ_Matrix {output.bsj}; fi
        if [ ! -e {output.fsj} ] && [ -e {params.outprefix}.ciri3.FSJ_Matrix ]; then cp {params.outprefix}.ciri3.FSJ_Matrix {output.fsj}; fi

        test -s {output.result}
        test -s {output.bsj}
        test -s {output.fsj}
        """


rule ciri3_detect_unique:
    input:
        tsv=f"{OUTDIR}/ciri3/samples.unique.tsv",
        fasta=config["reference"]["fasta"],
        gtf=config["reference"]["gtf"],
        bais=expand(f"{OUTDIR}/star/{{sample}}/{{sample}}.unique.mapq11.sorted.bam.bai", sample=SAMPLES)
    output:
        result=f"{OUTDIR}/ciri3/unique.all_samples.ciri3",
        bsj=f"{OUTDIR}/ciri3/unique.all_samples.ciri3.BSJ_Matrix",
        fsj=f"{OUTDIR}/ciri3/unique.all_samples.ciri3.FSJ_Matrix"
    log:
        "logs/ciri3.detect.unique.log"
    conda:
        "envs/ciri3.yaml"
    params:
        jar=config["ciri3"]["jar"],
        outprefix=f"{OUTDIR}/ciri3/unique.all_samples.ciri3",
        Ma=int(config.get("ciri3", {}).get("Ma", 1)),
        W=int(config.get("ciri3", {}).get("W", 1))
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.result}) $(dirname {log})
        java -jar {params.jar} \
          -A {input.gtf} \
          -Ma {params.Ma} \
          -W {params.W} \
          -I {input.tsv} \
          -O {params.outprefix} \
          -F {input.fasta} \
          > {log} 2>&1

        # CIRI3 naming differs across releases; normalize to expected targets.
        if [ ! -e {output.result} ] && [ -e {params.outprefix}.ciri3 ]; then cp {params.outprefix}.ciri3 {output.result}; fi
        if [ ! -e {output.bsj} ] && [ -e {params.outprefix}.ciri3.BSJ_Matrix ]; then cp {params.outprefix}.ciri3.BSJ_Matrix {output.bsj}; fi
        if [ ! -e {output.fsj} ] && [ -e {params.outprefix}.ciri3.FSJ_Matrix ]; then cp {params.outprefix}.ciri3.FSJ_Matrix {output.fsj}; fi

        test -s {output.result}
        test -s {output.bsj}
        test -s {output.fsj}
        """
