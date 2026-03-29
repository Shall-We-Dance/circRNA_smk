# workflow/rules/align_star.smk
import os

OUTDIR = config["output"]["dir"]



rule star_align_ciri3:
    input:
        idx_sa=f"{config['reference']['star_index']}/SA",
        r1=f"{OUTDIR}/tmp/fastp/{{sample}}_R1.fastq.gz",
        r2=f"{OUTDIR}/tmp/fastp/{{sample}}_R2.fastq.gz"
    output:
        bam=temp(f"{OUTDIR}/tmp/star/{{sample}}/{{sample}}.Aligned.sortedByCoord.out.bam"),
        unmapped_r1=temp(f"{OUTDIR}/tmp/star/{{sample}}/{{sample}}.Unmapped.out.mate1"),
        unmapped_r2=temp(f"{OUTDIR}/tmp/star/{{sample}}/{{sample}}.Unmapped.out.mate2"),
        log_final=f"{OUTDIR}/star/{{sample}}/{{sample}}.Log.final.out",
        log_final_qc=f"{OUTDIR}/qc/star/{{sample}}/{{sample}}.Log.final.out",
        sj=f"{OUTDIR}/star/{{sample}}/{{sample}}.SJ.out.tab"
    log:
        f"logs/star/{{sample}}.log"
    threads: int(config["threads"]["star"])
    conda:
        "envs/star.yaml"
    params:
        index=config["reference"]["star_index"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.bam}) $(dirname {output.log_final}) $(dirname {log})

        STAR \
          --runThreadN {threads} \
          --genomeDir {params.index} \
          --readFilesIn {input.r1} {input.r2} \
          --readFilesCommand zcat \
          --outFileNamePrefix {OUTDIR}/tmp/star/{wildcards.sample}/{wildcards.sample}. \
          --outSAMtype BAM SortedByCoordinate \
          --outReadsUnmapped Fastx \
          --outSJfilterOverhangMin 15 12 12 12 \
          --alignSJoverhangMin 15 \
          --alignSJDBoverhangMin 15 \
          --outFilterMultimapNmax 20 \
          --outFilterScoreMin 1 \
          --outFilterMatchNmin 1 \
          --outFilterMismatchNmax 2 \
          --chimSegmentMin 15 \
          --chimScoreMin 15 \
          --chimJunctionOverhangMin 15 \
          > {log} 2>&1

        cp {OUTDIR}/tmp/star/{wildcards.sample}/{wildcards.sample}.Log.final.out {output.log_final}
        cp {OUTDIR}/tmp/star/{wildcards.sample}/{wildcards.sample}.Log.final.out {output.log_final_qc}
        cp {OUTDIR}/tmp/star/{wildcards.sample}/{wildcards.sample}.SJ.out.tab {output.sj}
        """
