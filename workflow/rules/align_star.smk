# workflow/rules/align_star.smk
OUTDIR = config["output"]["dir"]
BWA_INDEX_EXTS = ("amb", "ann", "bwt", "pac", "sa")
BWA_INDEXED_FASTA = config["reference"]["bwa_indexed_fasta"]


def maybe_temp(path):
    return path if KEEP_BAM else temp(path)



rule star_align_ciri3:
    input:
        idx_sa=f"{config['reference']['star_index']}/SA",
        r1=f"{OUTDIR}/tmp/fastp/{{sample}}_R1.fastq.gz",
        r2=f"{OUTDIR}/tmp/fastp/{{sample}}_R2.fastq.gz"
    output:
        chimeric=f"{OUTDIR}/star/{{sample}}/{{sample}}.Chimeric.out.junction",
        bam=maybe_temp(f"{OUTDIR}/star/{{sample}}/{{sample}}.Aligned.sortedByCoord.out.bam"),
        unmapped_r1=temp(f"{OUTDIR}/star/{{sample}}/{{sample}}.Unmapped.out.mate1"),
        unmapped_r2=temp(f"{OUTDIR}/star/{{sample}}/{{sample}}.Unmapped.out.mate2"),
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
        mkdir -p $(dirname {output.bam}) $(dirname {output.log_final}) $(dirname {output.log_final_qc}) $(dirname {log})

        STAR \
          --runThreadN {threads} \
          --genomeDir {params.index} \
          --readFilesIn {input.r1} {input.r2} \
          --readFilesCommand zcat \
          --outFileNamePrefix {OUTDIR}/star/{wildcards.sample}/{wildcards.sample}. \
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

        cp {output.log_final} {output.log_final_qc}
        """


rule bwa_remap_ciri3:
    input:
        bam=f"{OUTDIR}/star/{{sample}}/{{sample}}.Aligned.sortedByCoord.out.bam",
        unmapped_r1=f"{OUTDIR}/star/{{sample}}/{{sample}}.Unmapped.out.mate1",
        unmapped_r2=f"{OUTDIR}/star/{{sample}}/{{sample}}.Unmapped.out.mate2",
        bwa_idx=lambda wildcards: [f"{BWA_INDEXED_FASTA}.{ext}" for ext in BWA_INDEX_EXTS]
    output:
        bwa_bam=maybe_temp(f"{OUTDIR}/star/{{sample}}/{{sample}}.bwa.bam")
    log:
        f"logs/star/{{sample}}.bwa.log"
    threads: int(config["threads"]["star"])
    conda:
        "envs/star.yaml"
    params:
        bwa_index=BWA_INDEXED_FASTA
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.bwa_bam}) $(dirname {log})

        # STAR mode of CIRI3 requires BWA remapping for unmapped reads.
        # If both unmapped FASTQ files are empty, emit an empty BAM with the
        # same header as STAR's coordinate-sorted BAM.
        if [[ ! -s {input.unmapped_r1} && ! -s {input.unmapped_r2} ]]; then
          samtools view -H {input.bam} \
            | samtools view -@ {threads} -b - \
            > {output.bwa_bam} 2>> {log}
        else
          bwa mem -T 19 -t {threads} {params.bwa_index} {input.unmapped_r1} {input.unmapped_r2} 2>> {log} \
            | samtools view -@ {threads} -b - \
            > {output.bwa_bam} 2>> {log}
        fi
        """
