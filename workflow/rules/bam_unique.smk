# workflow/rules/bam_unique.smk
import os

OUTDIR = config["output"]["dir"]
MIN_MAPQ = int(config["filtering"]["min_mapq"])


rule filter_unique_primary_and_index_star:
    input:
        bam=f"{OUTDIR}/star/{{sample}}/{{sample}}.Aligned.sortedByCoord.out.bam"
    output:
        bam=f"{OUTDIR}/star/{{sample}}/{{sample}}.unique.mapq11.sorted.bam",
        bai=f"{OUTDIR}/star/{{sample}}/{{sample}}.unique.mapq11.sorted.bam.bai"
    log:
        f"logs/samtools/{{sample}}.mapq_unique_filter.log"
    threads: config["threads"]["samtools"]
    conda:
        "envs/samtools.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.bam}) $(dirname {log})

        # Keep uniquely mapped STAR reads by combining:
        # 1) MAPQ threshold (default: >=11)
        # 2) primary mapped alignments only (-F 0x904 removes unmapped/secondary/supplementary)
        # 3) STAR NH tag equals 1 (exclude multimappers)
        samtools view -@ {threads} -h -q {MIN_MAPQ} -F 0x904 {input.bam} \
          | awk 'BEGIN{{OFS="\t"}} /^@/ {{print; next}} /(^|\t)NH:i:1(\t|$)/ {{print}}' \
          | samtools sort -@ {threads} -o {output.bam} - \
          > {log} 2>&1

        samtools index -@ {threads} {output.bam} >> {log} 2>&1
        """


rule filter_unique_primary_and_index_bwa:
    input:
        bam=f"{OUTDIR}/star/{{sample}}/{{sample}}.bwa.bam"
    output:
        bam=f"{OUTDIR}/star/{{sample}}/{{sample}}.bwa.unique.mapq11.sorted.bam",
        bai=f"{OUTDIR}/star/{{sample}}/{{sample}}.bwa.unique.mapq11.sorted.bam.bai"
    log:
        f"logs/samtools/{{sample}}.bwa.mapq_unique_filter.log"
    threads: config["threads"]["samtools"]
    conda:
        "envs/samtools.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.bam}) $(dirname {log})

        # Keep uniquely mapped BWA reads using MAPQ + primary mapped alignments.
        # BWA output typically has no NH tag, so only mapq/flag filtering is applied.
        samtools view -@ {threads} -h -q {MIN_MAPQ} -F 0x904 {input.bam} \
          | samtools sort -@ {threads} -o {output.bam} - \
          > {log} 2>&1

        samtools index -@ {threads} {output.bam} >> {log} 2>&1
        """
