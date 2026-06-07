rule export_rnaseq_bigwig_tracks:
    input:
        bam=lambda wc: star_bam_path(wc.sample),
        bai=lambda wc: star_bai_path(wc.sample)
    output:
        raw=f"{OUTDIR}/igv/rnaseq/raw/{{sample}}.rnaseq.raw.bw",
        normalized=f"{OUTDIR}/igv/rnaseq/normalized/{{sample}}.rnaseq.CPM.bw"
    log:
        "logs/igv/rnaseq/{sample}.bigwig.log"
    threads: int(config["threads"].get("igv_bigwig", 2))
    conda:
        "envs/deeptools.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.raw}")" "$(dirname "{output.normalized}")" "$(dirname "{log}")"

        echo "Sample: {wildcards.sample}" > "{log}"
        echo "Input BAM: {input.bam}" >> "{log}"
        echo "Raw BigWig: {output.raw}" >> "{log}"
        echo "Normalized BigWig: {output.normalized}" >> "{log}"
        echo "bamCoverage binSize=1 samFlagExclude=3844" >> "{log}"

        bamCoverage \
          --bam "{input.bam}" \
          --outFileName "{output.raw}" \
          --outFileFormat bigwig \
          --binSize 1 \
          --normalizeUsing None \
          --samFlagExclude 3844 \
          --numberOfProcessors {threads} \
          >> "{log}" 2>&1

        bamCoverage \
          --bam "{input.bam}" \
          --outFileName "{output.normalized}" \
          --outFileFormat bigwig \
          --binSize 1 \
          --normalizeUsing CPM \
          --samFlagExclude 3844 \
          --numberOfProcessors {threads} \
          >> "{log}" 2>&1
        """
