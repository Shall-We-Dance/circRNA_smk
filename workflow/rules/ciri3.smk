# workflow/rules/ciri3.smk

rule star_bam_index_for_ciri3:
    input:
        bam=f"{OUTDIR}/star/{{sample}}/{{sample}}.Aligned.sortedByCoord.out.bam"
    output:
        bai=maybe_temp(f"{OUTDIR}/star/{{sample}}/{{sample}}.Aligned.sortedByCoord.out.bam.bai")
    log:
        "logs/ciri3/{sample}.star_bam_index.log"
    threads: int(config["threads"].get("samtools_index", config["threads"].get("samtools", 1)))
    conda:
        "envs/samtools.yaml"
    shell:
        r"""
        set -euo pipefail
        log_dir=$(dirname "{log}")
        mkdir -p "$log_dir"
        samtools index -@ {threads} "{input.bam}" > "{log}" 2>&1
        """


rule ciri3_detect_star:
    input:
        chimeric=f"{OUTDIR}/star/{{sample}}/{{sample}}.Chimeric.out.junction",
        bam=f"{OUTDIR}/star/{{sample}}/{{sample}}.Aligned.sortedByCoord.out.bam",
        bwa=f"{OUTDIR}/star/{{sample}}/{{sample}}.bwa.bam",
        fasta=config["reference"]["fasta"],
        gtf=config["reference"]["gtf"],
        bai=f"{OUTDIR}/star/{{sample}}/{{sample}}.Aligned.sortedByCoord.out.bam.bai"
    output:
        result=f"{OUTDIR}/ciri3/per_sample/{{sample}}.ciri3",
        bsj=f"{OUTDIR}/ciri3/per_sample/{{sample}}.ciri3.BSJ_Matrix",
        fsj=f"{OUTDIR}/ciri3/per_sample/{{sample}}.ciri3.FSJ_Matrix"
    log:
        "logs/ciri3/{sample}.detect.star.log"
    threads: int(config["threads"].get("samtools_view", config["threads"].get("samtools", 1)))
    conda:
        "envs/ciri3.yaml"
    params:
        jar=config["ciri3"]["jar"],
        Ma=int(config.get("ciri3", {}).get("Ma", 1)),
        W=int(config.get("ciri3", {}).get("W", 1)),
        normalize_script=CIRI3_NORMALIZE_SCRIPT
    shell:
        r"""
        set -euo pipefail
        output_dir=$(dirname "{output.result}")
        log_dir=$(dirname "{log}")
        mkdir -p "$output_dir" "$log_dir"
        tmpdir=$(mktemp -d "$output_dir/ciri3_{wildcards.sample}.XXXXXX")
        trap 'rm -rf "$tmpdir"' EXIT

        star_sam="$tmpdir/{wildcards.sample}.Aligned.sortedByCoord.out.sam"
        bwa_sam="$tmpdir/{wildcards.sample}.bwa.sam"
        samples_tsv="$tmpdir/samples_{wildcards.sample}.tsv"

        samtools view -@ {threads} -h -o "$star_sam" "{input.bam}"
        samtools view -@ {threads} -h -o "$bwa_sam" "{input.bwa}"

        chimeric_abs=$(realpath "{input.chimeric}")
        star_sam_abs=$(realpath "$star_sam")
        bwa_sam_abs=$(realpath "$bwa_sam")
        printf "%s,%s,%s\n" "$chimeric_abs" "$star_sam_abs" "$bwa_sam_abs" > "$samples_tsv"

        java -jar "{params.jar}" \
          -A "{input.gtf}" \
          -Ma {params.Ma} \
          -W {params.W} \
          -I "$samples_tsv" \
          -O "{output.result}" \
          -F "{input.fasta}" \
          > "{log}" 2>&1

        python "{params.normalize_script}" \
          --sample "{wildcards.sample}" \
          --result "{output.result}" \
          --bsj "{output.bsj}" \
          --fsj "{output.fsj}" \
          --tmpdir "$tmpdir" \
          --log "{log}"

        test -s "{output.result}"
        test -s "{output.bsj}"
        test -s "{output.fsj}"
        """


rule merge_ciri3_outputs:
    input:
        ciri3=expand(f"{OUTDIR}/ciri3/per_sample/{{sample}}.ciri3", sample=SAMPLES),
        bsj=expand(f"{OUTDIR}/ciri3/per_sample/{{sample}}.ciri3.BSJ_Matrix", sample=SAMPLES),
        fsj=expand(f"{OUTDIR}/ciri3/per_sample/{{sample}}.ciri3.FSJ_Matrix", sample=SAMPLES)
    output:
        ciri3=f"{OUTDIR}/ciri3/all_samples.ciri3",
        bsj=f"{OUTDIR}/ciri3/all_samples.ciri3.BSJ_Matrix",
        fsj=f"{OUTDIR}/ciri3/all_samples.ciri3.FSJ_Matrix"
    log:
        "logs/ciri3/merge_outputs.log"
    params:
        samples=SAMPLES
    conda:
        "envs/ciri3.yaml"
    script:
        "scripts/merge_ciri3_outputs.py"
