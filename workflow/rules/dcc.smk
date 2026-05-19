# workflow/rules/dcc.smk

rule dcc_detect_star:
    input:
        chimeric=expand(f"{OUTDIR}/star/{{sample}}/{{sample}}.Chimeric.out.junction", sample=SAMPLES),
        bam=expand(f"{OUTDIR}/star/{{sample}}/{{sample}}.Aligned.sortedByCoord.out.bam", sample=SAMPLES),
        bai=expand(f"{OUTDIR}/star/{{sample}}/{{sample}}.Aligned.sortedByCoord.out.bam.bai", sample=SAMPLES),
        sj=expand(f"{OUTDIR}/star/{{sample}}/{{sample}}.SJ.out.tab", sample=SAMPLES),
        fasta=config["reference"]["fasta"],
        gtf=config["reference"]["gtf"],
        dcc_runner=DCC_RUNNER_SCRIPT
    output:
        circ_counts=f"{DCC_RAW_OUTDIR}/CircRNACount",
        circ_coordinates=f"{DCC_RAW_OUTDIR}/CircCoordinates",
        linear_counts=f"{DCC_RAW_OUTDIR}/LinearCount",
        circ_skip=f"{DCC_RAW_OUTDIR}/CircSkipJunctions",
        annotation=f"{DCC_OUTDIR}/all_samples.dcc",
        bsj=f"{DCC_OUTDIR}/all_samples.dcc.BSJ_Matrix",
        fsj=f"{DCC_OUTDIR}/all_samples.dcc.FSJ_Matrix",
        per_sample_annotation=expand(f"{DCC_OUTDIR}/per_sample/{{sample}}.dcc", sample=SAMPLES),
        per_sample_bsj=expand(f"{DCC_OUTDIR}/per_sample/{{sample}}.dcc.BSJ_Matrix", sample=SAMPLES),
        per_sample_fsj=expand(f"{DCC_OUTDIR}/per_sample/{{sample}}.dcc.FSJ_Matrix", sample=SAMPLES)
    log:
        "logs/dcc/dcc_detect_star.log"
    threads: int(config["threads"].get("dcc", config["threads"].get("star", 1)))
    conda:
        "envs/dcc.yaml"
    params:
        samples=SAMPLES,
        normalize_script=DCC_NORMALIZE_SCRIPT,
        stranded_flag="" if DCC_STRANDED else "-N",
        secondstrand_flag="-ss" if DCC_SECONDSTRAND else "",
        filter_flag="-F" if DCC_FILTER else "",
        chrM_flag="-M" if DCC_FILTER_CHRM else "",
        filter_gene_flag="-fg" if DCC_FILTER_BY_GENE else "",
        keep_temp_flag="-k" if DCC_KEEP_TEMP else "",
        gene_flag="-G" if DCC_RUN_GENE_COUNTS else "",
        repeat_file=DCC_REPEAT_GTF,
        end_tolerance=DCC_END_TOL,
        min_length=DCC_MIN_LENGTH,
        max_length=DCC_MAX_LENGTH,
        min_count=DCC_MIN_COUNT,
        min_replicates=DCC_MIN_REPLICATES,
        extra_args=DCC_EXTRA_ARGS
    shell:
        r"""
        set -euo pipefail
        raw_dir="{DCC_RAW_OUTDIR}"
        out_dir="{DCC_OUTDIR}"
        tmp_dir="$out_dir/tmp"
        input_dir="$out_dir/inputs"
        log_dir=$(dirname "{log}")
        mkdir -p "$raw_dir" "$tmp_dir" "$input_dir" "$log_dir"

        samplesheet="$input_dir/samplesheet.txt"
        bam_list="$input_dir/bams.txt"
        : > "$samplesheet"
        : > "$bam_list"
        repeat_arg=()
        if [ -n "{params.repeat_file}" ]; then
          repeat_arg=(-R "{params.repeat_file}")
        fi

        for path in {input.chimeric:q}; do
          realpath "$path" >> "$samplesheet"
        done
        for path in {input.bam:q}; do
          realpath "$path" >> "$bam_list"
        done

        rm -rf "$tmp_dir"
        mkdir -p "$tmp_dir"
        rm -f \
          "{output.circ_counts}" \
          "{output.circ_coordinates}" \
          "{output.linear_counts}" \
          "{output.circ_skip}"

        python "{input.dcc_runner}" "@$samplesheet" \
          -T {threads} \
          -O "$raw_dir" \
          -t "$tmp_dir" \
          -D \
          -an "{input.gtf}" \
          -B "@$bam_list" \
          -A "{input.fasta}" \
          -E {params.end_tolerance} \
          -n {params.min_length} \
          -m {params.max_length} \
          -Nr {params.min_count} {params.min_replicates} \
          {params.gene_flag} \
          {params.stranded_flag} \
          {params.secondstrand_flag} \
          {params.filter_flag} \
          {params.chrM_flag} \
          {params.filter_gene_flag} \
          {params.keep_temp_flag} \
          "${{repeat_arg[@]}}" \
          {params.extra_args} \
          > "{log}" 2>&1

        python "{params.normalize_script}" \
          --samples {params.samples:q} \
          --circ-counts "{output.circ_counts}" \
          --circ-coordinates "{output.circ_coordinates}" \
          --linear-counts "{output.linear_counts}" \
          --circ-skip "{output.circ_skip}" \
          --annotation "{output.annotation}" \
          --bsj "{output.bsj}" \
          --fsj "{output.fsj}" \
          --per-sample-dir "{DCC_OUTDIR}/per_sample" \
          --log "{log}"
        """
