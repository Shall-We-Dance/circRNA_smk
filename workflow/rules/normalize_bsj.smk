rule normalize_bsj_by_uniquely_mapped:
    input:
        bsj_matrix=f"{OUTDIR}/ciri3/all_samples.ciri3.BSJ_Matrix",
        star_logs=expand(
            f"{OUTDIR}/star/{{sample}}/{{sample}}.Log.final.out",
            sample=SAMPLES,
        )
    output:
        normalized=f"{OUTDIR}/ciri3/all_samples.ciri3.BSJ_Matrix.CPM_uniqueMapped.tsv"
    params:
        star_dir=f"{OUTDIR}/star",
        per_million=1e6
    conda:
        "envs/circrna_plots.yaml"
    log:
        "logs/ciri3/normalize_bsj_by_uniquely_mapped.log"
    script:
        "scripts/normalize_bsj_by_uniquely_mapped.R"
