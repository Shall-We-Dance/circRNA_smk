rule normalize_fsj_by_uniquely_mapped:
    input:
        fsj_matrix=f"{OUTDIR}/ciri3/all_samples.ciri3.FSJ_Matrix",
        star_logs=expand(
            f"{OUTDIR}/star/{{sample}}/{{sample}}.Log.final.out",
            sample=SAMPLES,
        )
    output:
        normalized=f"{OUTDIR}/ciri3/all_samples.ciri3.FSJ_Matrix.CPM_uniqueMapped.tsv"
    params:
        star_dir=f"{OUTDIR}/star",
        per_million=1e6
    conda:
        "envs/circrna_plots.yaml"
    log:
        "logs/ciri3/normalize_fsj_by_uniquely_mapped.log"
    script:
        "scripts/normalize_fsj_by_uniquely_mapped.R"


rule normalize_totalRNA_by_uniquely_mapped:
    input:
        gene_counts=f"{OUTDIR}/featurecount/totalRNA.counts.txt",
        star_logs=expand(
            f"{OUTDIR}/star/{{sample}}/{{sample}}.Log.final.out",
            sample=SAMPLES,
        )
    output:
        cleaned=f"{OUTDIR}/featurecount/totalRNA.counts.cleaned.tsv",
        normalized=f"{OUTDIR}/featurecount/totalRNA.counts.CPM_uniqueMapped.tsv"
    params:
        star_dir=f"{OUTDIR}/star",
        per_million=1e6
    conda:
        "envs/circrna_plots.yaml"
    log:
        "logs/featurecount/normalize_totalRNA_by_uniquely_mapped.log"
    script:
        "scripts/normalize_totalRNA_by_uniquely_mapped.R"
