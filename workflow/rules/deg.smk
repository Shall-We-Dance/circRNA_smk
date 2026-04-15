OUTDIR = config["output"]["dir"]


rule bsj_deg_analysis:
    input:
        bsj=f"{OUTDIR}/ciri3/all_samples.ciri3.BSJ_Matrix"
    output:
        metadata=f"{OUTDIR}/deg/bsj/sample_metadata.tsv",
        all_results=f"{OUTDIR}/deg/bsj/all_groups/deseq2_results.tsv",
        all_volcano=f"{OUTDIR}/deg/bsj/all_groups/volcano.pdf",
        all_heatmap=f"{OUTDIR}/deg/bsj/all_groups/heatmap_top50.pdf",
        pca=f"{OUTDIR}/deg/bsj/all_groups/pca.pdf",
        vst_counts=f"{OUTDIR}/deg/bsj/all_groups/vst_counts.tsv",
        pairwise_results=expand(
            f"{OUTDIR}/deg/bsj/pairwise/{{pair}}/deseq2_results.tsv",
            pair=[f"{a}_vs_{b}" for a, b in DEG_PAIRWISE],
        ),
        pairwise_volcano=expand(
            f"{OUTDIR}/deg/bsj/pairwise/{{pair}}/volcano.pdf",
            pair=[f"{a}_vs_{b}" for a, b in DEG_PAIRWISE],
        ),
        pairwise_heatmap=expand(
            f"{OUTDIR}/deg/bsj/pairwise/{{pair}}/heatmap_top50.pdf",
            pair=[f"{a}_vs_{b}" for a, b in DEG_PAIRWISE],
        )
    params:
        groups=DEG_GROUPS
    conda:
        "envs/deg.yaml"
    log:
        "logs/deg/bsj_deg_analysis.log"
    script:
        "scripts/run_bsj_deg.R"
