OUTDIR = config["output"]["dir"]


rule bsj_deg_analysis:
    input:
        bsj=f"{OUTDIR}/ciri3/all_samples.ciri3.BSJ_Matrix",
        ciri3=f"{OUTDIR}/ciri3/all_samples.ciri3"
    output:
        metadata=f"{OUTDIR}/deg/bsj/sample_metadata.tsv",
        all_results=f"{OUTDIR}/deg/bsj/all_groups/deseq2_results.tsv",
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
        groups=DEG_GROUPS,
        min_total_count=DEG_MIN_TOTAL_COUNT,
        min_samples_detected=DEG_MIN_SAMPLES_DETECTED,
        padj_cutoff=DEG_PADJ_CUTOFF,
        lfc_cutoff=DEG_LFC_CUTOFF
    conda:
        "envs/deg.yaml"
    log:
        "logs/deg/bsj_deg_analysis.log"
    script:
        "scripts/run_bsj_deg.R"
