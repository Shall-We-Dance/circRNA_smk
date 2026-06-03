rule circrna_plots_metadata:
    output:
        metadata=f"{CIRCRNA_PLOTS_OUTDIR}/metadata/sample_metadata.tsv"
    params:
        groups=DEG_GROUPS
    log:
        "logs/circrna_plots/metadata.log"
    script:
        "scripts/generate_circrna_metadata.py"


rule circrna_plots_per_comparison:
    input:
        result=f"{CIRI3_DE_OUTDIR}/{{comparison}}/de_bsj/result.txt",
        infor=f"{CIRI3_DE_OUTDIR}/{{comparison}}/de_bsj/infor.tsv"
    output:
        cleaned=f"{CIRCRNA_PLOTS_OUTDIR}/comparisons/{{comparison}}/cleaned_result.tsv",
        candidates=f"{CIRCRNA_PLOTS_OUTDIR}/comparisons/{{comparison}}/significant_candidates.tsv",
        ma_unlabeled_png=f"{CIRCRNA_PLOTS_OUTDIR}/comparisons/{{comparison}}/MA_unlabeled.png",
        ma_unlabeled_pdf=f"{CIRCRNA_PLOTS_OUTDIR}/comparisons/{{comparison}}/MA_unlabeled.pdf",
        ma_labeled_png=f"{CIRCRNA_PLOTS_OUTDIR}/comparisons/{{comparison}}/MA_labeled.png",
        ma_labeled_pdf=f"{CIRCRNA_PLOTS_OUTDIR}/comparisons/{{comparison}}/MA_labeled.pdf",
        volcano_unlabeled_png=f"{CIRCRNA_PLOTS_OUTDIR}/comparisons/{{comparison}}/volcano_unlabeled.png",
        volcano_unlabeled_pdf=f"{CIRCRNA_PLOTS_OUTDIR}/comparisons/{{comparison}}/volcano_unlabeled.pdf",
        volcano_labeled_png=f"{CIRCRNA_PLOTS_OUTDIR}/comparisons/{{comparison}}/volcano_labeled.png",
        volcano_labeled_pdf=f"{CIRCRNA_PLOTS_OUTDIR}/comparisons/{{comparison}}/volcano_labeled.pdf"
    wildcard_constraints:
        comparison=CIRI3_DE_COMPARISON_REGEX
    params:
        comparison=lambda wc: wc.comparison,
        case_label=lambda wc: DEG_COMPARISONS[wc.comparison]["case_group"],
        control_label=lambda wc: DEG_COMPARISONS[wc.comparison]["control_group"],
        fdr_cutoff=DEG_PADJ_CUTOFF,
        logfc_cutoff=DEG_LFC_CUTOFF,
        label_top_n=CIRCRNA_PLOTS_LABEL_TOP_N,
        reverse_ciri3_logfc_direction=CIRCRNA_PLOTS_REVERSE_LFC,
        fdr_floor_for_plot=CIRCRNA_PLOTS_FDR_FLOOR
    conda:
        "envs/circrna_plots.yaml"
    log:
        "logs/circrna_plots/comparisons/{comparison}.log"
    script:
        "scripts/plot_circrna_de_bsj_per_comparison.R"


rule circrna_plots_dataset_summary:
    input:
        bsj_matrix=f"{OUTDIR}/ciri3/all_samples.ciri3.BSJ_Matrix",
        metadata=f"{CIRCRNA_PLOTS_OUTDIR}/metadata/sample_metadata.tsv",
        star_logs=expand(
            f"{OUTDIR}/star/{{sample}}/{{sample}}.Log.final.out",
            sample=CIRCRNA_PLOTS_SAMPLES,
        )
    output:
        # Part 2: normalization tables and plots
        sample_summary=f"{CIRCRNA_PLOTS_OUTDIR}/dataset/normalization/sample_summary.tsv",
        bsj_cpm_matrix=f"{CIRCRNA_PLOTS_OUTDIR}/dataset/normalization/BSJ_matrix_CPM.tsv",
        condition_summary=f"{CIRCRNA_PLOTS_OUTDIR}/dataset/normalization/condition_summary.tsv",
        sample_barplot_png=f"{CIRCRNA_PLOTS_OUTDIR}/dataset/normalization/all_samples_barplot.png",
        sample_barplot_pdf=f"{CIRCRNA_PLOTS_OUTDIR}/dataset/normalization/all_samples_barplot.pdf",
        conditions_boxplot_png=f"{CIRCRNA_PLOTS_OUTDIR}/dataset/normalization/all_conditions_boxplot.png",
        conditions_boxplot_pdf=f"{CIRCRNA_PLOTS_OUTDIR}/dataset/normalization/all_conditions_boxplot.pdf",
        # Part 3: heatmaps - all circRNAs
        all_cpm_matrix=f"{CIRCRNA_PLOTS_OUTDIR}/dataset/heatmaps/all_circRNAs/absolute_CPM_matrix.tsv",
        all_zscore_matrix=f"{CIRCRNA_PLOTS_OUTDIR}/dataset/heatmaps/all_circRNAs/log2CPM_zscore_matrix.tsv",
        all_cpm_heatmap_png=f"{CIRCRNA_PLOTS_OUTDIR}/dataset/heatmaps/all_circRNAs/absolute_CPM_heatmap.png",
        all_cpm_heatmap_pdf=f"{CIRCRNA_PLOTS_OUTDIR}/dataset/heatmaps/all_circRNAs/absolute_CPM_heatmap.pdf",
        all_zscore_heatmap_png=f"{CIRCRNA_PLOTS_OUTDIR}/dataset/heatmaps/all_circRNAs/zscore_heatmap.png",
        all_zscore_heatmap_pdf=f"{CIRCRNA_PLOTS_OUTDIR}/dataset/heatmaps/all_circRNAs/zscore_heatmap.pdf",
        # Part 3: heatmaps - top-N variable (declared via expand for each top_n in config)
        topN_cpm_matrices=expand(
            f"{CIRCRNA_PLOTS_OUTDIR}/dataset/heatmaps/top{{n}}_variable/absolute_CPM_matrix.tsv",
            n=CIRCRNA_PLOTS_HEATMAP_TOP_N_VALUES,
        ),
        topN_zscore_matrices=expand(
            f"{CIRCRNA_PLOTS_OUTDIR}/dataset/heatmaps/top{{n}}_variable/log2CPM_zscore_matrix.tsv",
            n=CIRCRNA_PLOTS_HEATMAP_TOP_N_VALUES,
        ),
        topN_cpm_heatmaps_png=expand(
            f"{CIRCRNA_PLOTS_OUTDIR}/dataset/heatmaps/top{{n}}_variable/absolute_CPM_heatmap.png",
            n=CIRCRNA_PLOTS_HEATMAP_TOP_N_VALUES,
        ),
        topN_cpm_heatmaps_pdf=expand(
            f"{CIRCRNA_PLOTS_OUTDIR}/dataset/heatmaps/top{{n}}_variable/absolute_CPM_heatmap.pdf",
            n=CIRCRNA_PLOTS_HEATMAP_TOP_N_VALUES,
        ),
        topN_zscore_heatmaps_png=expand(
            f"{CIRCRNA_PLOTS_OUTDIR}/dataset/heatmaps/top{{n}}_variable/zscore_heatmap.png",
            n=CIRCRNA_PLOTS_HEATMAP_TOP_N_VALUES,
        ),
        topN_zscore_heatmaps_pdf=expand(
            f"{CIRCRNA_PLOTS_OUTDIR}/dataset/heatmaps/top{{n}}_variable/zscore_heatmap.pdf",
            n=CIRCRNA_PLOTS_HEATMAP_TOP_N_VALUES,
        )
    params:
        star_dir=f"{OUTDIR}/star",
        condition_order=CIRCRNA_PLOTS_CONDITION_ORDER,
        condition_colors=CIRCRNA_PLOTS_CONDITION_COLORS_RESOLVED,
        per_million=CIRCRNA_PLOTS_PER_MILLION,
        heatmap_top_n_values=CIRCRNA_PLOTS_HEATMAP_TOP_N_VALUES,
        heatmap_cluster_rows=CIRCRNA_PLOTS_HEATMAP_CLUSTER_ROWS,
        heatmap_cluster_cols=CIRCRNA_PLOTS_HEATMAP_CLUSTER_COLS,
        heatmap_show_rownames_all=CIRCRNA_PLOTS_HEATMAP_SHOW_ROWNAMES_ALL,
        heatmap_show_rownames_top=CIRCRNA_PLOTS_HEATMAP_SHOW_ROWNAMES_TOP,
        heatmap_absolute_upper_clip_quantile=CIRCRNA_PLOTS_HEATMAP_CLIP_Q,
        heatmap_absolute_use_log_for_colormap=CIRCRNA_PLOTS_HEATMAP_USE_LOG_COLORMAP,
        heatmap_log_pseudocount=CIRCRNA_PLOTS_HEATMAP_LOG_PSEUDOCOUNT,
        heatmap_zscore_breaks_limit=CIRCRNA_PLOTS_HEATMAP_ZSCORE_BREAKS,
        topN_outdir_template=lambda wildcards: f"{CIRCRNA_PLOTS_OUTDIR}/dataset/heatmaps/top{{n}}_variable"
    conda:
        "envs/circrna_plots.yaml"
    log:
        "logs/circrna_plots/dataset_summary.log"
    script:
        "scripts/plot_circrna_dataset_summary.R"


rule circrna_plots_linear_counterpart_summary:
    input:
        fsj_cpm_matrix=f"{OUTDIR}/ciri3/all_samples.ciri3.FSJ_Matrix.CPM_uniqueMapped.tsv",
        totalrna_cpm_matrix=f"{OUTDIR}/featurecount/totalRNA.counts.CPM_uniqueMapped.tsv",
        metadata=f"{CIRCRNA_PLOTS_OUTDIR}/metadata/sample_metadata.tsv"
    output:
        fsj_barplot_png=f"{CIRCRNA_PLOTS_OUTDIR}/linear_counterpart/fsj/all_samples_barplot.png",
        fsj_barplot_pdf=f"{CIRCRNA_PLOTS_OUTDIR}/linear_counterpart/fsj/all_samples_barplot.pdf",
        fsj_boxplot_png=f"{CIRCRNA_PLOTS_OUTDIR}/linear_counterpart/fsj/all_conditions_boxplot.png",
        fsj_boxplot_pdf=f"{CIRCRNA_PLOTS_OUTDIR}/linear_counterpart/fsj/all_conditions_boxplot.pdf",
        totalrna_barplot_png=f"{CIRCRNA_PLOTS_OUTDIR}/linear_counterpart/totalRNA/all_samples_barplot.png",
        totalrna_barplot_pdf=f"{CIRCRNA_PLOTS_OUTDIR}/linear_counterpart/totalRNA/all_samples_barplot.pdf",
        totalrna_boxplot_png=f"{CIRCRNA_PLOTS_OUTDIR}/linear_counterpart/totalRNA/all_conditions_boxplot.png",
        totalrna_boxplot_pdf=f"{CIRCRNA_PLOTS_OUTDIR}/linear_counterpart/totalRNA/all_conditions_boxplot.pdf"
    params:
        condition_order=CIRCRNA_PLOTS_CONDITION_ORDER,
        condition_colors=CIRCRNA_PLOTS_CONDITION_COLORS_RESOLVED
    conda:
        "envs/circrna_plots.yaml"
    log:
        "logs/circrna_plots/linear_counterpart_summary.log"
    script:
        "scripts/plot_circrna_linear_counterpart.R"
