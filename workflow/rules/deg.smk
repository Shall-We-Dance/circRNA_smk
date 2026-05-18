rule bsj_deg_analysis:
    input:
        bsj=f"{OUTDIR}/ciri3/all_samples.ciri3.BSJ_Matrix",
        ciri3=f"{OUTDIR}/ciri3/all_samples.ciri3"
    output:
        metadata=f"{OUTDIR}/deg/deseq2/sample_metadata.tsv",
        all_results=f"{OUTDIR}/deg/deseq2/all_groups/deseq2_results.tsv",
        all_heatmap=f"{OUTDIR}/deg/deseq2/all_groups/heatmap_top50.pdf",
        pca=f"{OUTDIR}/deg/deseq2/all_groups/pca.pdf",
        vst_counts=f"{OUTDIR}/deg/deseq2/all_groups/vst_counts.tsv",
        pairwise_results=expand(
            f"{OUTDIR}/deg/deseq2/pairwise/{{pair}}/deseq2_results.tsv",
            pair=DEG_COMPARISON_NAMES,
        ),
        pairwise_volcano=expand(
            f"{OUTDIR}/deg/deseq2/pairwise/{{pair}}/volcano.pdf",
            pair=DEG_COMPARISON_NAMES,
        ),
        pairwise_volcano_labeled=expand(
            f"{OUTDIR}/deg/deseq2/pairwise/{{pair}}/volcano_labeled.pdf",
            pair=DEG_COMPARISON_NAMES,
        ),
        pairwise_heatmap=expand(
            f"{OUTDIR}/deg/deseq2/pairwise/{{pair}}/heatmap_top50.pdf",
            pair=DEG_COMPARISON_NAMES,
        ),
        pairwise_pca=expand(
            f"{OUTDIR}/deg/deseq2/pairwise/{{pair}}/pca.pdf",
            pair=DEG_COMPARISON_NAMES,
        )
    params:
        groups=DEG_GROUPS,
        comparisons=DEG_COMPARISONS,
        min_total_count=DEG_MIN_TOTAL_COUNT,
        min_samples_detected=DEG_MIN_SAMPLES_DETECTED,
        padj_cutoff=DEG_PADJ_CUTOFF,
        lfc_cutoff=DEG_LFC_CUTOFF
    conda:
        "envs/deg.yaml"
    log:
        "logs/deg/deseq2_analysis.log"
    script:
        "scripts/run_bsj_deg.R"


rule dcc_bsj_deg_analysis:
    input:
        bsj=f"{DCC_OUTDIR}/all_samples.dcc.BSJ_Matrix",
        annotation=f"{DCC_OUTDIR}/all_samples.dcc"
    output:
        metadata=f"{DCC_DEG_OUTDIR}/sample_metadata.tsv",
        all_results=f"{DCC_DEG_OUTDIR}/all_groups/deseq2_results.tsv",
        all_heatmap=f"{DCC_DEG_OUTDIR}/all_groups/heatmap_top50.pdf",
        pca=f"{DCC_DEG_OUTDIR}/all_groups/pca.pdf",
        vst_counts=f"{DCC_DEG_OUTDIR}/all_groups/vst_counts.tsv",
        pairwise_results=expand(
            f"{DCC_DEG_OUTDIR}/pairwise/{{pair}}/deseq2_results.tsv",
            pair=DCC_DEG_COMPARISON_NAMES,
        ),
        pairwise_volcano=expand(
            f"{DCC_DEG_OUTDIR}/pairwise/{{pair}}/volcano.pdf",
            pair=DCC_DEG_COMPARISON_NAMES,
        ),
        pairwise_volcano_labeled=expand(
            f"{DCC_DEG_OUTDIR}/pairwise/{{pair}}/volcano_labeled.pdf",
            pair=DCC_DEG_COMPARISON_NAMES,
        ),
        pairwise_heatmap=expand(
            f"{DCC_DEG_OUTDIR}/pairwise/{{pair}}/heatmap_top50.pdf",
            pair=DCC_DEG_COMPARISON_NAMES,
        ),
        pairwise_pca=expand(
            f"{DCC_DEG_OUTDIR}/pairwise/{{pair}}/pca.pdf",
            pair=DCC_DEG_COMPARISON_NAMES,
        )
    params:
        groups=DEG_GROUPS,
        comparisons=DEG_COMPARISONS,
        min_total_count=DEG_MIN_TOTAL_COUNT,
        min_samples_detected=DEG_MIN_SAMPLES_DETECTED,
        padj_cutoff=DEG_PADJ_CUTOFF,
        lfc_cutoff=DEG_LFC_CUTOFF,
        source_name="DCC"
    conda:
        "envs/deg.yaml"
    log:
        "logs/dcc_deg/deseq2_analysis.log"
    script:
        "scripts/run_bsj_deg.R"


rule dcc_circtest_all_groups:
    input:
        circ_counts=f"{DCC_RAW_OUTDIR}/CircRNACount",
        linear_counts=f"{DCC_RAW_OUTDIR}/LinearCount",
        circ_coordinates=f"{DCC_RAW_OUTDIR}/CircCoordinates"
    output:
        results=f"{DCC_CIRCTEST_OUTDIR}/all_groups/circtest_results.tsv",
        summary=f"{DCC_CIRCTEST_OUTDIR}/all_groups/circtest_summary.tsv",
        filtered_circ=f"{DCC_CIRCTEST_OUTDIR}/all_groups/filtered_circ_counts.tsv",
        filtered_linear=f"{DCC_CIRCTEST_OUTDIR}/all_groups/filtered_linear_counts.tsv",
        ratio_plots=f"{DCC_CIRCTEST_OUTDIR}/all_groups/ratio_plots.pdf"
    params:
        all_samples=SAMPLES,
        selected_samples=DEG_SELECTED_SAMPLES,
        group_labels=[
            DEG_SAMPLE_TO_GROUP[sample]
            for sample in DEG_SELECTED_SAMPLES
        ] if DCC_DEG_RUN_CIRCTEST else [],
        comparison="all_groups",
        case_group="",
        control_group="",
        padj_cutoff=DEG_PADJ_CUTOFF,
        filter_sample=DCC_CIRCTEST_FILTER_SAMPLE,
        filter_count=DCC_CIRCTEST_FILTER_COUNT,
        percentage=DCC_CIRCTEST_PERCENTAGE,
        max_plots=DCC_CIRCTEST_MAX_PLOTS
    conda:
        "envs/dcc.yaml"
    log:
        "logs/dcc_deg/circtest.all_groups.log"
    script:
        "scripts/run_dcc_circtest.R"


rule dcc_circtest_pairwise:
    input:
        circ_counts=f"{DCC_RAW_OUTDIR}/CircRNACount",
        linear_counts=f"{DCC_RAW_OUTDIR}/LinearCount",
        circ_coordinates=f"{DCC_RAW_OUTDIR}/CircCoordinates"
    output:
        results=f"{DCC_CIRCTEST_OUTDIR}/pairwise/{{comparison}}/circtest_results.tsv",
        summary=f"{DCC_CIRCTEST_OUTDIR}/pairwise/{{comparison}}/circtest_summary.tsv",
        filtered_circ=f"{DCC_CIRCTEST_OUTDIR}/pairwise/{{comparison}}/filtered_circ_counts.tsv",
        filtered_linear=f"{DCC_CIRCTEST_OUTDIR}/pairwise/{{comparison}}/filtered_linear_counts.tsv",
        ratio_plots=f"{DCC_CIRCTEST_OUTDIR}/pairwise/{{comparison}}/ratio_plots.pdf"
    wildcard_constraints:
        comparison=DEG_COMPARISON_REGEX
    params:
        all_samples=SAMPLES,
        selected_samples=lambda wc: (
            list(DEG_COMPARISONS[wc.comparison]["case"])
            + list(DEG_COMPARISONS[wc.comparison]["control"])
        ),
        group_labels=lambda wc: (
            [DEG_COMPARISONS[wc.comparison]["case_group"]]
            * len(DEG_COMPARISONS[wc.comparison]["case"])
            + [DEG_COMPARISONS[wc.comparison]["control_group"]]
            * len(DEG_COMPARISONS[wc.comparison]["control"])
        ),
        comparison=lambda wc: wc.comparison,
        case_group=lambda wc: DEG_COMPARISONS[wc.comparison]["case_group"],
        control_group=lambda wc: DEG_COMPARISONS[wc.comparison]["control_group"],
        padj_cutoff=DEG_PADJ_CUTOFF,
        filter_sample=DCC_CIRCTEST_FILTER_SAMPLE,
        filter_count=DCC_CIRCTEST_FILTER_COUNT,
        percentage=DCC_CIRCTEST_PERCENTAGE,
        max_plots=DCC_CIRCTEST_MAX_PLOTS
    conda:
        "envs/dcc.yaml"
    log:
        "logs/dcc_deg/circtest.{comparison}.log"
    script:
        "scripts/run_dcc_circtest.R"


rule ciri3_de_prepare_bsj:
    input:
        bsj=f"{OUTDIR}/ciri3/all_samples.ciri3.BSJ_Matrix",
        featurecounts=f"{OUTDIR}/featurecount/totalRNA.counts.txt",
        ciri3=lambda wc: ciri3_per_sample_result_paths(wc.comparison)
    output:
        info=f"{CIRI3_DE_OUTDIR}/{{comparison}}/de_bsj/infor.tsv",
        gene_expression=f"{CIRI3_DE_OUTDIR}/{{comparison}}/de_bsj/Gene_Expression.txt",
        bsj_matrix=f"{CIRI3_DE_OUTDIR}/{{comparison}}/de_bsj/BSJ_Matrix.txt"
    log:
        "logs/ciri3_de/{comparison}.prepare_bsj.log"
    params:
        selected_samples=lambda wc: ciri3_de_all_samples(wc.comparison),
        sample_classes=lambda wc: ciri3_de_sample_class_map(wc.comparison)
    conda:
        "envs/ciri3.yaml"
    script:
        "scripts/prepare_ciri3_de_bsj_inputs.py"


rule ciri3_de_prepare_ratio:
    input:
        bsj=f"{OUTDIR}/ciri3/all_samples.ciri3.BSJ_Matrix",
        fsj=f"{OUTDIR}/ciri3/all_samples.ciri3.FSJ_Matrix"
    output:
        info=f"{CIRI3_DE_OUTDIR}/{{comparison}}/de_ratio/infor.tsv",
        bsj_matrix=f"{CIRI3_DE_OUTDIR}/{{comparison}}/de_ratio/BSJ_Matrix.txt",
        fsj_matrix=f"{CIRI3_DE_OUTDIR}/{{comparison}}/de_ratio/FSJ_Matrix.txt"
    log:
        "logs/ciri3_de/{comparison}.prepare_ratio.log"
    params:
        selected_samples=lambda wc: ciri3_de_all_samples(wc.comparison),
        sample_classes=lambda wc: ciri3_de_sample_class_map(wc.comparison)
    conda:
        "envs/ciri3.yaml"
    script:
        "scripts/prepare_ciri3_de_ratio_inputs.py"


rule ciri3_de_prepare_relative:
    input:
        bsj=f"{OUTDIR}/ciri3/all_samples.ciri3.BSJ_Matrix",
        ciri3_all=f"{OUTDIR}/ciri3/all_samples.ciri3",
        ciri3=lambda wc: ciri3_per_sample_result_paths(wc.comparison)
    output:
        info=f"{CIRI3_DE_OUTDIR}/{{comparison}}/de_relative/infor.tsv",
        bsj_matrix=f"{CIRI3_DE_OUTDIR}/{{comparison}}/de_relative/BSJ_Matrix.txt",
        circ_gene=f"{CIRI3_DE_OUTDIR}/{{comparison}}/de_relative/circ_Gene.txt"
    log:
        "logs/ciri3_de/{comparison}.prepare_relative.log"
    params:
        selected_samples=lambda wc: ciri3_de_all_samples(wc.comparison),
        sample_classes=lambda wc: ciri3_de_sample_class_map(wc.comparison)
    conda:
        "envs/ciri3.yaml"
    script:
        "scripts/prepare_ciri3_de_relative_inputs.py"


rule ciri3_run_de_bsj:
    input:
        jar=CIRI3_JAR,
        ciri3_ready=CIRI3_READY,
        bsj_yes=CIRI3_BSJ_YES,
        info=f"{CIRI3_DE_OUTDIR}/{{comparison}}/de_bsj/infor.tsv",
        gene_expression=f"{CIRI3_DE_OUTDIR}/{{comparison}}/de_bsj/Gene_Expression.txt",
        bsj_matrix=f"{CIRI3_DE_OUTDIR}/{{comparison}}/de_bsj/BSJ_Matrix.txt"
    output:
        result=f"{CIRI3_DE_OUTDIR}/{{comparison}}/de_bsj/result.txt"
    conda:
        "envs/ciri3.yaml"
    log:
        "logs/ciri3_de/{comparison}.de_bsj.log"
    shell:
        r"""
        set -uo pipefail
        output_dir=$(dirname "{output.result}")
        log_dir=$(dirname "{log}")
        mkdir -p "$output_dir" "$log_dir"

        set +e
        java -jar "{input.jar}" DE_BSJ \
          -I "{input.info}" \
          -G "{input.gene_expression}" \
          -M "{input.bsj_matrix}" \
          -O "{output.result}" \
          > "{log}" 2>&1
        status=$?
        set -e

        if [ ! -s "{output.result}" ] && [ -s "{output.result}_Control_Case" ]; then
          cp "{output.result}_Control_Case" "{output.result}"
        fi

        if [ ! -s "{output.result}" ]; then
          echo "" >> "{log}"
          echo "DE_BSJ did not create a non-empty result. Re-running the official BSJ_yes.R command to capture the hidden R error:" >> "{log}"
          diag_out="{output.result}.diagnostic"
          rm -f "$diag_out"
          set +e
          Rscript "{input.bsj_yes}" \
            "{input.info}" \
            "{input.bsj_matrix}" \
            "{input.gene_expression}" \
            "$diag_out" \
            >> "{log}" 2>&1
          r_status=$?
          set -e
          rm -f "$diag_out"
          echo "Official BSJ_yes.R diagnostic exit code: $r_status" >> "{log}"
        fi

        test -s "{output.result}" || (echo "DE_BSJ failed and did not create output. Exit code: $status" >&2; exit 1)
        """


rule ciri3_run_de_ratio:
    input:
        jar=CIRI3_JAR,
        ciri3_ready=CIRI3_READY,
        info=f"{CIRI3_DE_OUTDIR}/{{comparison}}/de_ratio/infor.tsv",
        bsj_matrix=f"{CIRI3_DE_OUTDIR}/{{comparison}}/de_ratio/BSJ_Matrix.txt",
        fsj_matrix=f"{CIRI3_DE_OUTDIR}/{{comparison}}/de_ratio/FSJ_Matrix.txt"
    output:
        result=f"{CIRI3_DE_OUTDIR}/{{comparison}}/de_ratio/result.txt"
    params:
        finalize_script=CIRI3_RATIO_RELATIVE_DEG_SCRIPT
    conda:
        "envs/ciri3.yaml"
    log:
        "logs/ciri3_de/{comparison}.de_ratio.log"
    shell:
        r"""
        set -uo pipefail
        output_dir=$(dirname "{output.result}")
        log_dir=$(dirname "{log}")
        mkdir -p "$output_dir" "$log_dir"
        rm -f "{output.result}" "{output.result}_"*

        set +e
        java -jar "{input.jar}" DE_Ratio \
          -I "{input.info}" \
          -BM "{input.bsj_matrix}" \
          -FM "{input.fsj_matrix}" \
          -O "{output.result}" \
          > "{log}" 2>&1
        status=$?
        set -e

        if [ "$status" -ne 0 ]; then
          echo "" >> "{log}"
          echo "CIRI3 DE_Ratio exited with code $status; attempting native finalization from prepared matrices." >> "{log}"
        fi

        Rscript "{params.finalize_script}" \
          --method de_ratio \
          --info "{input.info}" \
          --bsj "{input.bsj_matrix}" \
          --fsj "{input.fsj_matrix}" \
          --out "{output.result}" \
          >> "{log}" 2>&1

        test -s "{output.result}" || (echo "DE_Ratio finalization failed and did not create output. CIRI3 exit code: $status" >&2; exit 1)
        """


rule ciri3_run_de_relative:
    input:
        jar=CIRI3_JAR,
        ciri3_ready=CIRI3_READY,
        info=f"{CIRI3_DE_OUTDIR}/{{comparison}}/de_relative/infor.tsv",
        bsj_matrix=f"{CIRI3_DE_OUTDIR}/{{comparison}}/de_relative/BSJ_Matrix.txt",
        circ_gene=f"{CIRI3_DE_OUTDIR}/{{comparison}}/de_relative/circ_Gene.txt"
    output:
        result=f"{CIRI3_DE_OUTDIR}/{{comparison}}/de_relative/result.txt"
    params:
        finalize_script=CIRI3_RATIO_RELATIVE_DEG_SCRIPT
    conda:
        "envs/ciri3.yaml"
    log:
        "logs/ciri3_de/{comparison}.de_relative.log"
    shell:
        r"""
        set -uo pipefail
        output_dir=$(dirname "{output.result}")
        log_dir=$(dirname "{log}")
        mkdir -p "$output_dir" "$log_dir"
        rm -f "{output.result}" "{output.result}_"*

        set +e
        java -jar "{input.jar}" DE_Relative \
          -I "{input.info}" \
          -M "{input.bsj_matrix}" \
          -GC "{input.circ_gene}" \
          -O "{output.result}" \
          > "{log}" 2>&1
        status=$?
        set -e

        if [ "$status" -ne 0 ]; then
          echo "" >> "{log}"
          echo "CIRI3 DE_Relative exited with code $status; attempting native finalization from prepared matrices." >> "{log}"
        fi

        Rscript "{params.finalize_script}" \
          --method de_relative \
          --info "{input.info}" \
          --bsj "{input.bsj_matrix}" \
          --circ-gene "{input.circ_gene}" \
          --out "{output.result}" \
          >> "{log}" 2>&1

        test -s "{output.result}" || (echo "DE_Relative finalization failed and did not create output. CIRI3 exit code: $status" >&2; exit 1)
        """


rule ciri3_de_bsj_pairwise_plots:
    input:
        result=f"{CIRI3_DE_OUTDIR}/{{comparison}}/de_bsj/result.txt",
        bsj_matrix=f"{CIRI3_DE_OUTDIR}/{{comparison}}/de_bsj/BSJ_Matrix.txt",
        ciri3=f"{OUTDIR}/ciri3/all_samples.ciri3"
    output:
        heatmap=f"{CIRI3_DE_OUTDIR}/{{comparison}}/de_bsj/heatmap_top50.pdf",
        pca=f"{CIRI3_DE_OUTDIR}/{{comparison}}/de_bsj/pca.pdf",
        volcano=f"{CIRI3_DE_OUTDIR}/{{comparison}}/de_bsj/volcano.pdf",
        volcano_labeled=f"{CIRI3_DE_OUTDIR}/{{comparison}}/de_bsj/volcano_labeled.pdf"
    wildcard_constraints:
        comparison=CIRI3_DE_COMPARISON_REGEX
    params:
        method="de_bsj",
        mode="pairwise",
        comparison=lambda wc: wc.comparison,
        groups=DEG_GROUPS,
        comparisons=DEG_COMPARISONS,
        padj_cutoff=DEG_PADJ_CUTOFF,
        lfc_cutoff=DEG_LFC_CUTOFF
    conda:
        "envs/deg.yaml"
    log:
        "logs/ciri3_de/{comparison}.de_bsj.plots.log"
    script:
        "scripts/plot_ciri3_deg.R"


rule ciri3_de_ratio_pairwise_plots:
    input:
        result=f"{CIRI3_DE_OUTDIR}/{{comparison}}/de_ratio/result.txt",
        bsj_matrix=f"{CIRI3_DE_OUTDIR}/{{comparison}}/de_ratio/BSJ_Matrix.txt",
        fsj_matrix=f"{CIRI3_DE_OUTDIR}/{{comparison}}/de_ratio/FSJ_Matrix.txt",
        ciri3=f"{OUTDIR}/ciri3/all_samples.ciri3"
    output:
        heatmap=f"{CIRI3_DE_OUTDIR}/{{comparison}}/de_ratio/heatmap_top50.pdf",
        pca=f"{CIRI3_DE_OUTDIR}/{{comparison}}/de_ratio/pca.pdf",
        volcano=f"{CIRI3_DE_OUTDIR}/{{comparison}}/de_ratio/volcano.pdf",
        volcano_labeled=f"{CIRI3_DE_OUTDIR}/{{comparison}}/de_ratio/volcano_labeled.pdf"
    wildcard_constraints:
        comparison=CIRI3_DE_COMPARISON_REGEX
    params:
        method="de_ratio",
        mode="pairwise",
        comparison=lambda wc: wc.comparison,
        groups=DEG_GROUPS,
        comparisons=DEG_COMPARISONS,
        padj_cutoff=DEG_PADJ_CUTOFF,
        lfc_cutoff=DEG_LFC_CUTOFF
    conda:
        "envs/deg.yaml"
    log:
        "logs/ciri3_de/{comparison}.de_ratio.plots.log"
    script:
        "scripts/plot_ciri3_deg.R"


rule ciri3_de_relative_pairwise_plots:
    input:
        result=f"{CIRI3_DE_OUTDIR}/{{comparison}}/de_relative/result.txt",
        bsj_matrix=f"{CIRI3_DE_OUTDIR}/{{comparison}}/de_relative/BSJ_Matrix.txt",
        ciri3=f"{OUTDIR}/ciri3/all_samples.ciri3"
    output:
        heatmap=f"{CIRI3_DE_OUTDIR}/{{comparison}}/de_relative/heatmap_top50.pdf",
        pca=f"{CIRI3_DE_OUTDIR}/{{comparison}}/de_relative/pca.pdf",
        volcano=f"{CIRI3_DE_OUTDIR}/{{comparison}}/de_relative/volcano.pdf",
        volcano_labeled=f"{CIRI3_DE_OUTDIR}/{{comparison}}/de_relative/volcano_labeled.pdf"
    wildcard_constraints:
        comparison=CIRI3_DE_COMPARISON_REGEX
    params:
        method="de_relative",
        mode="pairwise",
        comparison=lambda wc: wc.comparison,
        groups=DEG_GROUPS,
        comparisons=DEG_COMPARISONS,
        padj_cutoff=DEG_PADJ_CUTOFF,
        lfc_cutoff=DEG_LFC_CUTOFF
    conda:
        "envs/deg.yaml"
    log:
        "logs/ciri3_de/{comparison}.de_relative.plots.log"
    script:
        "scripts/plot_ciri3_deg.R"


rule ciri3_de_bsj_all_groups_plots:
    input:
        results=expand(
            f"{CIRI3_DE_OUTDIR}/{{comparison}}/de_bsj/result.txt",
            comparison=CIRI3_DE_COMPARISON_NAMES,
        ),
        bsj_matrix=f"{OUTDIR}/ciri3/all_samples.ciri3.BSJ_Matrix",
        ciri3=f"{OUTDIR}/ciri3/all_samples.ciri3"
    output:
        heatmap=f"{CIRI3_DE_OUTDIR}/all_groups/de_bsj/heatmap_top50.pdf",
        pca=f"{CIRI3_DE_OUTDIR}/all_groups/de_bsj/pca.pdf"
    params:
        method="de_bsj",
        mode="all_groups",
        groups=DEG_GROUPS,
        comparisons=DEG_COMPARISONS,
        padj_cutoff=DEG_PADJ_CUTOFF,
        lfc_cutoff=DEG_LFC_CUTOFF
    conda:
        "envs/deg.yaml"
    log:
        "logs/ciri3_de/all_groups.de_bsj.plots.log"
    script:
        "scripts/plot_ciri3_deg.R"


rule ciri3_de_ratio_all_groups_plots:
    input:
        results=expand(
            f"{CIRI3_DE_OUTDIR}/{{comparison}}/de_ratio/result.txt",
            comparison=CIRI3_DE_COMPARISON_NAMES,
        ),
        bsj_matrix=f"{OUTDIR}/ciri3/all_samples.ciri3.BSJ_Matrix",
        fsj_matrix=f"{OUTDIR}/ciri3/all_samples.ciri3.FSJ_Matrix",
        ciri3=f"{OUTDIR}/ciri3/all_samples.ciri3"
    output:
        heatmap=f"{CIRI3_DE_OUTDIR}/all_groups/de_ratio/heatmap_top50.pdf",
        pca=f"{CIRI3_DE_OUTDIR}/all_groups/de_ratio/pca.pdf"
    params:
        method="de_ratio",
        mode="all_groups",
        groups=DEG_GROUPS,
        comparisons=DEG_COMPARISONS,
        padj_cutoff=DEG_PADJ_CUTOFF,
        lfc_cutoff=DEG_LFC_CUTOFF
    conda:
        "envs/deg.yaml"
    log:
        "logs/ciri3_de/all_groups.de_ratio.plots.log"
    script:
        "scripts/plot_ciri3_deg.R"


rule ciri3_de_relative_all_groups_plots:
    input:
        results=expand(
            f"{CIRI3_DE_OUTDIR}/{{comparison}}/de_relative/result.txt",
            comparison=CIRI3_DE_COMPARISON_NAMES,
        ),
        bsj_matrix=f"{OUTDIR}/ciri3/all_samples.ciri3.BSJ_Matrix",
        ciri3=f"{OUTDIR}/ciri3/all_samples.ciri3"
    output:
        heatmap=f"{CIRI3_DE_OUTDIR}/all_groups/de_relative/heatmap_top50.pdf",
        pca=f"{CIRI3_DE_OUTDIR}/all_groups/de_relative/pca.pdf"
    params:
        method="de_relative",
        mode="all_groups",
        groups=DEG_GROUPS,
        comparisons=DEG_COMPARISONS,
        padj_cutoff=DEG_PADJ_CUTOFF,
        lfc_cutoff=DEG_LFC_CUTOFF
    conda:
        "envs/deg.yaml"
    log:
        "logs/ciri3_de/all_groups.de_relative.plots.log"
    script:
        "scripts/plot_ciri3_deg.R"
