OUTDIR = config["output"]["dir"]
SAMPLES = list(config["samples"].keys())

PAIRWISE_NAMES = [f"{a}_vs_{b}" for a, b in DEG_PAIRWISE]
PAIRWISE_TO_GROUPS = {
    f"{a}_vs_{b}": (a, b)
    for a, b in DEG_PAIRWISE
}


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
        pairwise_volcano_labeled=expand(
            f"{OUTDIR}/deg/bsj/pairwise/{{pair}}/volcano_labeled.pdf",
            pair=[f"{a}_vs_{b}" for a, b in DEG_PAIRWISE],
        ),
        pairwise_heatmap=expand(
            f"{OUTDIR}/deg/bsj/pairwise/{{pair}}/heatmap_top50.pdf",
            pair=[f"{a}_vs_{b}" for a, b in DEG_PAIRWISE],
        ),
        pairwise_pca=expand(
            f"{OUTDIR}/deg/bsj/pairwise/{{pair}}/pca.pdf",
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


rule prepare_ciri3_de_bsj_inputs:
    input:
        bsj_matrix=f"{OUTDIR}/ciri3/all_samples.ciri3.BSJ_Matrix",
        gene_counts=f"{OUTDIR}/featurecount/totalRNA.counts.txt",
        ciri3=expand(f"{OUTDIR}/ciri3/per_sample/{{sample}}.ciri3", sample=SAMPLES)
    output:
        infor=f"{OUTDIR}/deg/ciri3/de_bsj/pairwise/{{pair}}/infor.tsv",
        bsj_matrix=f"{OUTDIR}/deg/ciri3/de_bsj/pairwise/{{pair}}/BSJ_Matrix.txt",
        gene_expression=f"{OUTDIR}/deg/ciri3/de_bsj/pairwise/{{pair}}/Gene_Expression.txt"
    params:
        pair=lambda wc: wc.pair,
        group_a=lambda wc: PAIRWISE_TO_GROUPS[wc.pair][0],
        group_b=lambda wc: PAIRWISE_TO_GROUPS[wc.pair][1],
        groups=DEG_GROUPS,
        sample_names=SAMPLES
    script:
        "scripts/prepare_ciri3_de_bsj_inputs.py"


rule run_ciri3_de_bsj_pairwise:
    input:
        infor=f"{OUTDIR}/deg/ciri3/de_bsj/pairwise/{{pair}}/infor.tsv",
        bsj_matrix=f"{OUTDIR}/deg/ciri3/de_bsj/pairwise/{{pair}}/BSJ_Matrix.txt",
        gene_expression=f"{OUTDIR}/deg/ciri3/de_bsj/pairwise/{{pair}}/Gene_Expression.txt"
    output:
        result=f"{OUTDIR}/deg/ciri3/de_bsj/pairwise/{{pair}}/result.txt"
    params:
        jar=config["ciri3"]["jar"]
    conda:
        "envs/ciri3.yaml"
    log:
        "logs/deg/ciri3/de_bsj/{pair}.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.result}) $(dirname {log})
        java -jar {params.jar} DE_BSJ \
          -I {input.infor} \
          -G {input.gene_expression} \
          -M {input.bsj_matrix} \
          -O {output.result} \
          > {log} 2>&1
        test -s {output.result}
        """


rule merge_ciri3_de_bsj_all_samples:
    input:
        results=expand(
            f"{OUTDIR}/deg/ciri3/de_bsj/pairwise/{{pair}}/result.txt",
            pair=PAIRWISE_NAMES,
        )
    output:
        all_samples=f"{OUTDIR}/deg/ciri3/de_bsj/all_samples/result.txt"
    params:
        pairs=PAIRWISE_NAMES
    script:
        "scripts/merge_ciri3_pairwise_deg.py"


rule prepare_ciri3_de_ratio_inputs:
    input:
        bsj_matrix=f"{OUTDIR}/ciri3/all_samples.ciri3.BSJ_Matrix",
        fsj_matrix=f"{OUTDIR}/ciri3/all_samples.ciri3.FSJ_Matrix"
    output:
        infor=f"{OUTDIR}/deg/ciri3/de_ratio/pairwise/{{pair}}/infor.tsv",
        bsj_matrix=f"{OUTDIR}/deg/ciri3/de_ratio/pairwise/{{pair}}/BSJ_Matrix.txt",
        fsj_matrix=f"{OUTDIR}/deg/ciri3/de_ratio/pairwise/{{pair}}/FSJ_Matrix.txt"
    params:
        pair=lambda wc: wc.pair,
        group_a=lambda wc: PAIRWISE_TO_GROUPS[wc.pair][0],
        group_b=lambda wc: PAIRWISE_TO_GROUPS[wc.pair][1],
        groups=DEG_GROUPS
    script:
        "scripts/prepare_ciri3_de_ratio_inputs.py"


rule run_ciri3_de_ratio_pairwise:
    input:
        infor=f"{OUTDIR}/deg/ciri3/de_ratio/pairwise/{{pair}}/infor.tsv",
        bsj_matrix=f"{OUTDIR}/deg/ciri3/de_ratio/pairwise/{{pair}}/BSJ_Matrix.txt",
        fsj_matrix=f"{OUTDIR}/deg/ciri3/de_ratio/pairwise/{{pair}}/FSJ_Matrix.txt"
    output:
        result=f"{OUTDIR}/deg/ciri3/de_ratio/pairwise/{{pair}}/result.txt"
    params:
        jar=config["ciri3"]["jar"]
    threads: int(config["threads"].get("deg_ratio", 1))
    conda:
        "envs/ciri3.yaml"
    log:
        "logs/deg/ciri3/de_ratio/{pair}.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.result}) $(dirname {log})
        java -jar {params.jar} DE_Ratio \
          -I {input.infor} \
          -BM {input.bsj_matrix} \
          -FM {input.fsj_matrix} \
          -T {threads} \
          -O {output.result} \
          > {log} 2>&1
        test -s {output.result}
        """


rule merge_ciri3_de_ratio_all_samples:
    input:
        results=expand(
            f"{OUTDIR}/deg/ciri3/de_ratio/pairwise/{{pair}}/result.txt",
            pair=PAIRWISE_NAMES,
        )
    output:
        all_samples=f"{OUTDIR}/deg/ciri3/de_ratio/all_samples/result.txt"
    params:
        pairs=PAIRWISE_NAMES
    script:
        "scripts/merge_ciri3_pairwise_deg.py"


rule prepare_ciri3_de_relative_inputs:
    input:
        ciri3=expand(f"{OUTDIR}/ciri3/per_sample/{{sample}}.ciri3", sample=SAMPLES)
    output:
        infor=f"{OUTDIR}/deg/ciri3/de_relative/pairwise/{{pair}}/infor.tsv"
    params:
        pair=lambda wc: wc.pair,
        group_a=lambda wc: PAIRWISE_TO_GROUPS[wc.pair][0],
        group_b=lambda wc: PAIRWISE_TO_GROUPS[wc.pair][1],
        groups=DEG_GROUPS,
        sample_names=SAMPLES
    script:
        "scripts/prepare_ciri3_de_relative_inputs.py"


rule run_ciri3_de_relative_pairwise:
    input:
        infor=f"{OUTDIR}/deg/ciri3/de_relative/pairwise/{{pair}}/infor.tsv"
    output:
        result=f"{OUTDIR}/deg/ciri3/de_relative/pairwise/{{pair}}/result.txt"
    params:
        jar=config["ciri3"]["jar"]
    threads: int(config["threads"].get("deg_relative", 1))
    conda:
        "envs/ciri3.yaml"
    log:
        "logs/deg/ciri3/de_relative/{pair}.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.result}) $(dirname {log})
        java -jar {params.jar} DE_Relative \
          -I {input.infor} \
          -T {threads} \
          -O {output.result} \
          > {log} 2>&1
        test -s {output.result}
        """


rule merge_ciri3_de_relative_all_samples:
    input:
        results=expand(
            f"{OUTDIR}/deg/ciri3/de_relative/pairwise/{{pair}}/result.txt",
            pair=PAIRWISE_NAMES,
        )
    output:
        all_samples=f"{OUTDIR}/deg/ciri3/de_relative/all_samples/result.txt"
    params:
        pairs=PAIRWISE_NAMES
    script:
        "scripts/merge_ciri3_pairwise_deg.py"
