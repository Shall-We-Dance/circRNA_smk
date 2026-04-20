import csv
import os
import re

OUTDIR = config["output"]["dir"]
SAMPLES = list(config["samples"].keys())

PAIRWISE_NAMES = [f"{a}_vs_{b}" for a, b in DEG_PAIRWISE]
PAIRWISE_TO_GROUPS = {
    f"{a}_vs_{b}": (a, b)
    for a, b in DEG_PAIRWISE
}


def ciri3_de_case_samples(comparison):
    return list(CIRI3_DE_COMPARISONS[comparison]["case"])


def ciri3_de_control_samples(comparison):
    return list(CIRI3_DE_COMPARISONS[comparison]["control"])


def ciri3_de_all_samples(comparison):
    return ciri3_de_case_samples(comparison) + ciri3_de_control_samples(comparison)


def ciri3_de_sample_class_map(comparison):
    out = {}
    for sample in ciri3_de_case_samples(comparison):
        out[sample] = "Case"
    for sample in ciri3_de_control_samples(comparison):
        out[sample] = "Control"
    return out


def ciri3_per_sample_result_paths(comparison):
    return [
        f"{OUTDIR}/ciri3/per_sample/{sample}.ciri3"
        for sample in ciri3_de_all_samples(comparison)
    ]


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


rule ciri3_de_prepare_bsj:
    input:
        bsj=f"{OUTDIR}/ciri3/all_samples.ciri3.BSJ_Matrix",
        featurecounts=f"{OUTDIR}/featurecount/totalRNA.counts.txt",
        ciri3=lambda wc: ciri3_per_sample_result_paths(wc.comparison)
    output:
        info=f"{CIRI3_DE_OUTDIR}/{{comparison}}/de_bsj/infor.tsv",
        gene_expression=f"{CIRI3_DE_OUTDIR}/{{comparison}}/de_bsj/Gene_Expression.txt",
        bsj_matrix=f"{CIRI3_DE_OUTDIR}/{{comparison}}/de_bsj/BSJ_Matrix.txt"
    run:
        comparison = wildcards.comparison
        selected_samples = ciri3_de_all_samples(comparison)
        sample_classes = ciri3_de_sample_class_map(comparison)

        os.makedirs(os.path.dirname(output.info), exist_ok=True)
        os.makedirs(os.path.dirname(output.gene_expression), exist_ok=True)
        os.makedirs(os.path.dirname(output.bsj_matrix), exist_ok=True)

        with open(output.info, "w", newline="") as fh:
            writer = csv.writer(fh, delimiter="\t")
            writer.writerow(["Sample", "Path", "Class"])
            for sample in selected_samples:
                writer.writerow([
                    sample,
                    os.path.abspath(f"{OUTDIR}/ciri3/per_sample/{sample}.ciri3"),
                    sample_classes[sample],
                ])

        with open(input.bsj) as fh:
            reader = csv.reader(fh, delimiter="\t")
            rows = list(reader)
        if not rows:
            raise ValueError(f"Empty BSJ matrix: {input.bsj}")

        header = rows[0]
        missing = [s for s in selected_samples if s not in header]
        if missing:
            raise ValueError(
                "Selected samples missing from merged BSJ matrix for comparison "
                f"{comparison}: " + ", ".join(missing)
            )

        keep_idx = [0] + [header.index(s) for s in selected_samples]
        with open(output.bsj_matrix, "w", newline="") as fh:
            writer = csv.writer(fh, delimiter="\t")
            writer.writerow([header[i] for i in keep_idx])
            for row in rows[1:]:
                padded = row + [""] * (len(header) - len(row))
                writer.writerow([padded[i] if padded[i] != "" else "0" for i in keep_idx])

        with open(input.featurecounts) as fh:
            raw_lines = [line.rstrip("\n") for line in fh if line.strip()]

        fc_lines = [line for line in raw_lines if not line.startswith("#")]
        if len(fc_lines) < 2:
            raise ValueError(f"featureCounts file appears malformed: {input.featurecounts}")

        fc_rows = list(csv.reader(fc_lines, delimiter="\t"))
        fc_header = fc_rows[0]

        if "Geneid" not in fc_header:
            raise ValueError("featureCounts header does not contain 'Geneid'.")

        geneid_idx = fc_header.index("Geneid")

        count_start_idx = 6
        if len(fc_header[count_start_idx:]) != len(SAMPLES):
            raise ValueError(
                "Unexpected number of featureCounts sample columns. "
                f"Found {len(fc_header[count_start_idx:])}, expected {len(SAMPLES)}."
            )

        sample_to_fc_idx = {}
        for i, sample in enumerate(SAMPLES):
            sample_to_fc_idx[sample] = count_start_idx + i

        missing_fc = [s for s in selected_samples if s not in sample_to_fc_idx]
        if missing_fc:
            raise ValueError(
                "Selected samples missing from featureCounts sample mapping for comparison "
                f"{comparison}: " + ", ".join(missing_fc)
            )

        with open(output.gene_expression, "w", newline="") as fh:
            writer = csv.writer(fh, delimiter="\t")
            writer.writerow(["Geneid"] + selected_samples)
            for row in fc_rows[1:]:
                padded = row + ["0"] * (len(fc_header) - len(row))
                out_row = [padded[geneid_idx]]
                for sample in selected_samples:
                    value = padded[sample_to_fc_idx[sample]]
                    out_row.append(value if value != "" else "0")
                writer.writerow(out_row)


rule ciri3_de_prepare_ratio:
    input:
        bsj=f"{OUTDIR}/ciri3/all_samples.ciri3.BSJ_Matrix",
        fsj=f"{OUTDIR}/ciri3/all_samples.ciri3.FSJ_Matrix"
    output:
        info=f"{CIRI3_DE_OUTDIR}/{{comparison}}/de_ratio/infor.tsv",
        bsj_matrix=f"{CIRI3_DE_OUTDIR}/{{comparison}}/de_ratio/BSJ_Matrix.txt",
        fsj_matrix=f"{CIRI3_DE_OUTDIR}/{{comparison}}/de_ratio/FSJ_Matrix.txt"
    run:
        comparison = wildcards.comparison
        selected_samples = ciri3_de_all_samples(comparison)
        sample_classes = ciri3_de_sample_class_map(comparison)

        os.makedirs(os.path.dirname(output.info), exist_ok=True)
        os.makedirs(os.path.dirname(output.bsj_matrix), exist_ok=True)
        os.makedirs(os.path.dirname(output.fsj_matrix), exist_ok=True)

        with open(output.info, "w", newline="") as fh:
            writer = csv.writer(fh, delimiter="\t")
            writer.writerow(["Sample", "Class"])
            for sample in selected_samples:
                writer.writerow([sample, sample_classes[sample]])

        for in_path, out_path in [(input.bsj, output.bsj_matrix), (input.fsj, output.fsj_matrix)]:
            with open(in_path) as fh:
                reader = csv.reader(fh, delimiter="\t")
                rows = list(reader)
            if not rows:
                raise ValueError(f"Empty matrix: {in_path}")

            header = rows[0]
            missing = [s for s in selected_samples if s not in header]
            if missing:
                raise ValueError(
                    f"Selected samples missing from matrix {in_path} for comparison "
                    f"{comparison}: " + ", ".join(missing)
                )

            keep_idx = [0] + [header.index(s) for s in selected_samples]
            with open(out_path, "w", newline="") as fh:
                writer = csv.writer(fh, delimiter="\t")
                writer.writerow([header[i] for i in keep_idx])
                for row in rows[1:]:
                    padded = row + [""] * (len(header) - len(row))
                    writer.writerow([padded[i] if padded[i] != "" else "0" for i in keep_idx])


rule ciri3_de_prepare_relative:
    input:
        bsj=f"{OUTDIR}/ciri3/all_samples.ciri3.BSJ_Matrix",
        ciri3_all=f"{OUTDIR}/ciri3/all_samples.ciri3",
        ciri3=lambda wc: ciri3_per_sample_result_paths(wc.comparison)
    output:
        info=f"{CIRI3_DE_OUTDIR}/{{comparison}}/de_relative/infor.tsv",
        bsj_matrix=f"{CIRI3_DE_OUTDIR}/{{comparison}}/de_relative/BSJ_Matrix.txt",
        circ_gene=f"{CIRI3_DE_OUTDIR}/{{comparison}}/de_relative/circ_Gene.txt"
    run:
        comparison = wildcards.comparison
        selected_samples = ciri3_de_all_samples(comparison)
        sample_classes = ciri3_de_sample_class_map(comparison)

        os.makedirs(os.path.dirname(output.info), exist_ok=True)
        os.makedirs(os.path.dirname(output.bsj_matrix), exist_ok=True)
        os.makedirs(os.path.dirname(output.circ_gene), exist_ok=True)

        with open(output.info, "w", newline="") as fh:
            writer = csv.writer(fh, delimiter="\t")
            writer.writerow(["Sample", "Path", "Class"])
            for sample in selected_samples:
                writer.writerow([
                    sample,
                    os.path.abspath(f"{OUTDIR}/ciri3/per_sample/{sample}.ciri3"),
                    sample_classes[sample],
                ])

        with open(input.bsj) as fh:
            reader = csv.reader(fh, delimiter="\t")
            rows = list(reader)
        if not rows:
            raise ValueError(f"Empty BSJ matrix: {input.bsj}")

        header = rows[0]
        missing = [s for s in selected_samples if s not in header]
        if missing:
            raise ValueError(
                "Selected samples missing from merged BSJ matrix for comparison "
                f"{comparison}: " + ", ".join(missing)
            )

        keep_idx = [0] + [header.index(s) for s in selected_samples]
        with open(output.bsj_matrix, "w", newline="") as fh:
            writer = csv.writer(fh, delimiter="\t")
            writer.writerow([header[i] for i in keep_idx])
            for row in rows[1:]:
                padded = row + [""] * (len(header) - len(row))
                writer.writerow([padded[i] if padded[i] != "" else "0" for i in keep_idx])

        with open(input.ciri3_all) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            if reader.fieldnames is None:
                raise ValueError(f"Empty or malformed CIRI3 table: {input.ciri3_all}")

            if "circRNA_ID" not in reader.fieldnames or "gene_id" not in reader.fieldnames:
                raise ValueError(
                    "all_samples.ciri3 must contain 'circRNA_ID' and 'gene_id' columns."
                )

            pairs = set()
            for row in reader:
                circ_id = row["circRNA_ID"].strip()
                gene_field = row["gene_id"].strip()

                if not circ_id or not gene_field or gene_field == "NA":
                    continue

                gene_parts = [
                    g.strip()
                    for g in re.split(r"[;,]", gene_field)
                    if g.strip() and g.strip() != "NA"
                ]
                for gene_id in gene_parts:
                    pairs.add((circ_id, gene_id))

        with open(output.circ_gene, "w", newline="") as fh:
            writer = csv.writer(fh, delimiter="\t")
            writer.writerow(["circRNA_ID", "gene_id"])
            for circ_id, gene_id in sorted(pairs):
                writer.writerow([circ_id, gene_id])


rule ciri3_run_de_bsj:
    input:
        info=f"{CIRI3_DE_OUTDIR}/{{comparison}}/de_bsj/infor.tsv",
        gene_expression=f"{CIRI3_DE_OUTDIR}/{{comparison}}/de_bsj/Gene_Expression.txt",
        bsj_matrix=f"{CIRI3_DE_OUTDIR}/{{comparison}}/de_bsj/BSJ_Matrix.txt"
    output:
        result=f"{CIRI3_DE_OUTDIR}/{{comparison}}/de_bsj/result.txt"
    params:
        ld_path="/projects/cell-surface/CIRI3/lib/rmats_deps",
        bsj_script="/projects/cell-surface/CIRI3/scripts/BSJ_yes.R"
    conda:
        "envs/ciri3.yaml"
    log:
        "logs/ciri3_de/{comparison}.de_bsj.log"
    shell:
        r"""
        set -uo pipefail
        export LD_LIBRARY_PATH={params.ld_path}:${{LD_LIBRARY_PATH:-}}
        mkdir -p $(dirname {output.result}) $(dirname {log})

        set +e
        Rscript --vanilla {params.bsj_script} \
          {input.info} \
          {input.bsj_matrix} \
          {input.gene_expression} \
          {output.result} \
          > {log} 2>&1
        status=$?
        set -e

        test -s {output.result} || (echo "DE_BSJ failed and did not create output. Exit code: $status" >&2; exit 1)
        """


rule ciri3_run_de_ratio:
    input:
        info=f"{CIRI3_DE_OUTDIR}/{{comparison}}/de_ratio/infor.tsv",
        bsj_matrix=f"{CIRI3_DE_OUTDIR}/{{comparison}}/de_ratio/BSJ_Matrix.txt",
        fsj_matrix=f"{CIRI3_DE_OUTDIR}/{{comparison}}/de_ratio/FSJ_Matrix.txt"
    output:
        result=f"{CIRI3_DE_OUTDIR}/{{comparison}}/de_ratio/result.txt"
    params:
        jar=config["ciri3"]["jar"],
        ld_path="/projects/cell-surface/CIRI3/lib/rmats_deps"
    conda:
        "envs/ciri3.yaml"
    log:
        "logs/ciri3_de/{comparison}.de_ratio.log"
    shell:
        r"""
        set -uo pipefail
        export LD_LIBRARY_PATH={params.ld_path}:${{LD_LIBRARY_PATH:-}}
        mkdir -p $(dirname {output.result}) $(dirname {log})

        set +e
        java -jar {params.jar} DE_Ratio \
          -I {input.info} \
          -BM {input.bsj_matrix} \
          -FM {input.fsj_matrix} \
          -O {output.result} \
          > {log} 2>&1
        status=$?
        set -e

        if [ ! -s {output.result} ] && [ -s {output.result}_Control_Case ]; then
          cp {output.result}_Control_Case {output.result}
        fi

        test -s {output.result} || (echo "DE_Ratio failed and did not create output. Exit code: $status" >&2; exit 1)
        """


rule ciri3_run_de_relative:
    input:
        info=f"{CIRI3_DE_OUTDIR}/{{comparison}}/de_relative/infor.tsv",
        bsj_matrix=f"{CIRI3_DE_OUTDIR}/{{comparison}}/de_relative/BSJ_Matrix.txt",
        circ_gene=f"{CIRI3_DE_OUTDIR}/{{comparison}}/de_relative/circ_Gene.txt"
    output:
        result=f"{CIRI3_DE_OUTDIR}/{{comparison}}/de_relative/result.txt"
    params:
        jar=config["ciri3"]["jar"],
        ld_path="/projects/cell-surface/CIRI3/lib/rmats_deps"
    conda:
        "envs/ciri3.yaml"
    log:
        "logs/ciri3_de/{comparison}.de_relative.log"
    shell:
        r"""
        set -uo pipefail
        export LD_LIBRARY_PATH={params.ld_path}:${{LD_LIBRARY_PATH:-}}
        mkdir -p $(dirname {output.result}) $(dirname {log})

        set +e
        java -jar {params.jar} DE_Relative \
          -I {input.info} \
          -M {input.bsj_matrix} \
          -GC {input.circ_gene} \
          -O {output.result} \
          > {log} 2>&1
        status=$?
        set -e

        if [ ! -s {output.result} ] && [ -s {output.result}_Control_Case ]; then
          cp {output.result}_Control_Case {output.result}
        fi

        test -s {output.result} || (echo "DE_Relative failed and did not create output. Exit code: $status" >&2; exit 1)
        """
