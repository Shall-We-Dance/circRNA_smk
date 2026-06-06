# workflow/rules/circtools.smk
#
# circtools integration (https://github.com/jakobilab/circtools).
#
# STAR alignment matches the OFFICIAL DCC tutorial command exactly:
# https://github.com/dieterich-lab/DCC tutorial section "Paired-end mapping".
# This is the authoritative reference; we follow it parameter-for-parameter
# so circtools/DCC detect results are directly comparable to other groups
# using the same upstream protocol.
#
# Three STAR alignments per sample:
#   1. circtools_star_paired  - paired-end joint mapping (-> sorted BAM for -B)
#   2. circtools_star_mate1   - mate1 (R1) single-end mapping (no BAM)
#   3. circtools_star_mate2   - mate2 (R2) single-end mapping (no BAM)
#
# Sample order across paired/mate1/mate2/BAM sheets MUST be identical -
# DCC matches them positionally.


# ---------------------------------------------------------------------------
# Install marker.
# ---------------------------------------------------------------------------
rule circtools_install_marker:
    output:
        marker=CIRCTOOLS_INSTALL_MARKER
    log:
        "logs/circtools/install.log"
    conda:
        "envs/circtools.yaml"
    params:
        need_primex_r=int(CIRCTOOLS_PRIMEX_ENABLED)
    shell:
        r"""
        set -euo pipefail
        log_dir=$(dirname "{log}")
        marker_dir=$(dirname "{output.marker}")
        mkdir -p "$log_dir" "$marker_dir"

        circtools --version > "{log}" 2>&1 || {{
          echo "circtools --version failed; conda env may not have built correctly." >> "{log}"
          exit 1
        }}

        if [ "{params.need_primex_r}" = "1" ]; then
          primex_src=$(python -c 'import circtools, os; print(os.path.join(os.path.dirname(circtools.__file__), "contrib", "primex"))') 2>>"{log}"
          if [ ! -d "$primex_src" ]; then
            echo "Bundled primex R package not found at: $primex_src" >> "{log}"
            exit 1
          fi
          R CMD INSTALL "$primex_src" >> "{log}" 2>&1 || {{
            echo "R CMD INSTALL primex failed." >> "{log}"
            exit 1
          }}
          Rscript -e 'suppressMessages(library(primex)); cat("primex R package loads OK\n")' >> "{log}" 2>&1 || exit 1
        else
          echo "primex not enabled; skipping R primex install." >> "{log}"
        fi

        date -u +"circtools install verified at %Y-%m-%dT%H:%M:%SZ" > "{output.marker}"
        """


# ---------------------------------------------------------------------------
# STAR 1/3: Paired-end joint mapping.
# Matches the official DCC paired-end command exactly, with two minimal additions:
#   --outSAMtype BAM SortedByCoordinate  (instead of unspecified - we need
#     the sorted BAM for DCC's -B host gene counting flag)
#   (post-step) samtools index, to satisfy DCC -B
# Everything else is verbatim from the DCC tutorial.
# ---------------------------------------------------------------------------
rule circtools_star_paired:
    input:
        idx_sa=f"{config['reference']['star_index']}/SA",
        r1=f"{OUTDIR}/tmp/fastp/{{sample}}_R1.fastq.gz",
        r2=f"{OUTDIR}/tmp/fastp/{{sample}}_R2.fastq.gz"
    output:
        chimeric_junction=f"{OUTDIR}/circtools/star/{{sample}}/paired/{{sample}}_Chimeric.out.junction",
        log_final=f"{OUTDIR}/circtools/star/{{sample}}/paired/{{sample}}_Log.final.out",
        sj=f"{OUTDIR}/circtools/star/{{sample}}/paired/{{sample}}_SJ.out.tab",
        bam=f"{OUTDIR}/circtools/star/{{sample}}/paired/{{sample}}_Aligned.sortedByCoord.out.bam",
        bai=f"{OUTDIR}/circtools/star/{{sample}}/paired/{{sample}}_Aligned.sortedByCoord.out.bam.bai",
        unmapped_r1=f"{OUTDIR}/circtools/star/{{sample}}/paired/{{sample}}_Unmapped.out.mate1.gz",
        unmapped_r2=f"{OUTDIR}/circtools/star/{{sample}}/paired/{{sample}}_Unmapped.out.mate2.gz"
    log:
        "logs/circtools/star_paired/{sample}.log"
    threads: int(config["threads"]["star"])
    conda:
        "envs/circtools_star.yaml"
    params:
        index=lambda wc, input: str(input.idx_sa)[:-3],
        prefix=lambda wc, output: str(output.chimeric_junction).removesuffix("Chimeric.out.junction")
    shell:
        r"""
        set -euo pipefail
        out_dir=$(dirname "{output.chimeric_junction}")
        log_dir=$(dirname "{log}")
        mkdir -p "$out_dir" "$log_dir"

        STAR \
          --runThreadN {threads} \
          --genomeDir "{params.index}" \
          --outSAMtype BAM SortedByCoordinate \
          --readFilesIn "{input.r1}" "{input.r2}" \
          --readFilesCommand zcat \
          --outFileNamePrefix "{params.prefix}" \
          --outReadsUnmapped Fastx \
          --outSJfilterOverhangMin 15 15 15 15 \
          --alignSJoverhangMin 15 \
          --alignSJDBoverhangMin 15 \
          --outFilterMultimapNmax 20 \
          --outFilterScoreMin 1 \
          --outFilterMatchNmin 1 \
          --outFilterMismatchNmax 2 \
          --chimSegmentMin 15 \
          --chimScoreMin 15 \
          --chimScoreSeparation 10 \
          --chimJunctionOverhangMin 15 \
          > "{log}" 2>&1

        gzip -f "{params.prefix}Unmapped.out.mate1" 2>> "{log}"
        gzip -f "{params.prefix}Unmapped.out.mate2" 2>> "{log}"

        samtools index -@ {threads} "{output.bam}" 2>> "{log}"
        """


# ---------------------------------------------------------------------------
# STAR 2/3: mate1 single-end re-mapping. Matches official DCC mate1 command
# verbatim (including --seedSearchStartLmax 30 and --outSAMtype None).
# ---------------------------------------------------------------------------
rule circtools_star_mate1:
    input:
        idx_sa=f"{config['reference']['star_index']}/SA",
        unmapped_r1=f"{OUTDIR}/circtools/star/{{sample}}/paired/{{sample}}_Unmapped.out.mate1.gz"
    output:
        chimeric_junction=f"{OUTDIR}/circtools/star/{{sample}}/mate1/{{sample}}_mate1_Chimeric.out.junction",
        log_final=f"{OUTDIR}/circtools/star/{{sample}}/mate1/{{sample}}_mate1_Log.final.out"
    log:
        "logs/circtools/star_mate1/{sample}.log"
    threads: int(config["threads"]["star"])
    conda:
        "envs/circtools_star.yaml"
    params:
        index=lambda wc, input: str(input.idx_sa)[:-3],
        prefix=lambda wc, output: str(output.chimeric_junction).removesuffix("Chimeric.out.junction")
    shell:
        r"""
        set -euo pipefail
        out_dir=$(dirname "{output.chimeric_junction}")
        log_dir=$(dirname "{log}")
        mkdir -p "$out_dir" "$log_dir"

        STAR \
          --runThreadN {threads} \
          --genomeDir "{params.index}" \
          --outSAMtype None \
          --readFilesIn "{input.unmapped_r1}" \
          --readFilesCommand zcat \
          --outFileNamePrefix "{params.prefix}" \
          --outReadsUnmapped Fastx \
          --outSJfilterOverhangMin 15 15 15 15 \
          --alignSJoverhangMin 15 \
          --alignSJDBoverhangMin 15 \
          --seedSearchStartLmax 30 \
          --outFilterMultimapNmax 20 \
          --outFilterScoreMin 1 \
          --outFilterMatchNmin 1 \
          --outFilterMismatchNmax 2 \
          --chimSegmentMin 15 \
          --chimScoreMin 15 \
          --chimScoreSeparation 10 \
          --chimJunctionOverhangMin 15 \
          > "{log}" 2>&1
        """


# ---------------------------------------------------------------------------
# STAR 3/3: mate2 single-end re-mapping. Same as mate1, R2 input.
# ---------------------------------------------------------------------------
rule circtools_star_mate2:
    input:
        idx_sa=f"{config['reference']['star_index']}/SA",
        unmapped_r2=f"{OUTDIR}/circtools/star/{{sample}}/paired/{{sample}}_Unmapped.out.mate2.gz"
    output:
        chimeric_junction=f"{OUTDIR}/circtools/star/{{sample}}/mate2/{{sample}}_mate2_Chimeric.out.junction",
        log_final=f"{OUTDIR}/circtools/star/{{sample}}/mate2/{{sample}}_mate2_Log.final.out"
    log:
        "logs/circtools/star_mate2/{sample}.log"
    threads: int(config["threads"]["star"])
    conda:
        "envs/circtools_star.yaml"
    params:
        index=lambda wc, input: str(input.idx_sa)[:-3],
        prefix=lambda wc, output: str(output.chimeric_junction).removesuffix("Chimeric.out.junction")
    shell:
        r"""
        set -euo pipefail
        out_dir=$(dirname "{output.chimeric_junction}")
        log_dir=$(dirname "{log}")
        mkdir -p "$out_dir" "$log_dir"

        STAR \
          --runThreadN {threads} \
          --genomeDir "{params.index}" \
          --outSAMtype None \
          --readFilesIn "{input.unmapped_r2}" \
          --readFilesCommand zcat \
          --outFileNamePrefix "{params.prefix}" \
          --outReadsUnmapped Fastx \
          --outSJfilterOverhangMin 15 15 15 15 \
          --alignSJoverhangMin 15 \
          --alignSJDBoverhangMin 15 \
          --seedSearchStartLmax 30 \
          --outFilterMultimapNmax 20 \
          --outFilterScoreMin 1 \
          --outFilterMatchNmin 1 \
          --outFilterMismatchNmax 2 \
          --chimSegmentMin 15 \
          --chimScoreMin 15 \
          --chimScoreSeparation 10 \
          --chimJunctionOverhangMin 15 \
          > "{log}" 2>&1
        """


# ---------------------------------------------------------------------------
# circtools detect (DCC). Same as before - flags driven by helpers.smk config.
# ---------------------------------------------------------------------------
rule circtools_detect:
    input:
        install_marker=CIRCTOOLS_INSTALL_MARKER,
        paired_junctions=expand(
            f"{OUTDIR}/circtools/star/{{sample}}/paired/{{sample}}_Chimeric.out.junction",
            sample=CIRCTOOLS_SAMPLES,
        ),
        paired_bams=expand(
            f"{OUTDIR}/circtools/star/{{sample}}/paired/{{sample}}_Aligned.sortedByCoord.out.bam",
            sample=CIRCTOOLS_SAMPLES,
        ),
        paired_bais=expand(
            f"{OUTDIR}/circtools/star/{{sample}}/paired/{{sample}}_Aligned.sortedByCoord.out.bam.bai",
            sample=CIRCTOOLS_SAMPLES,
        ),
        mate1_junctions=expand(
            f"{OUTDIR}/circtools/star/{{sample}}/mate1/{{sample}}_mate1_Chimeric.out.junction",
            sample=CIRCTOOLS_SAMPLES,
        ),
        mate2_junctions=expand(
            f"{OUTDIR}/circtools/star/{{sample}}/mate2/{{sample}}_mate2_Chimeric.out.junction",
            sample=CIRCTOOLS_SAMPLES,
        ),
        gtf=config["reference"]["gtf"],
        fasta=config["reference"]["fasta"]
    output:
        circ_count=f"{OUTDIR}/circtools/detect/CircRNACount",
        circ_coord=f"{OUTDIR}/circtools/detect/CircCoordinates",
        linear_count=f"{OUTDIR}/circtools/detect/LinearCount"
    log:
        "logs/circtools/detect.log"
    threads: int(config["threads"].get("circtools_detect", config["threads"]["star"]))
    conda:
        "envs/circtools.yaml"
    params:
        outdir=lambda wc, output: os.path.dirname(str(output.circ_count)),
        non_strand_flag=("-N" if not CIRCTOOLS_DETECT_STRANDED else ""),
        secondstrand_flag=("-ss" if CIRCTOOLS_DETECT_SS else ""),
        filter_flag=("-F" if CIRCTOOLS_DETECT_FILTER else ""),
        chrm_flag=("-M" if CIRCTOOLS_DETECT_CHRM else ""),
        filterbygene_flag=("-fg" if CIRCTOOLS_DETECT_FILTER_BY_GENE else ""),
        gene_flag=("-G" if CIRCTOOLS_DETECT_HOST_GENE else ""),
        rep_flag=(
            f'-R "{CIRCTOOLS_DETECT_REPEATS_GTF}"'
            if CIRCTOOLS_DETECT_REPEATS_GTF else ""
        ),
        refseq_flag=(
            f'-A "{config["reference"]["fasta"]}"'
            if CIRCTOOLS_DETECT_HOST_GENE else ""
        ),
        count_threshold=CIRCTOOLS_DETECT_COUNT_THRESHOLD,
        replicate_threshold=CIRCTOOLS_DETECT_REPLICATE_THRESHOLD
    shell:
        r"""
        set -euo pipefail
        outdir="{params.outdir}"
        log_dir=$(dirname "{log}")
        mkdir -p "$outdir" "$log_dir"
        tmpdir=$(mktemp -d "$outdir/tmp.detect.XXXXXX")
        trap 'rm -rf "$tmpdir"' EXIT

        paired_sheet="$tmpdir/samplesheet_paired.txt"
        mate1_sheet="$tmpdir/samplesheet_mate1.txt"
        mate2_sheet="$tmpdir/samplesheet_mate2.txt"
        bam_sheet="$tmpdir/bam_file_list.txt"
        : > "$paired_sheet"; : > "$mate1_sheet"; : > "$mate2_sheet"; : > "$bam_sheet"
        for jct in {input.paired_junctions}; do realpath "$jct" >> "$paired_sheet"; done
        for jct in {input.mate1_junctions}; do realpath "$jct" >> "$mate1_sheet"; done
        for jct in {input.mate2_junctions}; do realpath "$jct" >> "$mate2_sheet"; done
        for bam in {input.paired_bams}; do realpath "$bam" >> "$bam_sheet"; done

        circtools detect "@$paired_sheet" -D -T {threads} -O "$outdir/" -t "$tmpdir/detect_work/" -Pi -mt1 "@$mate1_sheet" -mt2 "@$mate2_sheet" -B "@$bam_sheet" {params.gene_flag} {params.refseq_flag} -an "{input.gtf}" {params.rep_flag} {params.non_strand_flag} {params.secondstrand_flag} {params.filter_flag} {params.chrm_flag} {params.filterbygene_flag} -Nr {params.count_threshold} {params.replicate_threshold} > "{log}" 2>&1

        test -s "{output.circ_count}" || (echo "Missing CircRNACount" >> "{log}"; exit 1)
        test -s "{output.circ_coord}" || (echo "Missing CircCoordinates" >> "{log}"; exit 1)
        test -s "{output.linear_count}" || (echo "Missing LinearCount" >> "{log}"; exit 1)
        """


# ---------------------------------------------------------------------------
# CIRI3 -> DCC synthesis (only runs when detect arm is off).
# ---------------------------------------------------------------------------
rule ciri3_to_dcc_dataset:
    input:
        bsj=f"{OUTDIR}/ciri3/all_samples.ciri3.BSJ_Matrix",
        fsj=f"{OUTDIR}/ciri3/all_samples.ciri3.FSJ_Matrix",
        merged_ciri3=f"{OUTDIR}/ciri3/all_samples.ciri3",
        script=os.path.join(workflow.basedir, "rules", "scripts", "ciri3_to_dcc_format.py")
    output:
        circ_count=f"{OUTDIR}/circtools/dcc_from_ciri3/dataset/CircRNACount",
        circ_coord=f"{OUTDIR}/circtools/dcc_from_ciri3/dataset/CircCoordinates",
        linear_count=f"{OUTDIR}/circtools/dcc_from_ciri3/dataset/LinearCount"
    log:
        "logs/circtools/ciri3_to_dcc_dataset.log"
    conda:
        "envs/circtools.yaml"
    params:
        outdir=lambda wc, output: os.path.dirname(str(output.circ_count)),
        samples_csv=",".join(CIRCTOOLS_SAMPLES)
    shell:
        r"""
        set -euo pipefail
        outdir="{params.outdir}"
        log_dir=$(dirname "{log}")
        mkdir -p "$outdir" "$log_dir"
        python "{input.script}" \
          --bsj-matrix "{input.bsj}" \
          --fsj-matrix "{input.fsj}" \
          --ciri3-annot "{input.merged_ciri3}" \
          --samples "{params.samples_csv}" \
          --out-dir "$outdir" \
          > "{log}" 2>&1
        """


# ---------------------------------------------------------------------------
# Per-comparison DCC subset.
# ---------------------------------------------------------------------------
rule circtools_dcc_subset:
    input:
        circ_count=CIRCTOOLS_CANONICAL_DCC_CIRC_COUNT,
        circ_coord=CIRCTOOLS_CANONICAL_DCC_CIRC_COORD,
        linear_count=CIRCTOOLS_CANONICAL_DCC_LINEAR_COUNT,
        script=os.path.join(workflow.basedir, "rules", "scripts", "subset_dcc_for_comparison.py")
    output:
        circ_count=f"{OUTDIR}/circtools/dcc_subset/{{comparison}}/CircRNACount",
        circ_coord=f"{OUTDIR}/circtools/dcc_subset/{{comparison}}/CircCoordinates",
        linear_count=f"{OUTDIR}/circtools/dcc_subset/{{comparison}}/LinearCount"
    log:
        "logs/circtools/dcc_subset_{comparison}.log"
    conda:
        "envs/circtools.yaml"
    params:
        in_dir=lambda wc, input: os.path.dirname(str(input.circ_count)),
        out_dir=lambda wc, output: os.path.dirname(str(output.circ_count)),
        samples_csv=lambda wc: ",".join(
            circtools_circtest_case_samples(wc.comparison)
            + circtools_circtest_control_samples(wc.comparison)
        )
    wildcard_constraints:
        comparison=CIRCTOOLS_CIRCTEST_COMPARISON_REGEX
    shell:
        r"""
        set -euo pipefail
        out_dir="{params.out_dir}"
        log_dir=$(dirname "{log}")
        mkdir -p "$out_dir" "$log_dir"
        python "{input.script}" \
          --in-dir "{params.in_dir}" \
          --out-dir "$out_dir" \
          --samples "{params.samples_csv}" \
          > "{log}" 2>&1
        """


# ---------------------------------------------------------------------------
# STAR Log.final.out staging for quickcheck.
# ---------------------------------------------------------------------------
rule circtools_quickcheck_star_log:
    input:
        log_final=lambda wc: circtools_source_star_log(wc.sample)
    output:
        log_final=f"{OUTDIR}/circtools/quickcheck_staging/star/{{sample}}/Log.final.out"
    log:
        "logs/circtools/quickcheck_star_log_{sample}.log"
    threads: 1
    run:
        import os
        out_dir = os.path.dirname(output.log_final)
        os.makedirs(out_dir, exist_ok=True)
        os.makedirs(os.path.dirname(log[0]), exist_ok=True)
        target = os.path.relpath(os.path.abspath(input.log_final), start=out_dir)
        if os.path.lexists(output.log_final):
            os.remove(output.log_final)
        os.symlink(target, output.log_final)
        with open(log[0], "w") as fh:
            fh.write(f"linked {input.log_final} -> {output.log_final}\n")


# ---------------------------------------------------------------------------
# circtest per comparison.
# ---------------------------------------------------------------------------
rule circtools_circtest:
    input:
        install_marker=CIRCTOOLS_INSTALL_MARKER,
        circ_count=f"{OUTDIR}/circtools/dcc_subset/{{comparison}}/CircRNACount",
        circ_coord=f"{OUTDIR}/circtools/dcc_subset/{{comparison}}/CircCoordinates",
        linear_count=f"{OUTDIR}/circtools/dcc_subset/{{comparison}}/LinearCount"
    output:
        csv=f"{OUTDIR}/circtools/circtest/{{comparison}}/circtest.csv",
        xlsx=f"{OUTDIR}/circtools/circtest/{{comparison}}/circtest.xlsx",
        pdf=f"{OUTDIR}/circtools/circtest/{{comparison}}/circtest.pdf"
    log:
        "logs/circtools/circtest_{comparison}.log"
    threads: 1
    conda:
        "envs/circtools.yaml"
    params:
        detect_dir=lambda wc, input: os.path.dirname(str(input.circ_count)) + "/",
        outdir=lambda wc, output: os.path.dirname(str(output.csv)),
        condition_list=lambda wc: ",".join(circtools_circtest_condition_names(wc.comparison)),
        condition_columns=lambda wc: circtools_circtest_condition_columns(wc.comparison),
        grouping=lambda wc: circtools_circtest_grouping(wc.comparison),
        num_replicates=lambda wc: circtools_circtest_num_replicates(wc.comparison),
        max_fdr=CIRCTOOLS_CIRCTEST_MAX_FDR,
        percentage=CIRCTOOLS_CIRCTEST_PERCENTAGE,
        filter_sample=CIRCTOOLS_CIRCTEST_FILTER_SAMPLE,
        filter_count=CIRCTOOLS_CIRCTEST_FILTER_COUNT,
        max_plots=CIRCTOOLS_CIRCTEST_MAX_PLOTS
    wildcard_constraints:
        comparison=CIRCTOOLS_CIRCTEST_COMPARISON_REGEX
    shell:
        r"""
        set -euo pipefail
        outdir="{params.outdir}"
        log_dir=$(dirname "{log}")
        mkdir -p "$outdir" "$log_dir"

        set +e
        circtools circtest \
          -d "{params.detect_dir}" \
          -l "{params.condition_list}" \
          -c "{params.condition_columns}" \
          -g "{params.grouping}" \
          -r {params.num_replicates} \
          -f {params.max_fdr} \
          -p {params.percentage} \
          -s {params.filter_sample} \
          -C {params.filter_count} \
          -m {params.max_plots} \
          -o "$outdir" \
          -n circtest \
          > "{log}" 2>&1
        circtest_status=$?
        set -e

        if grep -Fq "No candidates to plot, exiting." "{log}"; then
          python "workflow/rules/scripts/write_empty_circtest_outputs.py" \
            --csv "{output.csv}" \
            --xlsx "{output.xlsx}" \
            --pdf "{output.pdf}" \
            --comparison "{wildcards.comparison}" \
            --log "{log}"
          exit 0
        fi

        exit "$circtest_status"
        """


# ---------------------------------------------------------------------------
# primex candidate IDs.
# ---------------------------------------------------------------------------
rule circtools_primex_ids:
    input:
        circ_count=CIRCTOOLS_CANONICAL_DCC_CIRC_COUNT,
        circ_coord=CIRCTOOLS_CANONICAL_DCC_CIRC_COORD,
        script=os.path.join(workflow.basedir, "rules", "scripts", "select_primex_candidates.py")
    output:
        ids=f"{OUTDIR}/circtools/primex/candidate_ids.txt"
    log:
        "logs/circtools/primex_ids.log"
    conda:
        "envs/circtools.yaml"
    params:
        top_n=CIRCTOOLS_PRIMEX_TOP_N
    shell:
        r"""
        set -euo pipefail
        log_dir=$(dirname "{log}")
        out_dir=$(dirname "{output.ids}")
        mkdir -p "$log_dir" "$out_dir"
        python "{input.script}" \
          --circ-count "{input.circ_count}" \
          --circ-coord "{input.circ_coord}" \
          --top-n {params.top_n} \
          --out "{output.ids}" \
          > "{log}" 2>&1
        """


# ---------------------------------------------------------------------------
# circtools primex.
# ---------------------------------------------------------------------------
rule circtools_primex:
    input:
        install_marker=CIRCTOOLS_INSTALL_MARKER,
        circ_coord=CIRCTOOLS_CANONICAL_DCC_CIRC_COORD,
        ids=f"{OUTDIR}/circtools/primex/candidate_ids.txt",
        gtf=config["reference"]["gtf"],
        fasta=config["reference"]["fasta"]
    output:
        marker=f"{OUTDIR}/circtools/primex/primex.done"
    log:
        "logs/circtools/primex.log"
    threads: 1
    conda:
        "envs/circtools.yaml"
    params:
        outdir=lambda wc, output: os.path.dirname(str(output.marker)),
        title=CIRCTOOLS_PRIMEX_TITLE,
        product_size_low=CIRCTOOLS_PRIMEX_PRODUCT_SIZE[0],
        product_size_high=CIRCTOOLS_PRIMEX_PRODUCT_SIZE[1],
        num_pairs=CIRCTOOLS_PRIMEX_NUM_PAIRS,
        junction=CIRCTOOLS_PRIMEX_JUNCTION,
        no_blast_flag=("-b" if CIRCTOOLS_PRIMEX_NO_BLAST else ""),
        organism_flag=(
            f"-O {CIRCTOOLS_PRIMEX_ORGANISM}" if CIRCTOOLS_PRIMEX_ORGANISM else ""
        )
    shell:
        r"""
        set -euo pipefail
        outdir="{params.outdir}"
        log_dir=$(dirname "{log}")
        mkdir -p "$outdir" "$log_dir"

        if [ ! -s "{input.ids}" ]; then
          echo "No primex candidate IDs were selected." > "{log}"
          exit 1
        fi
        ids_args=$(awk 'NF>0 {{printf "%s ", $1}}' "{input.ids}")

        tmpdir=$(mktemp -d "$outdir/tmp.primex.XXXXXX")
        trap 'rm -rf "$tmpdir"' EXIT

        circtools primex \
          -d "{input.circ_coord}" \
          -g "{input.gtf}" \
          -f "{input.fasta}" \
          {params.organism_flag} \
          -o "$outdir" \
          -T "{params.title}" \
          -t "$tmpdir" \
          -p {params.product_size_low} {params.product_size_high} \
          -j {params.junction} \
          -n {params.num_pairs} \
          {params.no_blast_flag} \
          -i $ids_args \
          > "{log}" 2>&1

        date -u +"%Y-%m-%dT%H:%M:%SZ" > "{output.marker}"
        """


# ---------------------------------------------------------------------------
# quickcheck per comparison.
# ---------------------------------------------------------------------------
rule circtools_quickcheck:
    input:
        install_marker=CIRCTOOLS_INSTALL_MARKER,
        circ_count=f"{OUTDIR}/circtools/dcc_subset/{{comparison}}/CircRNACount",
        circ_coord=f"{OUTDIR}/circtools/dcc_subset/{{comparison}}/CircCoordinates",
        linear_count=f"{OUTDIR}/circtools/dcc_subset/{{comparison}}/LinearCount",
        star_logs=lambda wc: [
            f"{OUTDIR}/circtools/quickcheck_staging/star/{sample}/Log.final.out"
            for sample in (
                circtools_circtest_case_samples(wc.comparison)
                + circtools_circtest_control_samples(wc.comparison)
            )
        ]
    output:
        pdf=f"{OUTDIR}/circtools/quickcheck/{{comparison}}/quickcheck.pdf"
    log:
        "logs/circtools/quickcheck_{comparison}.log"
    threads: 1
    conda:
        "envs/circtools.yaml"
    params:
        detect_dir=lambda wc, input: os.path.dirname(str(input.circ_count)) + "/",
        star_dir=f"{OUTDIR}/circtools/quickcheck_staging/star/",
        outdir=lambda wc, output: os.path.dirname(str(output.pdf)),
        condition_list=lambda wc: ",".join(circtools_circtest_condition_names(wc.comparison)),
        grouping=lambda wc: circtools_circtest_grouping(wc.comparison)
    wildcard_constraints:
        comparison=CIRCTOOLS_CIRCTEST_COMPARISON_REGEX
    shell:
        r"""
        set -euo pipefail
        outdir="{params.outdir}"
        log_dir=$(dirname "{log}")
        mkdir -p "$outdir" "$log_dir"

        circtools quickcheck \
          -d "{params.detect_dir}" \
          -s "{params.star_dir}" \
          -l "{params.condition_list}" \
          -g "{params.grouping}" \
          -o "$outdir" \
          -n quickcheck \
          > "{log}" 2>&1
        """
