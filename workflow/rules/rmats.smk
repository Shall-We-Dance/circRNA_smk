# workflow/rules/rmats.smk

rule fetch_rmats_turbo:
    output:
        script=RMATS_SCRIPT,
        ready=RMATS_READY
    params:
        repo_url=RMATS_REPO_URL,
        ref=RMATS_REPO_REF,
        build_args=(
            ("--no-darts-model" if not _get_bool(rmats_cfg, "darts_model", False) else "")
            + ("" if RMATS_PAIRED_STATS else " --no-paired-model")
        ).strip()
    log:
        "logs/rmats/fetch_rmats_turbo.log"
    conda:
        "envs/rmats_turbo.yaml"
    shell:
        r"""
        set -euo pipefail
        log_dir=$(dirname "{log}")
        repo_dir=$(dirname "{output.ready}")
        repo_parent=$(dirname "$repo_dir")
        mkdir -p "$log_dir" "$repo_parent"
        log_abs="$(cd "$log_dir" && pwd)/$(basename "{log}")"
        workdir_abs=$(pwd)

        if [ -z "$repo_dir" ] || [ "$repo_dir" = "." ] || [ "$repo_dir" = "$(dirname "$repo_dir")" ]; then
          echo "Unsafe rMATS-turbo install_dir: $repo_dir" > "{log}"
          exit 1
        fi

        if [ -d "$repo_dir/.git" ]; then
          echo "Updating rMATS-turbo repository in $repo_dir" > "{log}"
          git -C "$repo_dir" remote set-url origin "{params.repo_url}" >> "{log}" 2>&1
          git -C "$repo_dir" fetch --tags --force origin >> "{log}" 2>&1
        elif [ -e "$repo_dir" ]; then
          if [ -s "{output.script}" ] && [ -e "$repo_dir/.rmats_build_complete" ]; then
            echo "Using existing rMATS-turbo install_dir with required files: $repo_dir" > "{log}"
            touch "{output.ready}"
            exit 0
          fi

          echo "rMATS-turbo install_dir exists but is not a complete checkout: $repo_dir" > "{log}"
          stamp=$(date +%Y%m%d%H%M%S)
          backup_dir="${{repo_dir}}.stale.${{stamp}}"
          suffix=0
          while [ -e "$backup_dir" ]; do
            suffix=$((suffix + 1))
            backup_dir="${{repo_dir}}.stale.${{stamp}}.${{suffix}}"
          done
          mv "$repo_dir" "$backup_dir"
          echo "Moved existing rMATS-turbo install_dir to backup: $backup_dir" >> "{log}"
          echo "Cloning rMATS-turbo repository from {params.repo_url} into $repo_dir" >> "{log}"
          git clone "{params.repo_url}" "$repo_dir" >> "{log}" 2>&1
          git -C "$repo_dir" fetch --tags --force origin >> "{log}" 2>&1
        else
          echo "Cloning rMATS-turbo repository from {params.repo_url} into $repo_dir" > "{log}"
          git clone "{params.repo_url}" "$repo_dir" >> "{log}" 2>&1
          git -C "$repo_dir" fetch --tags --force origin >> "{log}" 2>&1
        fi

        if ! git -C "$repo_dir" checkout --force "{params.ref}" >> "{log}" 2>&1; then
          echo "Could not check out rMATS-turbo ref '{params.ref}'." >> "{log}"
          exit 1
        fi

        if [ ! -e "$repo_dir/.rmats_build_complete" ]; then
          echo "Building rMATS-turbo with bash ./build_rmats {params.build_args}" >> "{log}"
          repo_dir_abs=$(cd "$repo_dir" && pwd)
          cd "$repo_dir_abs"
          # Upstream rMATS-turbo ships build scripts with CRLF line endings on
          # some clones, which makes the kernel reject the shebang
          # ("/bin/bash^M: bad interpreter"). Strip carriage returns from the
          # build/helper shell scripts, then invoke through bash explicitly so
          # a corrupted shebang cannot break the build.
          for _f in build_rmats build_rmats.sh test_rmats setup_environment.sh; do
            if [ -f "$_f" ]; then
              sed -i 's/\r$//' "$_f" >> "$log_abs" 2>&1 || true
            fi
          done
          find . -type f -name "*.sh" -exec sed -i 's/\r$//' {{}} + >> "$log_abs" 2>&1 || true
          bash ./build_rmats {params.build_args} >> "$log_abs" 2>&1
          touch "$repo_dir_abs/.rmats_build_complete"
          cd "$workdir_abs"
        fi

        test -s "{output.script}" || (echo "Missing rmats.py after checkout/build: {output.script}" >> "{log}"; exit 1)
        touch "{output.ready}"
        """


rule fetch_rmats2sashimiplot:
    output:
        script=RMATS2SASHIMI_SCRIPT,
        ready=RMATS2SASHIMI_READY
    params:
        repo_url=RMATS2SASHIMI_REPO_URL,
        ref=RMATS2SASHIMI_REPO_REF
    log:
        "logs/rmats/fetch_rmats2sashimiplot.log"
    conda:
        "envs/rmats2sashimiplot.yaml"
    shell:
        r"""
        set -euo pipefail
        log_dir=$(dirname "{log}")
        repo_dir=$(dirname "{output.ready}")
        repo_parent=$(dirname "$repo_dir")
        mkdir -p "$log_dir" "$repo_parent"

        if [ -z "$repo_dir" ] || [ "$repo_dir" = "." ] || [ "$repo_dir" = "$(dirname "$repo_dir")" ]; then
          echo "Unsafe rmats2sashimiplot install_dir: $repo_dir" > "{log}"
          exit 1
        fi

        if [ -d "$repo_dir/.git" ]; then
          echo "Updating rmats2sashimiplot repository in $repo_dir" > "{log}"
          git -C "$repo_dir" remote set-url origin "{params.repo_url}" >> "{log}" 2>&1
          git -C "$repo_dir" fetch --tags --force origin >> "{log}" 2>&1
        elif [ -e "$repo_dir" ]; then
          if [ -s "{output.script}" ]; then
            echo "Using existing rmats2sashimiplot install_dir with required files: $repo_dir" > "{log}"
          else
            echo "rmats2sashimiplot install_dir exists but is not a complete checkout: $repo_dir" > "{log}"
            stamp=$(date +%Y%m%d%H%M%S)
            backup_dir="${{repo_dir}}.stale.${{stamp}}"
            suffix=0
            while [ -e "$backup_dir" ]; do
              suffix=$((suffix + 1))
              backup_dir="${{repo_dir}}.stale.${{stamp}}.${{suffix}}"
            done
            mv "$repo_dir" "$backup_dir"
            echo "Moved existing rmats2sashimiplot install_dir to backup: $backup_dir" >> "{log}"
            echo "Cloning rmats2sashimiplot repository from {params.repo_url} into $repo_dir" >> "{log}"
            git clone "{params.repo_url}" "$repo_dir" >> "{log}" 2>&1
            git -C "$repo_dir" fetch --tags --force origin >> "{log}" 2>&1
          fi
        else
          echo "Cloning rmats2sashimiplot repository from {params.repo_url} into $repo_dir" > "{log}"
          git clone "{params.repo_url}" "$repo_dir" >> "{log}" 2>&1
          git -C "$repo_dir" fetch --tags --force origin >> "{log}" 2>&1
        fi

        if [ -d "$repo_dir/.git" ] && ! git -C "$repo_dir" checkout --force "{params.ref}" >> "{log}" 2>&1; then
          echo "Could not check out rmats2sashimiplot ref '{params.ref}'." >> "{log}"
          exit 1
        fi

        echo "Converting bundled rmats2sashimiplot/MISO source for Python 3 compatibility" >> "{log}"
        python -m lib2to3 -w -n "$repo_dir/src/rmats2sashimiplot" "$repo_dir/src/MISO" >> "{log}" 2>&1

        test -s "{output.script}" || (echo "Missing rmats2sashimiplot.py after checkout: {output.script}" >> "{log}"; exit 1)
        touch "{output.ready}"
        """


rule rmats_prepare_pairwise_inputs:
    input:
        bams=rmats_comparison_bams
    output:
        b1=f"{OUTDIR}/rmats/{{comparison}}/inputs/b1.txt",
        b2=f"{OUTDIR}/rmats/{{comparison}}/inputs/b2.txt",
        samples=f"{OUTDIR}/rmats/{{comparison}}/inputs/samples.tsv"
    wildcard_constraints:
        comparison=RMATS_COMPARISON_REGEX
    params:
        case_samples=lambda wc: rmats_comparison_case_samples(wc.comparison),
        control_samples=lambda wc: rmats_comparison_control_samples(wc.comparison),
        case_group=lambda wc: DEG_COMPARISONS[wc.comparison]["case_group"],
        control_group=lambda wc: DEG_COMPARISONS[wc.comparison]["control_group"]
    log:
        "logs/rmats/{comparison}.prepare_inputs.log"
    conda:
        "envs/rmats2sashimiplot.yaml"
    script:
        "scripts/prepare_rmats_pairwise_inputs.py"


rule rmats_turbo_pairwise:
    input:
        rmats_script=RMATS_SCRIPT,
        rmats_ready=RMATS_READY,
        b1=f"{OUTDIR}/rmats/{{comparison}}/inputs/b1.txt",
        b2=f"{OUTDIR}/rmats/{{comparison}}/inputs/b2.txt",
        bams=rmats_comparison_bams,
        gtf=config["reference"]["gtf"]
    output:
        events=expand(
            f"{OUTDIR}/rmats/{{{{comparison}}}}/{{event_type}}.MATS.{{count_type}}.txt",
            event_type=RMATS_EVENT_TYPES,
            count_type=RMATS_COUNT_TYPES,
        ),
        summary=f"{OUTDIR}/rmats/{{comparison}}/summary.txt"
    wildcard_constraints:
        comparison=RMATS_COMPARISON_REGEX
    params:
        read_type=RMATS_READ_TYPE,
        read_length=RMATS_READ_LENGTH,
        lib_type=RMATS_LIB_TYPE,
        task=RMATS_TASK,
        flags=" ".join(
            [
                "--variable-read-length" if RMATS_VARIABLE_READ_LENGTH else "",
                "--allow-clipping" if RMATS_ALLOW_CLIPPING else "",
                "--novelSS" if RMATS_NOVEL_SS else "",
                "--paired-stats" if RMATS_PAIRED_STATS else "",
                "--statoff" if RMATS_STATOFF else "",
                "--individual-counts" if RMATS_INDIVIDUAL_COUNTS else "",
                RMATS_EXTRA_ARGS,
            ]
        ).strip()
    log:
        "logs/rmats/{comparison}.rmats_turbo.log"
    threads: int(config["threads"].get("rmats", config["threads"].get("star", 1)))
    conda:
        "envs/rmats_turbo.yaml"
    shell:
        r"""
        set -euo pipefail
        outdir=$(dirname "{output.summary}")
        log_dir=$(dirname "{log}")
        mkdir -p "$outdir" "$log_dir"
        workdir_abs=$(pwd)
        rmats_script_abs=$(realpath "{input.rmats_script}")
        rmats_dir=$(dirname "$rmats_script_abs")
        b1_abs=$(realpath "{input.b1}")
        b2_abs=$(realpath "{input.b2}")
        gtf_abs=$(realpath "{input.gtf}")
        outdir_abs=$(realpath "$outdir")
        log_abs="$(cd "$log_dir" && pwd)/$(basename "{log}")"
        tmpdir_abs="$outdir_abs/tmp"

        if [ -z "$outdir_abs" ] || [ "$outdir_abs" = "/" ]; then
          echo "Unsafe rMATS output directory: $outdir_abs" > "$log_abs"
          exit 1
        fi

        echo "Cleaning rMATS tmp directory before run: $tmpdir_abs" > "$log_abs"
        rm -rf "$tmpdir_abs"
        mkdir -p "$tmpdir_abs"

        cd "$rmats_dir"
        python "$rmats_script_abs" \
          --b1 "$b1_abs" \
          --b2 "$b2_abs" \
          --gtf "$gtf_abs" \
          -t "{params.read_type}" \
          --readLength {params.read_length} \
          --nthread {threads} \
          --od "$outdir_abs" \
          --tmp "$tmpdir_abs" \
          --task "{params.task}" \
          --libType "{params.lib_type}" \
          {params.flags} \
          >> "$log_abs" 2>&1

        cd "$workdir_abs"
        for event_file in {output.events:q}; do
          test -s "$event_file" || (echo "Missing rMATS event file: $event_file" >> "$log_abs"; exit 1)
        done
        test -s "{output.summary}" || (echo "Missing rMATS summary: {output.summary}" >> "$log_abs"; exit 1)
        """


rule rmats_prepare_group_inputs:
    input:
        bams=rmats_group_bams
    output:
        b1=f"{OUTDIR}/rmats/groups/{{group}}/inputs/b1.txt",
        samples=f"{OUTDIR}/rmats/groups/{{group}}/inputs/samples.tsv"
    wildcard_constraints:
        group=RMATS_GROUP_REGEX
    params:
        group=lambda wc: wc.group,
        samples=lambda wc: rmats_group_samples(wc.group)
    log:
        "logs/rmats/groups/{group}.prepare_inputs.log"
    conda:
        "envs/rmats2sashimiplot.yaml"
    script:
        "scripts/prepare_rmats_group_inputs.py"


rule rmats_turbo_group:
    input:
        rmats_script=RMATS_SCRIPT,
        rmats_ready=RMATS_READY,
        b1=f"{OUTDIR}/rmats/groups/{{group}}/inputs/b1.txt",
        bams=rmats_group_bams,
        gtf=config["reference"]["gtf"]
    output:
        events=expand(
            f"{OUTDIR}/rmats/groups/{{{{group}}}}/{{event_type}}.MATS.{{count_type}}.txt",
            event_type=RMATS_EVENT_TYPES,
            count_type=RMATS_COUNT_TYPES,
        ),
        summary=f"{OUTDIR}/rmats/groups/{{group}}/summary.txt"
    wildcard_constraints:
        group=RMATS_GROUP_REGEX
    params:
        read_type=RMATS_READ_TYPE,
        read_length=RMATS_READ_LENGTH,
        lib_type=RMATS_LIB_TYPE,
        task=RMATS_TASK,
        flags=" ".join(
            [
                "--variable-read-length" if RMATS_VARIABLE_READ_LENGTH else "",
                "--allow-clipping" if RMATS_ALLOW_CLIPPING else "",
                "--novelSS" if RMATS_NOVEL_SS else "",
                "--statoff",
                "--individual-counts" if RMATS_INDIVIDUAL_COUNTS else "",
                RMATS_EXTRA_ARGS,
            ]
        ).strip()
    log:
        "logs/rmats/groups/{group}.rmats_turbo.log"
    threads: int(config["threads"].get("rmats", config["threads"].get("star", 1)))
    conda:
        "envs/rmats_turbo.yaml"
    shell:
        r"""
        set -euo pipefail
        outdir=$(dirname "{output.summary}")
        log_dir=$(dirname "{log}")
        mkdir -p "$outdir" "$log_dir"
        workdir_abs=$(pwd)
        rmats_script_abs=$(realpath "{input.rmats_script}")
        rmats_dir=$(dirname "$rmats_script_abs")
        b1_abs=$(realpath "{input.b1}")
        gtf_abs=$(realpath "{input.gtf}")
        outdir_abs=$(realpath "$outdir")
        log_abs="$(cd "$log_dir" && pwd)/$(basename "{log}")"
        tmpdir_abs="$outdir_abs/tmp"

        if [ -z "$outdir_abs" ] || [ "$outdir_abs" = "/" ]; then
          echo "Unsafe rMATS output directory: $outdir_abs" > "$log_abs"
          exit 1
        fi

        echo "Cleaning rMATS tmp directory before run: $tmpdir_abs" > "$log_abs"
        rm -rf "$tmpdir_abs"
        mkdir -p "$tmpdir_abs"

        cd "$rmats_dir"
        python "$rmats_script_abs" \
          --b1 "$b1_abs" \
          --gtf "$gtf_abs" \
          -t "{params.read_type}" \
          --readLength {params.read_length} \
          --nthread {threads} \
          --od "$outdir_abs" \
          --tmp "$tmpdir_abs" \
          --task "{params.task}" \
          --libType "{params.lib_type}" \
          {params.flags} \
          >> "$log_abs" 2>&1

        cd "$workdir_abs"
        for event_file in {output.events:q}; do
          test -s "$event_file" || (echo "Missing rMATS event file: $event_file" >> "$log_abs"; exit 1)
        done
        test -s "{output.summary}" || (echo "Missing rMATS summary: {output.summary}" >> "$log_abs"; exit 1)
        """


if SASHIMI_GFF3_FROM_GTF:
    rule sashimi_gtf_to_gff3:
        input:
            gtf=config["reference"]["gtf"]
        output:
            gff3=SASHIMI_GFF3
        log:
            "logs/rmats/sashimi/prepare_gff3.log"
        conda:
            "envs/rmats2sashimiplot.yaml"
        script:
            "scripts/gtf_to_gff3_for_sashimi.py"


rule rmats_prepare_sashimi_all_groups:
    input:
        bams=sashimi_all_group_bams,
        bais=sashimi_all_group_bais
    output:
        b1=f"{OUTDIR}/rmats/sashimi/bsj/inputs/b1.txt",
        b2=f"{OUTDIR}/rmats/sashimi/bsj/inputs/b2.txt",
        group_info=f"{OUTDIR}/rmats/sashimi/bsj/inputs/grouping.gf",
        samples=f"{OUTDIR}/rmats/sashimi/bsj/inputs/samples.tsv"
    params:
        group_names=DEG_GROUP_NAMES,
        groups=DEG_GROUPS
    log:
        "logs/rmats/sashimi/bsj.prepare_groups.log"
    conda:
        "envs/rmats2sashimiplot.yaml"
    script:
        "scripts/prepare_rmats_sashimi_all_groups.py"


rule rmats2sashimi_plot_bsj_events:
    input:
        script=RMATS2SASHIMI_SCRIPT,
        ready=RMATS2SASHIMI_READY,
        result=bsj_sashimi_result_path,
        ciri3=f"{OUTDIR}/ciri3/all_samples.ciri3",
        bsj_matrix=f"{OUTDIR}/ciri3/all_samples.ciri3.BSJ_Matrix",
        fsj_matrix=f"{OUTDIR}/ciri3/all_samples.ciri3.FSJ_Matrix",
        gff3=SASHIMI_GFF3,
        b1=f"{OUTDIR}/rmats/sashimi/bsj/inputs/b1.txt",
        b2=f"{OUTDIR}/rmats/sashimi/bsj/inputs/b2.txt",
        group_info=f"{OUTDIR}/rmats/sashimi/bsj/inputs/grouping.gf",
        samples=f"{OUTDIR}/rmats/sashimi/bsj/inputs/samples.tsv",
        bams=sashimi_all_group_bams,
        bais=sashimi_all_group_bais
    output:
        manifest=f"{OUTDIR}/rmats/sashimi/bsj/{{method}}/{{comparison}}/manifest.tsv",
        done=f"{OUTDIR}/rmats/sashimi/bsj/{{method}}/{{comparison}}/plots.done",
        plots=directory(f"{OUTDIR}/rmats/sashimi/bsj/{{method}}/{{comparison}}/plots"),
        bsj_only_plots=directory(f"{OUTDIR}/rmats/sashimi/bsj/{{method}}/{{comparison}}/plots_bsj_only")
    wildcard_constraints:
        method=BSJ_SASHIMI_METHOD_REGEX,
        comparison=CIRI3_DE_COMPARISON_REGEX
    params:
        outdir=lambda wc: f"{OUTDIR}/rmats/sashimi/bsj/{wc.method}/{wc.comparison}",
        padj_cutoff=DEG_PADJ_CUTOFF,
        lfc_cutoff=DEG_LFC_CUTOFF,
        max_events=SASHIMI_MAX_EVENTS,
        bsj_flank=SASHIMI_BSJ_FLANK,
        min_counts=SASHIMI_MIN_COUNTS,
        exon_s=SASHIMI_EXON_SCALE,
        intron_s=SASHIMI_INTRON_SCALE,
        auto_scale=SASHIMI_AUTO_SCALE,
        font_size=SASHIMI_FONT_SIZE,
        min_font_size=SASHIMI_MIN_FONT_SIZE,
        fig_width=SASHIMI_FIG_WIDTH,
        fig_height=SASHIMI_FIG_HEIGHT,
        min_fig_width=SASHIMI_MIN_FIG_WIDTH,
        max_fig_width=SASHIMI_MAX_FIG_WIDTH,
        min_fig_height=SASHIMI_MIN_FIG_HEIGHT,
        max_fig_height=SASHIMI_MAX_FIG_HEIGHT,
        label_b1=DEG_GROUP_NAMES[0] if DEG_GROUP_NAMES else "group_1",
        label_b2="other_groups",
        colors=SASHIMI_COLORS,
        fail_on_error=SASHIMI_FAIL_ON_ERROR,
        extra_args=SASHIMI_EXTRA_ARGS
    log:
        "logs/rmats/sashimi/bsj/{method}.{comparison}.log"
    threads: int(config["threads"].get("sashimi", 1))
    conda:
        "envs/rmats2sashimiplot.yaml"
    script:
        "scripts/plot_bsj_sashimi_events.py"


rule rmats_prepare_sashimi_groups:
    input:
        bams=rmats_comparison_bams,
        bais=rmats_comparison_bais
    output:
        b1=f"{OUTDIR}/rmats/sashimi/{{comparison}}/inputs/b1.txt",
        b2=f"{OUTDIR}/rmats/sashimi/{{comparison}}/inputs/b2.txt",
        group_info=f"{OUTDIR}/rmats/sashimi/{{comparison}}/inputs/grouping.gf",
        samples=f"{OUTDIR}/rmats/sashimi/{{comparison}}/inputs/samples.tsv"
    wildcard_constraints:
        comparison=RMATS_COMPARISON_REGEX
    params:
        group_names=lambda wc: sashimi_comparison_group_names(wc.comparison),
        groups=lambda wc: sashimi_comparison_groups(wc.comparison),
        b1_group=lambda wc: DEG_COMPARISONS[wc.comparison]["case_group"],
        b2_group=lambda wc: DEG_COMPARISONS[wc.comparison]["control_group"],
        b1_samples=lambda wc: rmats_comparison_case_samples(wc.comparison),
        b2_samples=lambda wc: rmats_comparison_control_samples(wc.comparison)
    log:
        "logs/rmats/sashimi/{comparison}.prepare_groups.log"
    conda:
        "envs/rmats2sashimiplot.yaml"
    script:
        "scripts/prepare_rmats_sashimi_groups.py"


rule rmats2sashimi_plot_events:
    input:
        script=RMATS2SASHIMI_SCRIPT,
        ready=RMATS2SASHIMI_READY,
        event_file=f"{OUTDIR}/rmats/{{comparison}}/{{event_type}}.MATS.{{count_type}}.txt",
        b1=f"{OUTDIR}/rmats/sashimi/{{comparison}}/inputs/b1.txt",
        b2=f"{OUTDIR}/rmats/sashimi/{{comparison}}/inputs/b2.txt",
        group_info=f"{OUTDIR}/rmats/sashimi/{{comparison}}/inputs/grouping.gf",
        bams=rmats_comparison_bams,
        bais=rmats_comparison_bais
    output:
        manifest=f"{OUTDIR}/rmats/sashimi/{{comparison}}/{{event_type}}.{{count_type}}/manifest.tsv",
        done=f"{OUTDIR}/rmats/sashimi/{{comparison}}/{{event_type}}.{{count_type}}/plots.done"
    wildcard_constraints:
        comparison=RMATS_COMPARISON_REGEX,
        event_type=RMATS_EVENT_TYPE_REGEX,
        count_type=RMATS_COUNT_TYPE_REGEX
    params:
        outdir=lambda wc: f"{OUTDIR}/rmats/sashimi/{wc.comparison}/{wc.event_type}.{wc.count_type}",
        fdr_cutoff=SASHIMI_FDR_CUTOFF,
        pvalue_cutoff=SASHIMI_PVALUE_CUTOFF,
        inc_diff_cutoff=SASHIMI_INC_DIFF_CUTOFF,
        max_events=SASHIMI_MAX_EVENTS,
        min_counts=SASHIMI_MIN_COUNTS,
        exon_s=SASHIMI_EXON_SCALE,
        intron_s=SASHIMI_INTRON_SCALE,
        font_size=SASHIMI_FONT_SIZE,
        fig_width=SASHIMI_FIG_WIDTH,
        fig_height=SASHIMI_FIG_HEIGHT,
        label_b1=sashimi_label_b1,
        label_b2=sashimi_label_b2,
        colors=SASHIMI_COLORS,
        fail_on_error=SASHIMI_FAIL_ON_ERROR,
        keep_event_chr_prefix=SASHIMI_KEEP_EVENT_CHR_PREFIX,
        remove_event_chr_prefix=SASHIMI_REMOVE_EVENT_CHR_PREFIX,
        extra_args=SASHIMI_EXTRA_ARGS
    log:
        "logs/rmats/sashimi/{comparison}.{event_type}.{count_type}.log"
    threads: int(config["threads"].get("sashimi", 1))
    conda:
        "envs/rmats2sashimiplot.yaml"
    script:
        "scripts/plot_rmats_sashimi_events.py"
