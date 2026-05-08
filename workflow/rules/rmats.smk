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
          echo "Building rMATS-turbo with ./build_rmats {params.build_args}" >> "{log}"
          repo_dir_abs=$(cd "$repo_dir" && pwd)
          cd "$repo_dir_abs"
          ./build_rmats {params.build_args} >> "$log_abs" 2>&1
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

        echo "Converting rmats2sashimiplot source for Python 3 compatibility" >> "{log}"
        python -m lib2to3 -w -n "$repo_dir/src/rmats2sashimiplot" >> "{log}" 2>&1

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
        tmpdir="$outdir/tmp"
        log_dir=$(dirname "{log}")
        mkdir -p "$outdir" "$tmpdir" "$log_dir"
        workdir_abs=$(pwd)
        rmats_script_abs=$(realpath "{input.rmats_script}")
        rmats_dir=$(dirname "$rmats_script_abs")
        b1_abs=$(realpath "{input.b1}")
        b2_abs=$(realpath "{input.b2}")
        gtf_abs=$(realpath "{input.gtf}")
        outdir_abs=$(realpath "$outdir")
        tmpdir_abs=$(realpath "$tmpdir")
        log_abs="$(cd "$log_dir" && pwd)/$(basename "{log}")"

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
          > "$log_abs" 2>&1

        cd "$workdir_abs"
        for event_file in {output.events:q}; do
          test -s "$event_file" || (echo "Missing rMATS event file: $event_file" >> "$log_abs"; exit 1)
        done
        test -s "{output.summary}" || (echo "Missing rMATS summary: {output.summary}" >> "$log_abs"; exit 1)
        """


rule rmats_prepare_sashimi_groups:
    input:
        bams=sashimi_all_condition_bams,
        bais=sashimi_all_condition_bais
    output:
        b1=f"{OUTDIR}/rmats/sashimi/all_conditions/b1.txt",
        b2=f"{OUTDIR}/rmats/sashimi/all_conditions/b2.txt",
        group_info=f"{OUTDIR}/rmats/sashimi/all_conditions/grouping.gf",
        samples=f"{OUTDIR}/rmats/sashimi/all_conditions/samples.tsv"
    params:
        group_names=SASHIMI_GROUP_NAMES,
        groups=DEG_GROUPS,
        reference_group=SASHIMI_REFERENCE_GROUP,
        b1_samples=SASHIMI_B1_SAMPLES,
        b2_samples=SASHIMI_B2_SAMPLES
    log:
        "logs/rmats/sashimi.prepare_groups.log"
    script:
        "scripts/prepare_rmats_sashimi_groups.py"


rule rmats2sashimi_plot_events:
    input:
        script=RMATS2SASHIMI_SCRIPT,
        ready=RMATS2SASHIMI_READY,
        event_file=f"{OUTDIR}/rmats/{{comparison}}/{{event_type}}.MATS.{{count_type}}.txt",
        b1=f"{OUTDIR}/rmats/sashimi/all_conditions/b1.txt",
        b2=f"{OUTDIR}/rmats/sashimi/all_conditions/b2.txt",
        group_info=f"{OUTDIR}/rmats/sashimi/all_conditions/grouping.gf",
        bams=sashimi_all_condition_bams,
        bais=sashimi_all_condition_bais
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
        label_b1=SASHIMI_LABEL_B1,
        label_b2=SASHIMI_LABEL_B2,
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
