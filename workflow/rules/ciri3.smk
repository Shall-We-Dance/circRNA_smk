# workflow/rules/ciri3.smk

rule fetch_ciri3_release:
    output:
        jar=CIRI3_JAR,
        bsj_yes=CIRI3_BSJ_YES,
        ready=CIRI3_READY
    params:
        repo_url=CIRI3_REPO_URL,
        ref=CIRI3_REPO_REF
    log:
        "logs/ciri3/fetch_ciri3.log"
    conda:
        "envs/ciri3.yaml"
    shell:
        r"""
        set -euo pipefail
        log_dir=$(dirname "{log}")
        repo_dir=$(dirname "{output.ready}")
        repo_parent=$(dirname "$repo_dir")
        mkdir -p "$log_dir" "$repo_parent"

        if [ -z "$repo_dir" ] || [ "$repo_dir" = "." ] || [ "$repo_dir" = "$(dirname "$repo_dir")" ]; then
          echo "Unsafe CIRI3 install_dir: $repo_dir" > "{log}"
          exit 1
        fi

        if [ -d "$repo_dir/.git" ]; then
          echo "Updating CIRI3 repository in $repo_dir" > "{log}"
          git -C "$repo_dir" remote set-url origin "{params.repo_url}" >> "{log}" 2>&1
          git -C "$repo_dir" fetch --tags --force origin >> "{log}" 2>&1
        elif [ -e "$repo_dir" ]; then
          if [ -s "{output.jar}" ] && [ -s "{output.bsj_yes}" ]; then
            echo "Using existing non-git CIRI3 install_dir with required files: $repo_dir" > "{log}"
            touch "{output.ready}"
            exit 0
          fi

          echo "CIRI3 install_dir exists but is not a complete git checkout: $repo_dir" > "{log}"
          stamp=$(date +%Y%m%d%H%M%S)
          backup_dir="${{repo_dir}}.stale.${{stamp}}"
          suffix=0
          while [ -e "$backup_dir" ]; do
            suffix=$((suffix + 1))
            backup_dir="${{repo_dir}}.stale.${{stamp}}.${{suffix}}"
          done
          mv "$repo_dir" "$backup_dir"
          echo "Moved existing CIRI3 install_dir to backup: $backup_dir" >> "{log}"
          echo "Cloning CIRI3 repository from {params.repo_url} into $repo_dir" >> "{log}"
          git clone "{params.repo_url}" "$repo_dir" >> "{log}" 2>&1
          git -C "$repo_dir" fetch --tags --force origin >> "{log}" 2>&1
        else
          echo "Cloning CIRI3 repository from {params.repo_url} into $repo_dir" > "{log}"
          git clone "{params.repo_url}" "$repo_dir" >> "{log}" 2>&1
          git -C "$repo_dir" fetch --tags --force origin >> "{log}" 2>&1
        fi

        if ! git -C "$repo_dir" checkout --force "{params.ref}" >> "{log}" 2>&1; then
          echo "Could not check out CIRI3 ref '{params.ref}'." >> "{log}"
          echo "The official GitHub tag corresponding to CIRI-3.0.1 is 'v3.0.1'." >> "{log}"
          exit 1
        fi

        test -s "{output.jar}" || (echo "Missing CIRI3 jar after checkout: {output.jar}" >> "{log}"; exit 1)
        test -s "{output.bsj_yes}" || (echo "Missing CIRI3 DEG script after checkout: {output.bsj_yes}" >> "{log}"; exit 1)
        touch "{output.ready}"
        """


rule star_bam_index_for_ciri3:
    input:
        bam=f"{OUTDIR}/star/{{sample}}/{{sample}}.Aligned.sortedByCoord.out.bam"
    output:
        bai=maybe_temp(f"{OUTDIR}/star/{{sample}}/{{sample}}.Aligned.sortedByCoord.out.bam.bai")
    log:
        "logs/ciri3/{sample}.star_bam_index.log"
    threads: int(config["threads"].get("samtools_index", config["threads"].get("samtools", 1)))
    conda:
        "envs/samtools.yaml"
    shell:
        r"""
        set -euo pipefail
        log_dir=$(dirname "{log}")
        mkdir -p "$log_dir"
        samtools index -@ {threads} "{input.bam}" > "{log}" 2>&1
        """


rule ciri3_detect_star:
    input:
        jar=CIRI3_JAR,
        ciri3_ready=CIRI3_READY,
        chimeric=f"{OUTDIR}/star/{{sample}}/{{sample}}.Chimeric.out.junction",
        bam=f"{OUTDIR}/star/{{sample}}/{{sample}}.Aligned.sortedByCoord.out.bam",
        bwa=f"{OUTDIR}/star/{{sample}}/{{sample}}.bwa.bam",
        fasta=config["reference"]["fasta"],
        gtf=config["reference"]["gtf"],
        bai=f"{OUTDIR}/star/{{sample}}/{{sample}}.Aligned.sortedByCoord.out.bam.bai"
    output:
        result=f"{OUTDIR}/ciri3/per_sample/{{sample}}.ciri3",
        bsj=f"{OUTDIR}/ciri3/per_sample/{{sample}}.ciri3.BSJ_Matrix",
        fsj=f"{OUTDIR}/ciri3/per_sample/{{sample}}.ciri3.FSJ_Matrix"
    log:
        "logs/ciri3/{sample}.detect.star.log"
    threads: int(config["threads"].get("samtools_view", config["threads"].get("samtools", 1)))
    conda:
        "envs/ciri3.yaml"
    params:
        Ma=int(config.get("ciri3", {}).get("Ma", 1)),
        W=int(config.get("ciri3", {}).get("W", 1)),
        normalize_script=CIRI3_NORMALIZE_SCRIPT
    shell:
        r"""
        set -euo pipefail
        output_dir=$(dirname "{output.result}")
        log_dir=$(dirname "{log}")
        mkdir -p "$output_dir" "$log_dir"
        tmpdir=$(mktemp -d "$output_dir/ciri3_{wildcards.sample}.XXXXXX")
        trap 'rm -rf "$tmpdir"' EXIT

        star_sam="$tmpdir/{wildcards.sample}.Aligned.sortedByCoord.out.sam"
        bwa_sam="$tmpdir/{wildcards.sample}.bwa.sam"
        samples_tsv="$tmpdir/samples_{wildcards.sample}.tsv"

        samtools view -@ {threads} -h -o "$star_sam" "{input.bam}"
        samtools view -@ {threads} -h -o "$bwa_sam" "{input.bwa}"

        chimeric_abs=$(realpath "{input.chimeric}")
        star_sam_abs=$(realpath "$star_sam")
        bwa_sam_abs=$(realpath "$bwa_sam")
        printf "%s,%s,%s\n" "$chimeric_abs" "$star_sam_abs" "$bwa_sam_abs" > "$samples_tsv"

        java -jar "{input.jar}" \
          -A "{input.gtf}" \
          -Ma {params.Ma} \
          -W {params.W} \
          -I "$samples_tsv" \
          -O "{output.result}" \
          -F "{input.fasta}" \
          > "{log}" 2>&1

        python "{params.normalize_script}" \
          --sample "{wildcards.sample}" \
          --result "{output.result}" \
          --bsj "{output.bsj}" \
          --fsj "{output.fsj}" \
          --tmpdir "$tmpdir" \
          --log "{log}"

        test -s "{output.result}"
        test -s "{output.bsj}"
        test -s "{output.fsj}"
        """


rule merge_ciri3_outputs:
    input:
        ciri3=expand(f"{OUTDIR}/ciri3/per_sample/{{sample}}.ciri3", sample=SAMPLES),
        bsj=expand(f"{OUTDIR}/ciri3/per_sample/{{sample}}.ciri3.BSJ_Matrix", sample=SAMPLES),
        fsj=expand(f"{OUTDIR}/ciri3/per_sample/{{sample}}.ciri3.FSJ_Matrix", sample=SAMPLES)
    output:
        ciri3=f"{OUTDIR}/ciri3/all_samples.ciri3",
        bsj=f"{OUTDIR}/ciri3/all_samples.ciri3.BSJ_Matrix",
        fsj=f"{OUTDIR}/ciri3/all_samples.ciri3.FSJ_Matrix"
    log:
        "logs/ciri3/merge_outputs.log"
    params:
        samples=SAMPLES
    conda:
        "envs/ciri3.yaml"
    script:
        "scripts/merge_ciri3_outputs.py"
