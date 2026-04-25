# workflow/rules/qc_fastp.smk

rule merge_raw_fastq_per_sample:
    input:
        r1=lambda wc: config["samples"][wc.sample]["R1"],
        r2=lambda wc: config["samples"][wc.sample]["R2"]
    output:
        merged_r1=temp(f"{OUTDIR}/tmp/merged_raw/{{sample}}_R1.fastq.gz"),
        merged_r2=temp(f"{OUTDIR}/tmp/merged_raw/{{sample}}_R2.fastq.gz")
    log:
        "logs/merge_fastq/{sample}.log"
    threads: int(config["threads"].get("merge_fastq", 2))
    conda:
        "envs/qc.yaml"
    shell:
        r"""
        set -euo pipefail
        output_dir=$(dirname "{output.merged_r1}")
        log_dir=$(dirname "{log}")
        mkdir -p "$output_dir" "$log_dir"
        cat {input.r1:q} > "{output.merged_r1}" 2> "{log}"
        cat {input.r2:q} > "{output.merged_r2}" 2>> "{log}"
        """


rule fastp_sample_level:
    input:
        r1=f"{OUTDIR}/tmp/merged_raw/{{sample}}_R1.fastq.gz",
        r2=f"{OUTDIR}/tmp/merged_raw/{{sample}}_R2.fastq.gz"
    output:
        clean_r1=temp(f"{OUTDIR}/tmp/fastp/{{sample}}_R1.fastq.gz"),
        clean_r2=temp(f"{OUTDIR}/tmp/fastp/{{sample}}_R2.fastq.gz"),
        html=f"{OUTDIR}/qc/fastp/{{sample}}/fastp.html",
        json=f"{OUTDIR}/qc/fastp/{{sample}}/fastp.json"
    log:
        f"logs/fastp/{{sample}}.log"
    threads: int(config["threads"]["fastp"])
    conda:
        "envs/qc.yaml"
    params:
        dedup_arg=("--dedup" if bool(config.get("fastp", {}).get("dedup_adapter", {}).get("dedup", True)) else "")
    shell:
        r"""
        set -euo pipefail
        tmp_dir=$(dirname "{output.clean_r1}")
        html_dir=$(dirname "{output.html}")
        log_dir=$(dirname "{log}")
        mkdir -p "$tmp_dir" "$html_dir" "$log_dir"

        fastp \
          -i "{input.r1}" -I "{input.r2}" \
          -o "{output.clean_r1}" -O "{output.clean_r2}" \
          --thread {threads} \
          {params.dedup_arg} \
          --html "{output.html}" --json "{output.json}" \
          > "{log}" 2>&1
        """


rule multiqc:
    input:
        expand(f"{OUTDIR}/qc/fastp/{{sample}}/fastp.html", sample=list(config["samples"].keys())),
        expand(f"{OUTDIR}/qc/star/{{sample}}/{{sample}}.Log.final.out", sample=list(config["samples"].keys())),
    output:
        html=f"{OUTDIR}/qc/multiqc/multiqc_report.html"
    log:
        "logs/multiqc.log"
    conda:
        "envs/multiqc.yaml"
    params:
        outdir=lambda wc, output: os.path.dirname(output.html),
        qc_dir=lambda wc, output: os.path.dirname(os.path.dirname(output.html))
    shell:
        r"""
        set -euo pipefail
        html_dir=$(dirname "{output.html}")
        log_dir=$(dirname "{log}")
        mkdir -p "$html_dir" "$log_dir"
        multiqc --force -o "{params.outdir}" "{params.qc_dir}" > "{log}" 2>&1
        """
