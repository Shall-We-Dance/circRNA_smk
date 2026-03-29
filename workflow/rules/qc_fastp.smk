# workflow/rules/qc_fastp.smk
import os

OUTDIR = config["output"]["dir"]


def units(sample):
    return list(range(len(config["samples"][sample]["R1"])))


rule merge_raw_fastq_per_sample:
    input:
        r1=lambda wc: [config["samples"][wc.sample]["R1"][i] for i in units(wc.sample)],
        r2=lambda wc: [config["samples"][wc.sample]["R2"][i] for i in units(wc.sample)]
    output:
        merged_r1=temp(f"{OUTDIR}/tmp/merged_raw/{{sample}}_R1.fastq.gz"),
        merged_r2=temp(f"{OUTDIR}/tmp/merged_raw/{{sample}}_R2.fastq.gz")
    threads: 2
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.merged_r1})
        cat {input.r1} > {output.merged_r1}
        cat {input.r2} > {output.merged_r2}
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
        dedup_arg=("--dedup" if bool(config.get("fastp", {}).get("dedup_adapter", {}).get("dedup", True)) else ""),
        disable_adapter_arg=("--disable_adapter_trimming" if bool(config.get("fastp", {}).get("grip_trim", {}).get("disable_adapter_trimming", False)) else ""),
        trim_f=int(config.get("fastp", {}).get("grip_trim", {}).get("trim_front_r1", 0)),
        trim_t=int(config.get("fastp", {}).get("grip_trim", {}).get("trim_tail_r1", 0)),
        trim_F=int(config.get("fastp", {}).get("grip_trim", {}).get("trim_front_r2", 0)),
        trim_T=int(config.get("fastp", {}).get("grip_trim", {}).get("trim_tail_r2", 0))
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.clean_r1}) $(dirname {output.html}) $(dirname {log})

        fastp \
          -i {input.r1} -I {input.r2} \
          -o {output.clean_r1} -O {output.clean_r2} \
          --thread {threads} \
          {params.dedup_arg} \
          {params.disable_adapter_arg} \
          -f {params.trim_f} -t {params.trim_t} \
          -F {params.trim_F} -T {params.trim_T} \
          --html {output.html} --json {output.json} \
          > {log} 2>&1
        """


rule multiqc:
    input:
        expand(f"{OUTDIR}/qc/fastp/{{sample}}/fastp.html", sample=list(config["samples"].keys())),
        expand(f"{OUTDIR}/qc/star/{{sample}}/{{sample}}.Log.final.out", sample=list(config["samples"].keys()))
    output:
        html=f"{OUTDIR}/qc/multiqc/multiqc_report.html"
    log:
        "logs/multiqc.log"
    conda:
        "envs/multiqc.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.html}) $(dirname {log})
        multiqc --force -o {OUTDIR}/qc/multiqc {OUTDIR}/qc {OUTDIR}/star > {log} 2>&1
        """
