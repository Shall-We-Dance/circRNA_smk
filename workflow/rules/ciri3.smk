# workflow/rules/ciri3.smk

OUTDIR = config["output"]["dir"]
SAMPLES = list(config["samples"].keys())


def maybe_temp(path):
    return path if KEEP_BAM else temp(path)


rule star_bam_index_for_ciri3:
    input:
        bam=f"{OUTDIR}/star/{{sample}}/{{sample}}.Aligned.sortedByCoord.out.bam"
    output:
        bai=maybe_temp(f"{OUTDIR}/star/{{sample}}/{{sample}}.Aligned.sortedByCoord.out.bam.bai")
    threads: int(config["threads"]["samtools"])
    conda:
        "envs/samtools.yaml"
    shell:
        r"""
        set -euo pipefail
        samtools index -@ {threads} {input.bam}
        """


rule ciri3_detect_star:
    input:
        chimeric=f"{OUTDIR}/star/{{sample}}/{{sample}}.Chimeric.out.junction",
        bam=f"{OUTDIR}/star/{{sample}}/{{sample}}.Aligned.sortedByCoord.out.bam",
        bwa=f"{OUTDIR}/star/{{sample}}/{{sample}}.bwa.bam",
        fasta=config["reference"]["fasta"],
        gtf=config["reference"]["gtf"],
        bai=f"{OUTDIR}/star/{{sample}}/{{sample}}.Aligned.sortedByCoord.out.bam.bai"
    output:
        result=f"{OUTDIR}/ciri3/star/{{sample}}.ciri3",
        bsj=f"{OUTDIR}/ciri3/star/{{sample}}.ciri3.BSJ_Matrix",
        fsj=f"{OUTDIR}/ciri3/star/{{sample}}.ciri3.FSJ_Matrix"
    log:
        "logs/ciri3/{sample}.detect.star.log"
    threads: int(config["threads"].get("samtools", 1))
    conda:
        "envs/ciri3.yaml"
    params:
        jar=config["ciri3"]["jar"],
        Ma=int(config.get("ciri3", {}).get("Ma", 1)),
        W=int(config.get("ciri3", {}).get("W", 1))
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.result}) $(dirname {log})
        tmpdir=$(mktemp -d "$(dirname {output.result})/ciri3_{wildcards.sample}.XXXXXX")
        trap 'rm -rf "$tmpdir"' EXIT
        star_sam="$tmpdir/{wildcards.sample}.Aligned.sortedByCoord.out.sam"
        bwa_sam="$tmpdir/{wildcards.sample}.bwa.sam"
        samples_tsv="$tmpdir/samples_{wildcards.sample}.tsv"

        samtools view -@ {threads} -h -o "$star_sam" {input.bam}
        samtools view -@ {threads} -h -o "$bwa_sam" {input.bwa}

        chimeric_abs=$(realpath "{input.chimeric}")
        star_sam_abs=$(realpath "$star_sam")
        bwa_sam_abs=$(realpath "$bwa_sam")
        printf "%s,%s,%s\n" "$chimeric_abs" "$star_sam_abs" "$bwa_sam_abs" > "$samples_tsv"

        java -jar {params.jar} \
          -A {input.gtf} \
          -Ma {params.Ma} \
          -W {params.W} \
          -I "$samples_tsv" \
          -O {output.result} \
          -F {input.fasta} \
          > {log} 2>&1

        # CIRI3 naming differs across releases; normalize to expected targets.
        if [ ! -e {output.result} ] && [ -e {output.result}.ciri3 ]; then cp {output.result}.ciri3 {output.result}; fi
        if [ ! -e {output.bsj} ] && [ -e {output.result}.ciri3.BSJ_Matrix ]; then cp {output.result}.ciri3.BSJ_Matrix {output.bsj}; fi
        if [ ! -e {output.fsj} ] && [ -e {output.result}.ciri3.FSJ_Matrix ]; then cp {output.result}.ciri3.FSJ_Matrix {output.fsj}; fi

        test -s {output.result}
        test -s {output.bsj}
        test -s {output.fsj}
        """


rule merge_ciri3_outputs:
    input:
        ciri3=expand(f"{OUTDIR}/ciri3/star/{{sample}}.ciri3", sample=SAMPLES),
        bsj=expand(f"{OUTDIR}/ciri3/star/{{sample}}.ciri3.BSJ_Matrix", sample=SAMPLES),
        fsj=expand(f"{OUTDIR}/ciri3/star/{{sample}}.ciri3.FSJ_Matrix", sample=SAMPLES)
    output:
        ciri3=f"{OUTDIR}/ciri3/all_samples.ciri3",
        bsj=f"{OUTDIR}/ciri3/all_samples.ciri3.BSJ_Matrix",
        fsj=f"{OUTDIR}/ciri3/all_samples.ciri3.FSJ_Matrix"
    params:
        samples=SAMPLES
    script:
        "scripts/merge_ciri3_outputs.py"


rule bsj_motif_analysis:
    input:
        bsj=f"{OUTDIR}/ciri3/all_samples.ciri3.BSJ_Matrix",
        ciri3=f"{OUTDIR}/ciri3/all_samples.ciri3",
        fasta=config["reference"]["fasta"]
    output:
        site_table=f"{OUTDIR}/motif/bsj_sites.tsv",
        motif_summary=f"{OUTDIR}/motif/bsj_motif_summary.tsv"
    conda:
        "envs/py_signal.yaml"
    params:
        flank=int(config.get("motif", {}).get("flank", 30)),
        kmer=int(config.get("motif", {}).get("kmer", 4)),
        weight_mode=str(config.get("motif", {}).get("weight_mode", "sum"))
    script:
        "scripts/bsj_motif_analysis.py"
