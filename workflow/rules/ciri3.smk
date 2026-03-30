# workflow/rules/ciri3.smk

OUTDIR = config["output"]["dir"]


rule star_bam_index_for_ciri3:
    input:
        bam=f"{OUTDIR}/star/{{sample}}/{{sample}}.Aligned.sortedByCoord.out.bam"
    output:
        bai=f"{OUTDIR}/star/{{sample}}/{{sample}}.Aligned.sortedByCoord.out.bam.bai"
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
        bwa=f"{OUTDIR}/star/{{sample}}/{{sample}}.bwa.unique.mapq11.sorted.bam",
        fasta=config["reference"]["fasta"],
        gtf=config["reference"]["gtf"],
        bai=f"{OUTDIR}/star/{{sample}}/{{sample}}.Aligned.sortedByCoord.out.bam.bai"
    output:
        result=f"{OUTDIR}/ciri3/star/{{sample}}.ciri3",
        bsj=f"{OUTDIR}/ciri3/star/{{sample}}.ciri3.BSJ_Matrix",
        fsj=f"{OUTDIR}/ciri3/star/{{sample}}.ciri3.FSJ_Matrix"
    log:
        "logs/ciri3/{sample}.detect.star.log"
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
        input_list="{input.chimeric},{input.bam},{input.bwa}"
        java -jar {params.jar} \
          -A {input.gtf} \
          -Ma {params.Ma} \
          -W {params.W} \
          -I "$input_list" \
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


rule ciri3_detect_unique:
    input:
        bam=f"{OUTDIR}/star/{{sample}}/{{sample}}.unique.mapq11.sorted.bam",
        fasta=config["reference"]["fasta"],
        gtf=config["reference"]["gtf"],
        bai=f"{OUTDIR}/star/{{sample}}/{{sample}}.unique.mapq11.sorted.bam.bai"
    output:
        result=f"{OUTDIR}/ciri3/unique/{{sample}}.ciri3",
        bsj=f"{OUTDIR}/ciri3/unique/{{sample}}.ciri3.BSJ_Matrix",
        fsj=f"{OUTDIR}/ciri3/unique/{{sample}}.ciri3.FSJ_Matrix"
    log:
        "logs/ciri3/{sample}.detect.unique.log"
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
        input_bam="{input.bam}"
        java -jar {params.jar} \
          -A {input.gtf} \
          -Ma {params.Ma} \
          -W {params.W} \
          -I "$input_bam" \
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
