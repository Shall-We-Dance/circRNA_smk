# workflow/rules/motif.smk

OUTDIR = config["output"]["dir"]


rule prepare_bsj_motif_inputs:
    input:
        bsj=f"{OUTDIR}/ciri3/star/{{sample}}.ciri3.BSJ_Matrix",
        ciri3=f"{OUTDIR}/ciri3/star/{{sample}}.ciri3",
        fasta=config["reference"]["fasta"]
    output:
        site_table=f"{OUTDIR}/motif/{{sample}}/bsj_sites.tsv",
        site_fasta=f"{OUTDIR}/motif/{{sample}}/bsj_sites.fa"
    conda:
        "envs/motif.yaml"
    params:
        flank=int(config.get("motif", {}).get("flank", 30)),
        weight_mode=str(config.get("motif", {}).get("weight_mode", "sum"))
    script:
        "scripts/prepare_bsj_motif_inputs.py"


rule run_homer_bsj_motif:
    input:
        site_table=f"{OUTDIR}/motif/{{sample}}/bsj_sites.tsv",
        site_fasta=f"{OUTDIR}/motif/{{sample}}/bsj_sites.fa"
    output:
        motif_summary=f"{OUTDIR}/motif/{{sample}}/bsj_motif_summary.tsv",
        homer_dir=directory(f"{OUTDIR}/motif/{{sample}}/homer"),
        known_results=f"{OUTDIR}/motif/{{sample}}/homer/knownResults.txt"
    log:
        "logs/motif/{sample}.homer.log"
    threads: int(config["threads"].get("samtools", 1))
    conda:
        "envs/motif.yaml"
    params:
        homer_len=str(config.get("motif", {}).get("homer_len", "8,10,12"))
    shell:
        r"""
        set -euo pipefail
        mkdir -p {output.homer_dir} $(dirname {log})

        findMotifs.pl {input.site_fasta} fasta {output.homer_dir} \
          -len {params.homer_len} \
          -p {threads} \
          > {log} 2>&1

        if [ -s {output.known_results} ]; then
          cp {output.known_results} {output.motif_summary}
        elif [ -s {output.homer_dir}/homerMotifs.all.motifs ]; then
          printf "# HOMER completed, but knownResults.txt is missing.\n" > {output.known_results}
          printf "# Please inspect: %s\n" "{output.homer_dir}/homerMotifs.all.motifs" >> {output.known_results}
          printf "# HOMER completed, but knownResults.txt is missing.\n" > {output.motif_summary}
          printf "# Please inspect: %s\n" "{output.homer_dir}/homerMotifs.all.motifs" >> {output.motif_summary}
        else
          echo "HOMER did not generate expected result files." >&2
          exit 1
        fi
        """
