OUTDIR = config["output"]["dir"]
SAMPLES = list(config["samples"].keys())


rule compute_splicing_site_stats:
    input:
        ciri3=f"{OUTDIR}/ciri3/per_sample/{{sample}}.ciri3",
        bsj=f"{OUTDIR}/ciri3/per_sample/{{sample}}.ciri3.BSJ_Matrix",
        fsj=f"{OUTDIR}/ciri3/per_sample/{{sample}}.ciri3.FSJ_Matrix"
    output:
        circ_table=f"{OUTDIR}/splicing/{{sample}}/circ_splice_sites.tsv",
        summary=f"{OUTDIR}/splicing/{{sample}}/summary.tsv",
        dist=f"{OUTDIR}/splicing/{{sample}}/distributions.tsv",
        plot=f"{OUTDIR}/splicing/{{sample}}/distribution.png"
    threads: int(config["threads"].get("splicing_stats", 2))
    conda:
        "envs/splicing_stats.yaml"
    script:
        "scripts/compute_splicing_site_stats.py"


rule summarize_splicing_site_stats:
    input:
        summaries=expand(f"{OUTDIR}/splicing/{{sample}}/summary.tsv", sample=SAMPLES),
        distributions=expand(f"{OUTDIR}/splicing/{{sample}}/distributions.tsv", sample=SAMPLES)
    output:
        summary=f"{OUTDIR}/splicing/all_samples_summary.tsv",
        distributions=f"{OUTDIR}/splicing/all_samples_distributions.tsv",
        plot=f"{OUTDIR}/splicing/all_samples_overview.png"
    params:
        samples=SAMPLES
    threads: int(config["threads"].get("splicing_stats", 2))
    conda:
        "envs/splicing_stats.yaml"
    script:
        "scripts/summarize_splicing_site_stats.py"
