# Circular RNA Snakemake Analysis Pipeline

This repository provides a reproducible Snakemake workflow for analyzing circular RNA from paired-end ribo-depleted RNA-seq data. The pipeline performs per-lane QC with fastp, generates MultiQC summaries, aligns reads to a reference genome with a prebuilt STAR index, performs BWA remapping for CIRI3, and performs downstream circRNA detection and quantification.

## Overview

### Inputs
- Paired-end raw FASTQ files (`*.fastq.gz`) for each sample.
- Multiple FASTQ pairs (lanes / technical replicates) may be assigned to the same sample.

### Outputs (per sample)
- **QC**
  - `results/qc/fastp/<sample>/unit<k>.html/json` (per-lane)
  - `results/qc/fastp/<sample>/fastp.html/json` (per-sample; computed after merge)
  - `results/qc/multiqc/multiqc_report.html`
- **Alignment**
  - `results/star/<sample>/<sample>.Aligned.sortedByCoord.out.bam`
  - `results/star/<sample>/<sample>.Aligned.sortedByCoord.out.bam.bai`
  - `results/star/<sample>/<sample>.bwa.bam`
- **Quantification / circRNA**
  - `results/featurecount/totalRNA.counts.txt`
  - `results/ciri3/all_samples.ciri3`
  - `results/ciri3/all_samples.ciri3.BSJ_Matrix`
  - `results/ciri3/all_samples.ciri3.FSJ_Matrix`
### Intermediate file handling
To minimize storage footprint, intermediate FASTQs produced by fastp and merged FASTQs are marked as temporary and are removed automatically by Snakemake. Final deliverables include:
- QC reports (fastp + multiqc)
- Final STAR BAM + BAI
- CIRI3 outputs and featureCounts count matrices

## Pipeline steps

1. **fastp QC (per-lane)**  
   Each input FASTQ pair is processed by fastp. Per-lane HTML/JSON reports are kept; cleaned FASTQs are temporary.

2. **Merge cleaned FASTQs**  
   Cleaned R1 files are concatenated to a single sample-level R1; cleaned R2 files are concatenated to a single sample-level R2 (gzip stream concatenation is valid).

3. **fastp (per-sample “report-only” pass)**  
   The merged FASTQs are passed through fastp with trimming/filtering disabled to produce a **sample-level report** consistent with the data used for alignment. Output FASTQs from this step are also temporary.

4. **STAR alignment**  
   Reads are aligned with STAR and a coordinate-sorted BAM is generated directly by STAR.

5. **BWA remapping for CIRI3**  
   Unmapped STAR mates are remapped by BWA to generate the CIRI3-required BWA BAM.

6. **circRNA quantification and aggregation**  
   Gene-level quantification is generated with featureCounts and circular RNA detection is run with CIRI3, producing per-sample outputs and merged summary matrices.

## Requirements

- **Snakemake** (recommended: Snakemake >= 7)
- Conda / Mamba (recommended for environment management)
- STAR
- samtools
- fastp
- MultiQC
- CIRI3 (download from https://github.com/gyjames/CIRI3 ; this workflow uses the Java 18 build)
- Python packages: `pysam`

All dependencies are provided via the conda environments under `workflow/rules/envs/` (including OpenJDK 18 for CIRI3).

## Installation

Recommended: create environments on-the-fly via Snakemake.

Example:
```bash
snakemake -s workflow/Snakefile --use-conda --cores 16
```

To speed up conda solves, consider using mamba:

```bash
snakemake -s workflow/Snakefile --use-conda --conda-frontend mamba --cores 16
```

## Configuration

Edit `config.yaml`.

Key fields:

* `reference.star_index`: prebuilt STAR genome index directory (must already contain index files such as `SA`)
* `reference.fasta`: reference FASTA (used by downstream steps such as CIRI3 support files)
* `reference.bwa_indexed_fasta`: BWA-indexed FASTA path (BWA sidecar files must already exist with this path as prefix)
* `reference.gtf`: reference annotation GTF
* `samples`: mapping of sample name to lists of FASTQs for R1 and R2

Example:

```yaml
reference:
  star_index: "/path/to/STAR/index"
  fasta: "/path/to/genome.fa"
  bwa_indexed_fasta: "/path/to/genome.fa"
  gtf: "/path/to/genes.gtf"

samples:
  sampleA:
    R1:
      - "raw/sampleA_L001_R1.fastq.gz"
      - "raw/sampleA_L002_R1.fastq.gz"
    R2:
      - "raw/sampleA_L001_R2.fastq.gz"
      - "raw/sampleA_L002_R2.fastq.gz"
```

Notes:

* The workflow assumes paired-end reads and requires both R1 and R2 lists to be the same length per sample.
* STAR genome index construction is **not** performed in this workflow; build the STAR index in advance and point `reference.star_index` to that directory.
* BWA index construction is **not** performed in this workflow; provide prebuilt index files for `reference.bwa_indexed_fasta`.
  Example:
  `bwa index /path/to/genome.fa`

## Running the workflow

```bash
snakemake -s workflow/Snakefile --use-conda --cores 16
```

Dry run:

```bash
snakemake -s workflow/Snakefile -n
```

## Contact

For questions or extensions (e.g., track hubs, fragment-level scaling, additional QC metrics), please open an issue or contact the pipeline maintainer.
