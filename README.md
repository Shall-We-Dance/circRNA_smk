# Circular RNA Snakemake Analysis Pipeline

This repository provides a reproducible Snakemake workflow for analyzing circular RNA from paired-end ribo-depleted RNA-seq data. The pipeline merges lane-level FASTQs to sample-level pairs, performs sample-level QC with fastp, generates MultiQC summaries, aligns reads to a reference genome with a prebuilt STAR index, performs BWA remapping for CIRI3, and performs downstream circRNA detection and quantification.

## Overview

### Inputs
- Paired-end raw FASTQ files (`*.fastq.gz`) for each sample.
- Multiple FASTQ pairs (lanes / technical replicates) may be assigned to the same sample.

### Outputs
- **QC**
  - `results/qc/fastp/<sample>/fastp.html/json` (per-sample)
  - `results/qc/multiqc/multiqc_report.html`
- **Alignment**
  - `results/star/<sample>/<sample>.Aligned.sortedByCoord.out.bam` (optional, controlled by `output.keep_bam`)
  - `results/star/<sample>/<sample>.Aligned.sortedByCoord.out.bam.bai` (optional, controlled by `output.keep_bam`)
  - `results/star/<sample>/<sample>.bwa.bam` (optional, controlled by `output.keep_bam`)
- **Quantification / circRNA**
  - `results/featurecount/totalRNA.counts.txt`
  - `results/ciri3/star/<sample>.ciri3`
  - `results/ciri3/star/<sample>.ciri3.BSJ_Matrix`
  - `results/ciri3/star/<sample>.ciri3.FSJ_Matrix`
  - `results/ciri3/all_samples.ciri3`
  - `results/ciri3/all_samples.ciri3.BSJ_Matrix`
  - `results/ciri3/all_samples.ciri3.FSJ_Matrix`
- **BSJ differential expression (optional, DESeq2)**
  - `results/deg/bsj/sample_metadata.tsv`
  - `results/deg/bsj/all_groups/deseq2_results.tsv`
  - `results/deg/bsj/all_groups/volcano.pdf`
  - `results/deg/bsj/all_groups/heatmap_top50.pdf`
  - `results/deg/bsj/all_groups/pca.pdf`
  - `results/deg/bsj/all_groups/vst_counts.tsv`
  - `results/deg/bsj/pairwise/<GroupA>_vs_<GroupB>/deseq2_results.tsv`
  - `results/deg/bsj/pairwise/<GroupA>_vs_<GroupB>/volcano.pdf`
  - `results/deg/bsj/pairwise/<GroupA>_vs_<GroupB>/heatmap_top50.pdf`
### Intermediate file handling
To minimize storage footprint, intermediate FASTQs produced by fastp and merged FASTQs are marked as temporary and are removed automatically by Snakemake. In addition, STAR/BWA BAM files are temporary by default (`output.keep_bam: false`) and will be removed after downstream rules finish. Set `output.keep_bam: true` if you want to retain BAM/BAI files. Final deliverables include:
- QC reports (fastp + multiqc)
- CIRI3 per-sample outputs and all-sample merged outputs
- featureCounts count matrices

## Pipeline steps

1. **Merge FASTQs (per-sample)**  
   Raw R1 files are concatenated into a single sample-level R1 file, and raw R2 files are concatenated into a single sample-level R2 file. (For gzip-compressed FASTQ files, stream concatenation with `cat` is valid.)

2. **fastp QC (per-sample)**  
   The merged FASTQs are processed by fastp once per sample. This produces sample-level HTML/JSON reports that are consistent with the reads used for alignment. The cleaned FASTQs from this step are temporary.

3. **STAR alignment**  
   Reads are aligned with STAR and a coordinate-sorted BAM is generated directly by STAR.

4. **BWA remapping for CIRI3**  
   Unmapped STAR mates are remapped by BWA to generate the CIRI3-required BWA BAM.

5. **circRNA quantification and aggregation**  
   Gene-level quantification is generated with featureCounts and circular RNA detection is run with CIRI3. The workflow writes per-sample CIRI3 outputs first, then merges all samples into `all_samples.ciri3`, `all_samples.ciri3.BSJ_Matrix`, and `all_samples.ciri3.FSJ_Matrix` with sample names as matrix column names.

6. **BSJ differential expression (optional)**  
   When `deg.enabled: true`, the pipeline uses `all_samples.ciri3.BSJ_Matrix` as count input and performs:
   - a full multi-group DESeq2 model (`design = ~ group`) to obtain overall BSJ DEG signals,
   - all pairwise group comparisons,
   - visualization outputs for each analysis (volcano plots, heatmaps, PCA).

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
* `output.keep_bam`: keep STAR/BWA BAM files (`false` by default to save disk)
* `deg.enabled`: run optional BSJ-level DE analysis (`false` by default)
* `deg.groups`: group-to-sample mapping for DESeq2 analysis

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
  sampleB:
    R1:
      - "raw/sampleB_R1.fastq.gz"
    R2:
      - "raw/sampleB_R2.fastq.gz"

deg:
  enabled: true
  groups:
    GroupA:
      - sampleA
    GroupB:
      - sampleB
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

For questions, please open an issue or contact the pipeline maintainer.
