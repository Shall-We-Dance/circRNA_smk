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
  - `results/ciri3/per_sample/<sample>.ciri3`
  - `results/ciri3/per_sample/<sample>.ciri3.BSJ_Matrix`
  - `results/ciri3/per_sample/<sample>.ciri3.FSJ_Matrix`
  - `results/ciri3/all_samples.ciri3`
  - `results/ciri3/all_samples.ciri3.BSJ_Matrix`
  - `results/ciri3/all_samples.ciri3.FSJ_Matrix`
  - `results/splicing/<sample>/circ_splice_sites.tsv`
  - `results/splicing/<sample>/summary.tsv`
  - `results/splicing/<sample>/distributions.tsv`
  - `results/splicing/<sample>/distribution.png`
  - `results/splicing/all_samples_summary.tsv`
  - `results/splicing/all_samples_distributions.tsv`
  - `results/splicing/all_samples_overview.png`
  - `results/motif/<sample>/bsj_sites.tsv`
  - `results/motif/<sample>/bsj_motif_summary.tsv` (copied from HOMER `knownResults.txt`)
  - `results/motif/<sample>/homer/` (full HOMER output directory)
  - `results/motif/all_samples_site_stats.tsv`
  - `results/motif/all_samples_known_motif_summary.tsv`
  - `results/motif/all_samples_overview.png`
- **BSJ differential expression (optional, DESeq2 + CIRI3 DE modules)**
  - `results/deg/bsj/sample_metadata.tsv`
  - `results/deg/bsj/all_groups/deseq2_results.tsv`
  - `results/deg/bsj/all_groups/heatmap_top50.pdf`
  - `results/deg/bsj/all_groups/pca.pdf`
  - `results/deg/bsj/all_groups/vst_counts.tsv`
  - `results/deg/bsj/pairwise/<CaseGroup>_vs_<ControlGroup>/deseq2_results.tsv`
  - `results/deg/bsj/pairwise/<CaseGroup>_vs_<ControlGroup>/volcano.pdf`
  - `results/deg/bsj/pairwise/<CaseGroup>_vs_<ControlGroup>/heatmap_top50.pdf`
  - `results/deg/ciri3/<CaseGroup>_vs_<ControlGroup>/de_bsj/result.txt`
  - `results/deg/ciri3/<CaseGroup>_vs_<ControlGroup>/de_ratio/result.txt`
  - `results/deg/ciri3/<CaseGroup>_vs_<ControlGroup>/de_relative/result.txt`

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
   The unified `deg:` block controls both DESeq2 and CIRI3 differential modules:
   - `deg.run_deseq2: true` runs DESeq2 BSJ differential analysis from `all_samples.ciri3.BSJ_Matrix` (overall + pairwise with volcano/heatmap/PCA plots),
   - `deg.run_de_bsj: true` runs CIRI3 `DE_BSJ` pairwise analysis (with per-pair `infor.tsv`, BSJ matrix subset, and featureCounts-derived gene expression matrix),
   - `deg.run_de_ratio: true` runs CIRI3 `DE_Ratio` pairwise analysis (with per-pair `infor.tsv`, BSJ subset, and FSJ subset),
   - `deg.run_de_relative: true` runs CIRI3 `DE_Relative` pairwise analysis (with per-pair `infor.tsv` built from sample CIRI3 result paths).
   Pairwise comparisons are generated from `deg.groups` in group order. For `conditionX`, `conditionY`, `conditionZ`, outputs are named `conditionY_vs_conditionX`, `conditionZ_vs_conditionY`, and `conditionZ_vs_conditionX`, where the first group is the case group and the second group is the control group.

7. **Splicing-site feature statistics (enabled by default)**  
   The workflow computes per-sample circRNA splicing-site feature tables using each sample's CIRI3 BSJ/FSJ matrices plus the reference genome FASTA. It reports BSJ/FSJ counts, BSJ-vs-FSJ ratio, BSJ span, and splice-site dinucleotide classes (canonical `GU-AG`, semi-canonical `GC-AG`, minor `AU-AC`, non-canonical, unknown), with per-sample distribution plots, then merges all samples into unified summary/distribution tables and an overview figure.

8. **BSJ motif analysis (optional, enabled by default)**  
   The workflow exports BSJ donor/acceptor sequence windows for each sample into `results/motif/<sample>/bsj_sites.tsv` and `results/motif/<sample>/bsj_sites.fa`, then runs HOMER `findMotifs.pl` per sample. It also performs all-sample motif summary statistics/plots. Outputs include:
   - `results/motif/<sample>/bsj_sites.tsv` (per-BSJ sequence windows),
   - `results/motif/<sample>/bsj_motif_summary.tsv` (copied from HOMER `knownResults.txt`),
   - `results/motif/<sample>/homer/` (full HOMER result directory),
   - `results/motif/all_samples_site_stats.tsv`,
   - `results/motif/all_samples_known_motif_summary.tsv`,
   - `results/motif/all_samples_overview.png`.

## Requirements

- **Snakemake** (recommended: Snakemake >= 7)
- Conda / Mamba (recommended for environment management)
- STAR
- samtools
- fastp
- MultiQC
- CIRI3 (download from https://github.com/gyjames/CIRI3 ; this workflow uses the Java 18 build)
- Python packages: `pysam`
- HOMER (`findMotifs.pl`)

All dependencies are provided via the conda environments under `workflow/rules/envs/` (including OpenJDK for CIRI3).

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
* `ciri3.jar`: CIRI3 Java archive; CIRI3 differential-expression modules are invoked via `java -jar ... DE_BSJ`, `DE_Ratio`, and `DE_Relative`
* `samples`: mapping of sample name to lists of FASTQs for R1 and R2
* `threads`: module-level CPU thread settings (configure each module independently; e.g., `threads.homer` for HOMER, `threads.samtools_view` for CIRI3 SAM conversion, `threads.splicing_stats` for splicing summary scripts)
* `output.keep_bam`: keep STAR/BWA BAM files (`false` by default to save disk)
* `deg.run_deseq2`: run optional BSJ-level DESeq2 analysis (`false` by default)
* `deg.run_de_bsj`: run CIRI3 `DE_BSJ` pairwise analysis (`false` by default)
* `deg.run_de_ratio`: run CIRI3 `DE_Ratio` pairwise analysis (`false` by default)
* `deg.run_de_relative`: run CIRI3 `DE_Relative` pairwise analysis (`false` by default)
* `deg.ciri3_gene_expression_from_featurecounts`: build CIRI3 `DE_BSJ` gene-expression input from featureCounts (`true` by default; currently required for `run_de_bsj`)
* `deg.groups`: group-to-sample mapping for DESeq2 + CIRI3 pairwise DEG analysis; each sample must exist under top-level `samples`, and enabled comparisons require at least two samples per group
* `deg.comparisons`: optional explicit pairwise comparison map; prefer `deg.groups` auto-generation unless you need a custom subset/order
* `deg.min_total_count`: minimum total BSJ count across selected samples (default `10`)
* `deg.min_samples_detected`: minimum number of samples with BSJ count > 0 (default `2`)
* `deg.padj_cutoff`: adjusted p-value threshold used for significance labels/plots (default `0.05`)
* `deg.lfc_cutoff`: absolute log2 fold-change threshold used for significance labels/plots (default `1.0`)
* `motif.enabled`: run BSJ motif module (`true` by default)
* `motif.flank`: sequence window half-size around BSJ donor/acceptor (default `30`)
* `motif.homer_len`: HOMER motif lengths for `findMotifs.pl -len` (default `"8,10,12"`)
* `motif.weight_mode`: BSJ weighting strategy (`sum`, `mean`, or `none`)

Example:

```yaml
threads:
  merge_fastq: 2
  fastp: 8
  star: 16
  bwa: 16
  samtools_index: 8
  samtools_view: 8
  featurecounts: 8
  faidx: 2
  homer: 8

reference:
  star_index: "/path/to/STAR/index"
  fasta: "/path/to/genome.fa"
  bwa_indexed_fasta: "/path/to/genome.fa"
  gtf: "/path/to/genes.gtf"

ciri3:
  jar: "/path/to/CIRI3/CIRI3_v1.0.1.jar"

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
  sampleC:
    R1:
      - "raw/sampleC_R1.fastq.gz"
    R2:
      - "raw/sampleC_R2.fastq.gz"
  sampleD:
    R1:
      - "raw/sampleD_R1.fastq.gz"
    R2:
      - "raw/sampleD_R2.fastq.gz"

deg:
  run_deseq2: true
  run_de_bsj: true
  run_de_ratio: true
  run_de_relative: true
  ciri3_gene_expression_from_featurecounts: true
  min_total_count: 10
  min_samples_detected: 2
  padj_cutoff: 0.05
  lfc_cutoff: 1.0
  groups:
    GroupA:
      - sampleA
      - sampleB
    GroupB:
      - sampleC
      - sampleD
```

Optional explicit comparisons can be provided under `deg.comparisons` when you need a custom subset or order:

```yaml
deg:
  groups:
    GroupA: [sampleA, sampleB]
    GroupB: [sampleC, sampleD]
  comparisons:
    GroupB_vs_GroupA:
      case_group: GroupB
      control_group: GroupA
```

### Migrating old CIRI3 DE config

The old top-level `ciri3_de:` block and `deg.enabled` switch are no longer supported. Update old configs to the unified `deg:` schema before running the workflow.

| Old field | New field |
| --- | --- |
| `deg.enabled` | `deg.run_deseq2` |
| `ciri3_de.enabled + ciri3_de.run_de_bsj` | `deg.run_de_bsj` |
| `ciri3_de.enabled + ciri3_de.run_de_ratio` | `deg.run_de_ratio` |
| `ciri3_de.enabled + ciri3_de.run_de_relative` | `deg.run_de_relative` |
| `ciri3_de.gene_expression_from_featurecounts` | `deg.ciri3_gene_expression_from_featurecounts` |
| `ciri3_de.outdir` | optional `deg.ciri3_outdir`; default remains `results/deg/ciri3` |
| `ciri3_de.comparisons` | optional `deg.comparisons`, preferably using `case_group`/`control_group`; omitting this lets the workflow generate comparisons from `deg.groups` |

Prefer omitting `deg.comparisons` and letting the workflow generate pairwise comparisons from `deg.groups`. If explicit sample-list comparisons are used, each case/control list must exactly match one configured `deg.groups` entry. Unknown names such as `sample_A` when the configured sample is `sampleA` fail early with a list of valid sample names.

Notes:

* The workflow assumes paired-end reads and requires both R1 and R2 lists to be the same length per sample.
* STAR genome index construction is **not** performed in this workflow; build the STAR index in advance and point `reference.star_index` to that directory.
* BWA index construction is **not** performed in this workflow; provide prebuilt index files for `reference.bwa_indexed_fasta`.
  Example:
  `bwa index /path/to/genome.fa`
* CIRI3 built-in differential-expression modules are called through the Java jar. The workflow does not call `scripts/BSJ_yes.R`, `rMATSexe`, or `lib/rmats_deps` directly.

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
