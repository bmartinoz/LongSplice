# LongSplice

A comprehensive Nextflow pipeline for long-read RNA-seq analysis with Oxford Nanopore Technologies (ONT) and PacBio sequencing data.

## Overview

LongSplice is a bioinformatic workflow developed in Nextflow that automates the analysis of long-read RNA-seq data. It integrates modern tools for quality control, preprocessing, alignment, isoform discovery, quantification, and differential analysis (expression, transcript usage, and splicing).

## ✨ Key Features

- **Input Validation**: Automatic verification of samplesheet.csv and reference files
- **Quality Control & Preprocessing**:
  - Initial and final QC with NanoPlot
  - cDNA strand orientation with PyChopper (optional)
  - Read filtering and trimming with Chopper
  - Read reorientation with Restrander

- **Alignment**:
  - Splice-aware alignment with Minimap2
  - Conversion, sorting, and indexing with Samtools
  - Coverage file generation (BigWig) and annotation (BigBed)

- **Quantification & Discovery**:
  - **Bambu**: Novel isoform discovery and gene/transcript quantification
  - **FLAIR**: Full-Length Alternative Isoform Analysis of RNA

- **Differential Analysis**:
  - **DESeq2**: Gene and transcript-level differential expression analysis
  - **DEXSeq**: Differential transcript usage (DTU) analysis
  - **SUPPA2**: Differential splicing event analysis (PSI)

- **Reporting**:
  - Alignment statistics (Flagstat, Idxstats, Stats)
  - Consolidated quality report with MultiQC

## 📋 Requirements

- **Nextflow** (>=21.10.6)
- **Container Manager** (choose one):
  - Docker (recommended)
  - Singularity
  - Conda (available, but Docker is preferred for reproducibility)

## 📂 Data Preparation (Input)

The pipeline requires a sample sheet (`samplesheet.csv`) with comma-separated values and the following mandatory columns:

| Column | Description |
|--------|-------------|
| `group` | Experimental condition (e.g., WT, KO). Used for contrasts. |
| `replicate` | Biological replicate identifier (e.g., 1, 2, 3). |
| `input_file` | Absolute path to FASTQ file (can be compressed .gz). |
| `fasta` | Absolute path to reference genome (FASTA). |
| `gtf` | Absolute path to reference annotation (GTF). |
| `technology` | Sequencing technology: ONT or PacBio. |

### Example samplesheet.csv

```csv
group,replicate,input_file,fasta,gtf,technology
WT,1,/data/wt_rep1.fastq.gz,/ref/genome.fa,/ref/genes.gtf,ONT
WT,2,/data/wt_rep2.fastq.gz,/ref/genome.fa,/ref/genes.gtf,ONT
KO,1,/data/ko_rep1.fastq.gz,/ref/genome.fa,/ref/genes.gtf,ONT
KO,2,/data/ko_rep2.fastq.gz,/ref/genome.fa,/ref/genes.gtf,ONT
```

## ⚙️ Pipeline Parameters

Configure the pipeline execution using command-line flags (`--parameter value`) or by modifying `nextflow.config`.

### Input & Output

| Parameter | Default | Description |
|-----------|---------|-------------|
| `samplesheet` | `./samplesheet1.csv` | Path to CSV file with samples. |
| `outdir` | `results` | Output directory for results. |
| `publish_dir_mode` | `copy` | File publication mode (copy, symlink). |

### Workflow Control (Skip Steps)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `skip_preprocess_ont` | `false` | Skip preprocessing (NanoPlot, Chopper, etc.). |
| `skip_pychopper` | `true` | Skip full-length cDNA detection/orientation. |
| `skip_chopper` | `false` | Skip quality filtering and trimming. |
| `skip_restrander` | `false` | Skip read reorientation. |
| `skip_alignment` | `true` | Skip Minimap2 alignment. |
| `skip_quantification` | `true` | Skip Bambu quantification. |
| `skip_dtu` | `true` | Skip differential transcript usage analysis (DEXSeq). |
| `skip_coverage` | `true` | Skip BigWig generation. |
| `skip_bigbed` | `true` | Skip BigBed generation. |

### Analysis Activation

| Parameter | Default | Description |
|-----------|---------|-------------|
| `run_nanoplot` | `false` | Run NanoPlot quality reports. |
| `run_dexseq` | `false` | Run DEXSeq for DTU. |
| `run_deseq2` | `false` | Run DESeq2 for differential expression. |
| `run_suppa` | `false` | Run SUPPA2 for splicing events. |
| `run_suppa_diff` | `false` | Run differential splicing analysis with SUPPA2. |
| `run_flair` | `true` | Run complete FLAIR pipeline. |
| `run_flair_diffexp` | `true` | Run differential expression with FLAIR. |
| `run_flair_diffsplice` | `true` | Run differential splicing with FLAIR. |

### Preprocessing Tool Options

| Parameter | Default | Description |
|-----------|---------|-------------|
| `nanoplot_opts` | `--minlength 0 --minqual 0` | Extra arguments for NanoPlot. |
| `chopper_args` | `--quality 10 --minlength 200` | Arguments for Chopper (quality/length filter). |
| `cdna_kit` | `PCS109` | cDNA kit used (for PyChopper). |
| `minimap2_opts` | `-k14` | Extra options for Minimap2 (spliced alignment). |

### Quantification & Differential Expression (Bambu/DESeq2)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `bambu_discovery` | `false` | Enable novel isoform discovery in Bambu. |
| `bambu_ndr` | `0.2` | Novel Discovery Rate. |
| `deseq2_level` | `gene` | DESeq2 analysis level (gene or transcript). |
| `deseq2_contrast_variable` | `group` | Samplesheet column for contrasts. |
| `deseq2_reference_level` | `H1975` | Control/reference level for DE. |
| `deseq2_target_level` | `HCC827` | Treatment/test level for DE. |
| `alpha` | `0.05` | Significance threshold (adjusted P-value). |
| `lfc_threshold` | `1` | Log2 Fold Change threshold. |

### FLAIR Options

| Parameter | Default | Description |
|-----------|---------|-------------|
| `flair_transcriptome_stringent` | `true` | Require full-length reads (spanning >25bp of first/last exon). |
| `flair_transcriptome_check_splice` | `true` | Verify coverage of 4/6 bp around splice sites. |
| `flair_transcriptome_support` | `3` | Minimum reads to support an isoform. |
| `shortread_junctions` | `null` | Path to BED file of short-read junctions (recommended for hybrid). |

### SUPPA2 Options

| Parameter | Default | Description |
|-----------|---------|-------------|
| `suppa_method` | `empirical` | Statistical method for SUPPA (empirical or classical). |
| `diffsplice_method` | `empirical` | Method for diffSplice. |
| `generateevents_boundary_mode` | `S` | Boundary mode for event generation (S: Strict, V: Variable). |

## 🚀 Usage

### Basic execution with Docker

```bash
nextflow run main.nf \
    --samplesheet ./samplesheet.csv \
    --outdir ./resultados \
    -profile docker
```

### Complete execution (all analyses enabled)

```bash
nextflow run main.nf \
    --samplesheet ./samplesheet.csv \
    --outdir ./resultados \
    -profile docker \
    --run_deseq2 true \
    --run_dexseq true \
    --run_suppa true \
    --run_flair true
```

### Using Singularity

```bash
nextflow run main.nf \
    --samplesheet ./samplesheet.csv \
    --outdir ./resultados \
    -profile singularity
```

## 🐳 Docker Container

The pipeline uses a main Docker container defined in `Dockerfile` that includes all necessary dependencies (Minimap2, Samtools, FLAIR, Bambu, R packages, etc.). The image is automatically built when using `-profile docker` or you can use the pre-built image `longsplice:latest`.

## 📊 Output Structure

Results are organized in the output directory with the following structure:

```
resultados/
├── qc/                    # Quality control reports
├── aligned/               # BAM files and indices
├── quantification/        # Abundance matrices
├── differential_expression/  # DESeq2/FLAIR results
├── differential_splicing/ # DEXSeq/SUPPA2 results
├── coverage/              # BigWig files
└── multiqc_report.html    # Consolidated QC report
```


## 👤 Author

- **Benjamín Matín Albornoz**

## 🤝 Contributing

Contributions are welcome! Please feel free to submit a Pull Request.Removed author and license sections from README.

## ❓ Support

For issues, questions, or feature requests, please open an issue on the GitHub repository.Removed author and license sections from README.

## 🔗 References

- [Nextflow Documentation](https://www.nextflow.io/docs/latest/index.html)
- [Minimap2](https://github.com/lh3/minimap2)
- [FLAIR](https://github.com/BrooksLabUCSC/flair)
- [Bambu](https://github.com/GoekeLab/bambu)
- [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
- [SUPPA2](https://github.com/comprna/SUPPA)



## 🔗 References

- [Nextflow Documentation](https://www.nextflow.io/docs/latest/index.html)
- [Minimap2](https://github.com/lh3/minimap2)
- [FLAIR](https://github.com/BrooksLabUCSC/flair)
- [Bambu](https://github.com/GoekeLab/bambu)
- [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
- [SUPPA2](https://github.com/comprna/SUPPA)
