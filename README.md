# Karkinos HPV Genotyping Pipeline

This repository provides a reproducible, Nextflow-based pipeline for human papillomavirus (HPV) genotyping from next-generation sequencing (NGS) data. The workflow performs read quality control, alignment to an HPV reference genome, mapping-quality filtering, and generation of alignment and coverage summaries suitable for downstream analysis and reporting.[web:193]

Although the pipeline was developed for HPV genotyping in a diagnostic and research context, the same framework can be adapted to other viral or microbial targets by changing the reference genome and input configuration.

---

## Purpose

The primary goal of this pipeline is to transform raw FASTQ data into a set of standardized, reproducible outputs that summarize:

- Read quality and trimming statistics  
- Alignment performance against an HPV reference genome  
- Per-region mapping metrics (idxstats)  
- Coverage statistics and basic genotype-level summaries  

This makes the pipeline useful for:

- Evaluating HPV NGS data in research and diagnostic feasibility studies  
- Generating consistent alignment and coverage reports for individual samples  
- Comparing samples across cohorts, time points, or sequencing runs  
- Providing a transparent analysis chain for publications or internal validation work  

---

## Scientific background

HPV genotyping is central to virological research, vaccine surveillance, and diagnostic assay evaluation. Many commercial and laboratory-developed tests are validated on a limited number of reference sequences, whereas real-world infections may involve diverse lineages and sublineages that can affect detection performance.[web:188]

NGS-based approaches provide a richer view of viral diversity but require a robust computational pipeline to process raw reads in a consistent way. Common steps include:

1. **Quality control and trimming** to remove low-quality bases and adapters (e.g. using fastp).  
2. **Alignment** of reads to a reference genome (e.g. BWA-MEM).  
3. **Filtering** based on mapping quality to reduce noise from ambiguous alignments.  
4. **Computation of statistics** such as read counts per reference sequence, coverage depth, and basic alignment metrics (SAMtools flagstat, idxstats, coverage).[web:193]

This pipeline implements these standard steps in a controlled, reproducible environment using Nextflow, and is intended to be a foundation for more advanced downstream analyses such as variant calling or lineage assignment.

---

## Quick start

This repository is designed so that a user can clone the project, install a small set of dependencies, and run the HPV genotyping pipeline on their own FASTQ data with minimal configuration.[web:193]

### 1. Clone the repository

```bash
git clone https://github.com/Sudu-09-Nep/hpv_genotype_detection.git
cd hpv_genotype_detection
```

### 2. Install dependencies

You can either use a dedicated conda environment or system-wide tools.

**Option A – conda (recommended)**

```bash
conda create -n hpv_ngs fastp bwa samtools
conda activate hpv_ngs

# Install Nextflow
curl -s https://get.nextflow.io | bash
chmod +x nextflow
mv nextflow ~/bin/  # or any directory on your PATH
```

**Option B – system packages (example)**

```bash
# Ubuntu / Debian (alignment + BAM tools)
sudo apt install bwa samtools

# macOS (Homebrew)
brew install fastp bwa samtools nextflow
```

Confirm that the following commands are available on your `PATH`:

```bash
fastp -h
bwa
samtools
nextflow -version
```

### 3. Prepare input data and configuration

By default, the repository assumes the following layout (as shipped):

```text
data/
├── raw_fastq/
│   ├── 7936_S200_L008_R1_001.fastq
│   └── 7936_S200_L008_R2_001.fastq
├── reference/
│   ├── reference.fasta
│   ├── reference.fasta.amb
│   ├── reference.fasta.ann
│   ├── reference.fasta.bwt
│   ├── reference.fasta.pac
│   └── reference.fasta.sa
└── samplesheet.csv
```

1. **Raw FASTQ files**

   Place your own FASTQ files in `data/raw_fastq/`. For example:

   ```text
   data/raw_fastq/
   ├── sample1_L001_R1.fastq.gz
   └── sample1_L001_R2.fastq.gz
   ```

2. **Reference genome**

   The reference genome is stored in `data/reference/reference.fasta` and indexed for BWA. To use a different HPV type or another organism, replace this FASTA file (and regenerate the BWA index if necessary).

3. **Sample sheet**

   The file `data/samplesheet.csv` describes the samples and input files. A minimal example:

   ```text
   sample_id,mode,fastq_r1,fastq_r2
   7936_S200,PE,data/raw_fastq/7936_S200_L008_R1_001.fastq,data/raw_fastq/7936_S200_L008_R2_001.fastq
   ```

   - `sample_id` is an arbitrary but unique sample identifier.  
   - `mode` can be `PE` (paired-end) or `SE` (single-end).  
   - `fastq_r1` and `fastq_r2` point to the read files; for single-end data, `fastq_r2` should be left empty.

### 4. Run the pipeline

From the project root:

```bash
nextflow run workflows/hpv_ngs.nf
```

The default configuration uses `nextflow.config` and the paths under `data/` as shown above. You can adjust resources (e.g. CPU threads) and other parameters in `nextflow.config` or via command-line parameters.

### 5. Inspect outputs

After a successful run, the main outputs are organised under `data/`:

```text
data/
├── trimmed/
│   ├── 7936_S200_fastp.html
│   ├── 7936_S200_fastp.json
│   ├── 7936_S200_R1_paired.fastq
│   ├── 7936_S200_R1_unpaired.fastq
│   ├── 7936_S200_R2_paired.fastq
│   └── 7936_S200_R2_unpaired.fastq
├── alignment/
│   ├── 7936_S200.sorted.bam
│   ├── 7936_S200.sorted.bam.bai
│   ├── 7936_S200.filtered_MAPQ20.bam
│   ├── 7936_S200.filtered_MAPQ20.bam.bai
│   └── 7936_S200_idxstats.txt
├── coverage/
│   └── 7936_S200_coverage.txt
├── logs/
│   └── 7936_S200_flagstat.txt
└── reports/
    └── 7936_S200_genotype_report.txt
```

These files contain:

- fastp QC reports (`trimmed/`)  
- BAM alignments and indices (`alignment/`)  
- Per-reference read counts (`*_idxstats.txt`)  
- Coverage summaries (`coverage/`)  
- Alignment statistics (`logs/`)  
- A consolidated genotype/alignment report (`reports/`)  

---

## Workflow overview

Conceptually, the pipeline is structured as follows:

1. **Input definition** via `data/samplesheet.csv`.  
2. **Read quality control and trimming** using fastp, generating paired/unpaired FASTQ files and HTML/JSON QC reports.  
3. **Alignment** of trimmed reads to `data/reference/reference.fasta` using BWA-MEM.  
4. **Sorting and indexing** of BAM files.  
5. **Mapping-quality filtering** (e.g. MAPQ ≥ 20) to produce `*.filtered_MAPQ20.bam`.  
6. **Summary statistics** via SAMtools (flagstat, idxstats, coverage).  
7. **Report generation**, combining key metrics into `*_genotype_report.txt` per sample.[web:193]

The entire workflow is orchestrated with Nextflow and implemented primarily in `workflows/hpv_ngs.nf`, with helper logic in `bin/run_hpv_pipeline.sh`.

---

## Repository structure

The repository is organised as:

```text
hpv_ngs_pipeline/
├── bin
│   └── run_hpv_pipeline.sh
├── data
│   ├── alignment/
│   ├── coverage/
│   ├── fastqc/
│   ├── logs/
│   ├── raw_fastq/
│   ├── reference/
│   ├── reports/
│   ├── samplesheet.csv
│   └── trimmed/
├── nextflow.config
├── README.md
└── workflows
    └── hpv_ngs.nf
```

- `bin/run_hpv_pipeline.sh` – helper script called within the Nextflow workflow (if applicable).  
- `workflows/hpv_ngs.nf` – main Nextflow script describing the pipeline.  
- `nextflow.config` – configuration for resources, parameters, and default paths.  
- `data/` – input data and generated outputs (see Quick start section for details).

Large data products such as FASTQ and BAM files are generated by the pipeline and can be excluded from version control (see `.gitignore`) to keep the repository lightweight and portable.

---

## Input requirements

The pipeline expects three main inputs:

1. **Raw sequencing reads**

   - Single-end or paired-end FASTQ files (optionally compressed).  
   - Stored under `data/raw_fastq/` and referenced in `data/samplesheet.csv`.

2. **Reference genome**

   - A FASTA file at `data/reference/reference.fasta` representing the HPV genome (or another target).  
   - BWA index files (`*.amb`, `*.ann`, `*.bwt`, `*.pac`, `*.sa`) can be pre-generated or created automatically by the pipeline.

3. **Sample configuration**

   - A CSV file `data/samplesheet.csv` linking sample IDs to FASTQ file paths and specifying `PE` or `SE` mode.

This design makes the workflow generic: changing the target organism or assay involves replacing the reference FASTA and updating the samplesheet, without modifying the core pipeline.

---

## Output description

For each sample, the pipeline typically generates:

- **fastp outputs**  
  - `*_fastp.html` – interactive HTML QC report.  
  - `*_fastp.json` – machine-readable QC report.

- **Alignment outputs**  
  - `*.sorted.bam` – sorted aligned BAM file.  
  - `*.sorted.bam.bai` – BAM index.  
  - `*.filtered_MAPQ20.bam` – BAM after MAPQ-based filtering.  
  - `*.filtered_MAPQ20.bam.bai` – index for the filtered BAM.  
  - `*_idxstats.txt` – per-reference mapping statistics.

- **Coverage outputs**  
  - `*_coverage.txt` – SAMtools coverage summary.

- **Log and report files**  
  - `*_flagstat.txt` – SAMtools flagstat summary.  
  - `*_genotype_report.txt` – text report summarising key metrics and basic genotype-level interpretation.

These outputs can be used directly for exploratory analysis or integrated into downstream workflows (e.g. variant calling, lineage analysis, or reporting).

---

## Single-end versus paired-end support

The pipeline supports both single-end and paired-end sequencing data:

- **Single-end mode (`SE`)** is used when only one read per fragment is available. QC and alignment focus on per-read metrics.  
- **Paired-end mode (`PE`)** uses information from both mates, improving alignment accuracy and enabling additional metrics such as insert size distribution.[web:193]

The mode is specified in `data/samplesheet.csv` and handled automatically within the Nextflow workflow.

---

## Limitations

This pipeline focuses on **read processing, alignment, and descriptive statistics**. It does **not** currently implement:

- Variant calling or haplotype reconstruction  
- Lineage or sublineage assignment  
- Primer/probe-level analysis  
- Thermodynamic modelling or clinical performance prediction  

Wet-lab experiments and additional bioinformatics analyses remain necessary for full assay validation and clinical interpretation.

---

## Future directions

Potential future extensions include:

- Integration of variant calling workflows (e.g. FreeBayes, GATK)  
- Automated lineage assignment for HPV based on variant profiles  
- Standardised CSV/TSV outputs for multi-sample summarisation  
- HTML or notebook-based reports for interactive exploration  
- Containerisation (Docker/Singularity) for fully encapsulated deployments.[web:195][web:196]

---

## Requirements

### Software

- Operating system: Linux or macOS  
- Shell: `bash` (≥ 4.0)  
- Core tools:  
  - fastp (quality control and trimming)  
  - BWA (alignment)  
  - SAMtools (BAM processing and statistics)  
  - Nextflow (workflow orchestration)[web:193]  

### Hardware (typical)

- RAM: ≥ 8 GB (16 GB recommended for larger datasets)  
- CPU: ≥ 4 cores (the pipeline can exploit additional cores via Nextflow)  
- Storage: sufficient for raw FASTQ and BAM files (NGS-scale).  

---

## Installation notes

Tool installation strategies differ by environment. For many users, a conda-based setup (as described in the Quick start) provides a simple and reproducible route; others may prefer system packages or containerised deployments.[web:193]

Cluster users can adapt `nextflow.config` to their scheduler (e.g. SLURM, SGE) by defining appropriate executors and resource profiles.

---

## License

This project is distributed under the **MIT License**. See the `LICENSE` file for the full text.

---

## Citation

If you use this pipeline in research or publications, please cite:

- The HPV reference genome or dataset that you used  
- The sequencing dataset (if publicly available)  
- This repository, for example:

> Aryal S. *Karkinos HPV Genotyping Pipeline*. GitHub, 2026.  
> Repository: https://github.com/Sudu-09-Nep/hpv_ngs_pipeline

---

## Acknowledgments

This pipeline builds on widely used open-source tools (fastp, BWA, SAMtools, Nextflow) and follows best practices inspired by established workflow collections such as nf-core.[web:193][web:195]
---

## Repository structure

The repository is organised as:

```text
hpv_ngs_pipeline/
├── bin
│   └── run_hpv_pipeline.sh
├── data
│   ├── alignment/
│   ├── coverage/
│   ├── fastqc/
│   ├── logs/
│   ├── raw_fastq/
│   ├── reference/
│   ├── reports/
│   ├── samplesheet.csv
│   └── trimmed/
├── nextflow.config
├── README.md
└── workflows
    └── hpv_ngs.nf
```

- `bin/run_hpv_pipeline.sh` – helper script called within the Nextflow workflow (if applicable).  
- `workflows/hpv_ngs.nf` – main Nextflow script describing the pipeline.  
- `nextflow.config` – configuration for resources, parameters, and default paths.  
- `data/` – input data and generated outputs (see Quick start section for details).

Large data products such as FASTQ and BAM files are generated by the pipeline and can be excluded from version control (see `.gitignore`) to keep the repository lightweight and portable.

---

## Input requirements

The pipeline expects three main inputs:

1. **Raw sequencing reads**

   - Single-end or paired-end FASTQ files (optionally compressed).  
   - Stored under `data/raw_fastq/` and referenced in `data/samplesheet.csv`.

2. **Reference genome**

   - A FASTA file at `data/reference/reference.fasta` representing the HPV genome (or another target).  
   - BWA index files (`*.amb`, `*.ann`, `*.bwt`, `*.pac`, `*.sa`) can be pre-generated or created automatically by the pipeline.

3. **Sample configuration**

   - A CSV file `data/samplesheet.csv` linking sample IDs to FASTQ file paths and specifying `PE` or `SE` mode.

This design makes the workflow generic: changing the target organism or assay involves replacing the reference FASTA and updating the samplesheet, without modifying the core pipeline.

---

## Output description

For each sample, the pipeline typically generates:

- **fastp outputs**  
  - `*_fastp.html` – interactive HTML QC report.  
  - `*_fastp.json` – machine-readable QC report.

- **Alignment outputs**  
  - `*.sorted.bam` – sorted aligned BAM file.  
  - `*.sorted.bam.bai` – BAM index.  
  - `*.filtered_MAPQ20.bam` – BAM after MAPQ-based filtering.  
  - `*.filtered_MAPQ20.bam.bai` – index for the filtered BAM.  
  - `*_idxstats.txt` – per-reference mapping statistics.

- **Coverage outputs**  
  - `*_coverage.txt` – SAMtools coverage summary.

- **Log and report files**  
  - `*_flagstat.txt` – SAMtools flagstat summary.  
  - `*_genotype_report.txt` – text report summarising key metrics and basic genotype-level interpretation.

These outputs can be used directly for exploratory analysis or integrated into downstream workflows (e.g. variant calling, lineage analysis, or reporting).

---

## Single-end versus paired-end support

The pipeline supports both single-end and paired-end sequencing data:

- **Single-end mode (`SE`)** is used when only one read per fragment is available. QC and alignment focus on per-read metrics.  
- **Paired-end mode (`PE`)** uses information from both mates, improving alignment accuracy and enabling additional metrics such as insert size distribution.[web:193]

The mode is specified in `data/samplesheet.csv` and handled automatically within the Nextflow workflow.

---

## Limitations

This pipeline focuses on **read processing, alignment, and descriptive statistics**. It does **not** currently implement:

- Variant calling or haplotype reconstruction  
- Lineage or sublineage assignment  
- Primer/probe-level analysis  
- Thermodynamic modelling or clinical performance prediction  

Wet-lab experiments and additional bioinformatics analyses remain necessary for full assay validation and clinical interpretation.

---

## Future directions

Potential future extensions include:

- Integration of variant calling workflows (e.g. FreeBayes, GATK)  
- Automated lineage assignment for HPV based on variant profiles  
- Standardised CSV/TSV outputs for multi-sample summarisation  
- HTML or notebook-based reports for interactive exploration  
- Containerisation (Docker/Singularity) for fully encapsulated deployments.[web:195][web:196]

---

## Requirements

### Software

- Operating system: Linux or macOS  
- Shell: `bash` (≥ 4.0)  
- Core tools:  
  - fastp (quality control and trimming)  
  - BWA (alignment)  
  - SAMtools (BAM processing and statistics)  
  - Nextflow (workflow orchestration)[web:193]  

### Hardware (typical)

- RAM: ≥ 8 GB (16 GB recommended for larger datasets)  
- CPU: ≥ 4 cores (the pipeline can exploit additional cores via Nextflow)  
- Storage: sufficient for raw FASTQ and BAM files (NGS-scale).  

---

## Installation notes

Tool installation strategies differ by environment. For many users, a conda-based setup (as described in the Quick start) provides a simple and reproducible route; others may prefer system packages or containerised deployments.[web:193]

Cluster users can adapt `nextflow.config` to their scheduler (e.g. SLURM, SGE) by defining appropriate executors and resource profiles.

---

## License

This project is distributed under the **MIT License**. See the `LICENSE` file for the full text.

---

## Citation

If you use this pipeline in research or publications, please cite:

- The HPV reference genome or dataset that you used  
- The sequencing dataset (if publicly available)  
- This repository, for example:

> Aryal S. *Karkinos HPV Genotyping Pipeline*. GitHub, 2026.  
> Repository: https://github.com/Sudu-09-Nep/hpv_ngs_pipeline

---

## Acknowledgments

This pipeline builds on widely used open-source tools (fastp, BWA, SAMtools, Nextflow) and follows best practices inspired by established workflow collections such as nf-core.
