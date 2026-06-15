# Karkinos HPV Genotyping Pipeline

This repository provides a reproducible, Nextflow-based workflow for human papillomavirus (HPV) genotyping from next-generation sequencing (NGS) data. The pipeline performs read quality control, alignment to an HPV reference genome, mapping-quality filtering, and generation of alignment and coverage summaries suitable for downstream analysis and reporting.

The workflow was originally developed for HPV NGS data but can be adapted to other viral or microbial targets by changing the reference genome and input configuration.

---

## Purpose

The primary goal of this pipeline is to transform raw FASTQ data into a set of standardized, reproducible outputs that summarize:

- Read quality and trimming statistics  
- Alignment performance against an HPV reference genome  
- Per-region mapping metrics (idxstats)  
- Coverage statistics and genotype-level summary metrics  

This makes the pipeline useful for:

- Evaluating HPV NGS data in research and diagnostic feasibility studies  
- Generating consistent alignment and coverage reports for individual samples  
- Comparing samples across cohorts, time points, or sequencing runs  
- Providing a transparent analysis chain for publications or internal validation work  

---

## Scientific background

HPV genotyping is central to virological research, vaccine surveillance, and diagnostic assay evaluation. Many commercial and laboratory-developed tests are validated on a limited number of reference sequences, whereas real-world infections may involve diverse lineages and sublineages that can affect detection performance.

NGS-based approaches provide a richer view of viral diversity but require a robust computational pipeline to process raw reads in a consistent way. Common steps include:

1. Quality control and trimming to remove low-quality bases and adapters.  
2. Alignment of reads to a reference genome.  
3. Filtering based on mapping quality to reduce noise from ambiguous alignments.  
4. Computation of statistics such as read counts per reference sequence, coverage depth, and basic alignment metrics.

This pipeline implements these standard steps in a controlled, reproducible environment using fastp for trimming and QC, BWA-MEM for alignment, SAMtools for BAM processing and statistics, and Nextflow for workflow orchestration [1–4]. It is intended to serve as a foundation for more advanced downstream analyses such as variant calling or lineage assignment.

---

## Installation

This section describes how to obtain the code and install the minimal set of software required to run the pipeline.

### Prerequisites

The pipeline has been tested on Linux and macOS systems with the following tools:

- `bash` (≥ 4.0)  
- `fastp` (for read trimming and quality control) [1]  
- `bwa` (for alignment using the BWA-MEM algorithm) [2]  
- `samtools` (for BAM processing and summary statistics) [3]  
- `nextflow` (for workflow orchestration) [4]  

We recommend managing these dependencies with **conda**, but system packages (e.g. `apt`, `brew`) or containers can also be used.

### Clone the repository

```bash
git clone https://github.com/Sudu-09-Nep/hpv_ngs_pipeline.git
cd hpv_ngs_pipeline
```

### Install dependencies (conda, recommended)

```bash
conda create -n hpv_ngs fastp bwa samtools
conda activate hpv_ngs

# Install Nextflow
curl -s https://get.nextflow.io | bash
chmod +x nextflow
mv nextflow ~/bin/  # or any directory on your PATH
```

Alternatively, install the tools with your preferred package manager, ensuring that `fastp`, `bwa`, `samtools`, and `nextflow` are available on your `PATH`:

```bash
fastp -h
bwa
samtools
nextflow -version
```

---

## Directory layout and where to put your data

After cloning, the repository provides an empty `data/` structure (directories contain only `.gitkeep` files to preserve the layout under version control):

```text
data/
├── raw_fastq/   # Place your raw FASTQ files here
├── reference/   # Place your reference.fasta (and indices) here
├── trimmed/     # Created by the pipeline (fastp outputs)
├── alignment/   # Created by the pipeline (BAM files and idxstats)
├── coverage/    # Created by the pipeline (coverage summaries)
├── logs/        # Created by the pipeline (flagstat and logs)
├── reports/     # Created by the pipeline (genotype reports)
├── fastqc/      # Reserved for future QC outputs
└── samplesheet.csv
```

- **Input FASTQ files**: add to `data/raw_fastq/` and reference them in `data/samplesheet.csv`.  
- **Reference genome**: place as `data/reference/reference.fasta` (BWA indices are generated or reused in the same directory).  
- All other directories are filled automatically by the pipeline when you run it.

---

## Quick start

Once the repository is cloned and dependencies are installed, the pipeline can be run on any suitably formatted dataset.

### 1. Prepare input data and configuration

1. **Raw FASTQ files**

   Place your FASTQ files in `data/raw_fastq/`. For example:

   ```text
   data/raw_fastq/
   ├── 7936_S200_L008_R1_001.fastq
   └── 7936_S200_L008_R2_001.fastq
   ```

2. **Reference genome**

   Place your HPV reference genome in:

   ```text
   data/reference/reference.fasta
   ```

   BWA index files (`.amb`, `.ann`, `.bwt`, `.pac`, `.sa`) can either be pre-generated (by running `bwa index`) or generated by the pipeline if that logic is enabled.

3. **Sample sheet**

   The file `data/samplesheet.csv` describes the samples and input files. A minimal example:

   ```text
   sample_id,mode,fastq_r1,fastq_r2
   7936_S200,PE,data/raw_fastq/7936_S200_L008_R1_001.fastq,data/raw_fastq/7936_S200_L008_R2_001.fastq
   sample1,SE,data/raw_fastq/sample1_L001_R1.fastq,
   ```

   - `sample_id` – unique sample identifier  
   - `mode` – `PE` for paired-end, `SE` for single-end  
   - `fastq_r1`, `fastq_r2` – paths to the read files (leave `fastq_r2` empty for single-end)

### 2. Run the pipeline

From the project root:

```bash
nextflow run workflows/hpv_ngs.nf
```

**To Generate the Multiqc report
```bash
bin/generate_multiqc_all.sh
```

The default configuration (`nextflow.config`) uses:

- `data/samplesheet.csv` for sample metadata  
- `data/reference/reference.fasta` as the alignment reference  
- Standard fastp, BWA, and SAMtools settings tuned for HPV-scale genomes  

You can adjust resources (e.g. CPU threads, memory) and parameters in `nextflow.config` or via the Nextflow command line.

### 3. Inspect outputs

After a successful run, outputs are written under `data/` (for example, for sample `7936_S200`):

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

These include:

- fastp QC reports (`trimmed/`)  
- BAM alignments and indices (`alignment/`)  
- idxstats and coverage summaries (`alignment/`, `coverage/`)  
- flagstat logs (`logs/`)  
- genotype/alignment reports (`reports/`)  

---

## Workflow overview

Conceptually, the pipeline is structured as follows:

1. **Input definition** via `data/samplesheet.csv`.  
2. **Read quality control and trimming** using fastp, generating paired/unpaired FASTQ files and HTML/JSON QC reports [1].  
3. **Alignment** of trimmed reads to `data/reference/reference.fasta` using the BWA-MEM algorithm [2].  
4. **Sorting and indexing** of BAM files using SAMtools [3].  
5. **Mapping-quality filtering** (e.g. MAPQ ≥ 20) to produce `*.filtered_MAPQ20.bam`.  
6. **Summary statistics** via SAMtools (flagstat, idxstats, coverage) [3].  
7. **Report generation**, combining key metrics into `*_genotype_report.txt` per sample.

The entire workflow is orchestrated with Nextflow and implemented primarily in `workflows/hpv_ngs.nf`, with helper logic in `bin/run_hpv_pipeline.sh` [4].

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

- `bin/run_hpv_pipeline.sh` – helper script(s) used within the Nextflow workflow.  
- `workflows/hpv_ngs.nf` – main Nextflow script describing the pipeline.  
- `nextflow.config` – configuration for resources, parameters, and default paths.  
- `data/` – input data and generated outputs (see sections above for details).

Large data products such as FASTQ and BAM files are generated by the pipeline and excluded from version control to keep the repository lightweight and portable.

---

## Single-end versus paired-end support

The pipeline supports both single-end and paired-end sequencing data:

- **Single-end mode (`SE`)** is used when only one read per fragment is available; QC and alignment focus on per-read metrics.  
- **Paired-end mode (`PE`)** uses information from both mates, improving alignment accuracy and enabling additional metrics such as insert size distribution.

The mode is specified in `data/samplesheet.csv` and handled automatically within the Nextflow workflow.

---

## Limitations

This pipeline focuses on read processing, alignment, and descriptive statistics. It does **not** currently implement:

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
- Containerisation (Docker/Singularity) for fully encapsulated deployments  

---

## Requirements summary

**Software**

- Operating system: Linux or macOS  
- Shell: `bash` (≥ 4.0)  
- Core tools: fastp, BWA, SAMtools, Nextflow [1–4]  

**Hardware (typical)**

- RAM: ≥ 8 GB (16 GB recommended for larger datasets)  
- CPU: ≥ 4 cores (the pipeline can exploit additional cores via Nextflow)  
- Storage: sufficient for raw FASTQ and BAM files (NGS-scale)  

---

## Citation

If you use this pipeline in research or publications, please cite:

- The HPV reference genome or dataset used  
- The sequencing dataset (if publicly available)  
- This repository, for example:

> Aryal S. *Karkinos HPV Genotyping Pipeline*. GitHub, 2026.  
> Repository: https://github.com/Sudu-09-Nep/hpv_ngs_pipeline

In addition, please cite the core tools as appropriate (see References below).

---

## References

[1] Chen S, Zhou Y, Chen Y, Gu J. fastp: an ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics*. 2018;34(17):i884–i890. doi:10.1093/bioinformatics/bty560.  
[2] Li H. Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. *arXiv preprint* arXiv:1303.3997. 2013.  
[3] Li H, Handsaker B, Wysoker A, et al. The Sequence Alignment/Map format and SAMtools. *Bioinformatics*. 2009;25(16):2078–2079.  
[4] Di Tommaso P, Chatzou M, Floden EW, et al. Nextflow enables reproducible computational workflows. *Nature Biotechnology*. 2017;35(4):316–319. doi:10.1038/nbt.3820.  
