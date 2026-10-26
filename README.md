# Sc-BRCA-Nextflow: A Breast Cancer Single-Cell RNA-seq Case Study

<img width="2172" height="724" alt="Untitled-1 copy" src="https://github.com/user-attachments/assets/3c934be1-3c60-4371-886a-567a7f49b4a4" />

An automated and reproducible Nextflow workflow for processing 10x Genomics single-cell RNA-seq data from human breast cancer metastasis and mammary-gland samples. The pipeline coordinates FASTQ validation, read-level quality control, Cell Ranger processing, and per-sample exploratory analysis using Seurat and Scanpy, with support for Linux HPC environments using SLURM.

## Research Team

**Horizon Science Communication LLC — Bioinformatics Team**

- **Dr. Hesham Abdullah**
- **Mohamed Elsisi**

## Project Background

Understanding the cellular heterogeneity of metastatic breast cancer is an important challenge in cancer biology. The GEO study associated with this project investigates the cellular composition of human breast cancer metastases and the molecular programs associated with aggressive metastatic cell states.

The study describes metastatic breast cancer cells as displaying a combination of gene-expression features associated with embryonic development of tubular organs and transcriptional characteristics of their original adult mammary luminal lineage. This provides the biological context for exploring gene-expression patterns at single-cell resolution.

**Study accession:** [GEO: GSE281490](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE281490)

Additional records: [BioProject PRJNA1184144](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1184144) | [SRA SRP544221](https://www.ncbi.nlm.nih.gov/sra/?term=SRP544221)

## Project Objective

The objective of this project is to establish a reproducible computational workflow for processing and exploring single-cell RNA-seq data from breast cancer metastasis and mammary-gland samples.

The workflow is designed to:

- Validate paired-end FASTQ inputs and organize lane-specific reads by biological sample.
- Assess raw sequencing quality using FastQC and MultiQC.
- Process each sample with Cell Ranger to generate gene-expression matrices and quality metrics.
- Perform per-sample exploratory analysis using Seurat and Scanpy.
- Generate reports and visualizations to support quality assessment and downstream investigation.

## Dataset Overview

| Attribute | Description |
|---|---|
| Biological samples processed | 6 |
| Paired-end FASTQ files | 48 |
| Data type | 10x Genomics single-cell RNA-seq |
| Read QC | FastQC, MultiQC |
| Quantification and matrix generation | Cell Ranger |
| Downstream exploratory analysis | Seurat, Scanpy |
| Workflow engine | Nextflow |
| Compute environment | Linux HPC / SLURM |

The dataset includes breast cancer metastatic samples and human mammary-gland samples. The original study provides the biological context; this repository documents the computational processing and exploratory analysis performed on the six selected samples.

Raw FASTQ files are not stored in this repository.

> **Results scope:** The current Seurat and Scanpy scripts perform
> exploratory analysis separately for each sample. They do not yet
> perform joint multi-sample integration, definitive cell-type
> annotation, or differential-expression analysis. The results
> highlighted here therefore focus on workflow execution, QC, and
> generated analysis outputs rather than unsupported biological
> discovery claims.

## Why Single-Cell RNA-seq and Nextflow?

### Single-cell RNA sequencing

Single-cell RNA-seq enables gene-expression profiles to be examined at the individual-cell level rather than only as an aggregate across a tissue sample. This makes it suitable for exploratory investigation of cellular heterogeneity in metastatic breast cancer and comparison of expression patterns across the selected samples.

### Nextflow workflow automation

Processing multiple samples involves repeated quality-control, quantification, and downstream-analysis steps. Nextflow organizes these steps into a modular workflow with defined inputs, outputs, and execution dependencies.

In this project, Nextflow supports:

- **Workflow orchestration:** coordinates FASTQ validation, FastQC, MultiQC, Cell Ranger, and downstream analysis.
- **Multi-sample processing:** applies the workflow to the six samples and handles lane-specific paired FASTQ inputs.
- **Reproducibility:** keeps the computational steps and execution configuration in version-controlled workflow files.
- **HPC execution:** supports execution through a SLURM-based compute environment.
- **Execution management:** supports resuming eligible completed tasks when the workflow is rerun with Nextflow's resume functionality.

Separate Conda environment manifests are provided for MultiQC, Seurat, and Scanpy. FastQC and Cell Ranger are configured through the target HPC environment and may require cluster-specific setup.

## Workflow Overview

The workflow follows this sequence:

```text
Paired-end FASTQ files
          |
          v
Input validation and lane pairing
          |
          v
     FastQC per sample
          |
          v
    Raw-read MultiQC
          |
          v
 Cell Ranger count per sample
          |
          v
   +------+------+
   |             |
   v             v
Seurat         Scanpy
analysis       analysis
   |             |
   +------+------+
          |
          v
    Final MultiQC
```
The workflow produces read-quality reports, Cell Ranger metrics and expression matrices, per-sample Seurat and Scanpy analysis outputs, and a final MultiQC report.

## Repository structure

``` text
Sc-BRCA-Nextflow/
├── main.nf
├── README.md
├── .gitattributes
├── nextflow.config
├── .gitignore
├── scripts/
│   └── submit_pipeline.sh
├── bin/
│   ├── seurat_script.R
│   └── scanpy_script.py
├── envs/
│   ├── multiqc_env.yml
│   ├── seurat_env.yml
│   └── scanpy_env.yml
└── results/
    ├── BBM0_scanpy.pdf
    ├── BBM0_seurat.pdf
    ├── BBM5_scanpy.pdf
    ├── BBM5_seurat.pdf
    ├── BSM1_scanpy.pdf
    ├── BSM1_seurat.pdf
    ├── cellranger_metrics.tsv
    ├── cellranger_summary.html
    ├── hMG1_scanpy.pdf
    ├── hMG1_seurat.pdf
    ├── hMG2_scanpy.pdf
    ├── hMG2_seurat.pdf
    ├── hMG3_scanpy.pdf
    ├── hMG3_seurat.pdf
    └── multiqc_final.html          
```

The repository includes selected analysis reports and summary metrics
under `results/`. Raw FASTQ files, the Cell Ranger reference genome,
and other large intermediate outputs are not included and must be
provided or generated separately.

## Highlights

-   **End-to-end processing:** Organizes FASTQ validation, read QC, Cell
    Ranger processing, downstream analysis, and report generation into a
    Nextflow workflow.
-   **Multi-lane handling:** Groups lane-specific paired FASTQs by
    biological sample and passes them to Cell Ranger without manually
    concatenating lanes.
-   **Two analysis ecosystems:** Generates exploratory per-sample
    outputs with both Seurat (R) and Scanpy (Python).
-   **Reproducible environments:** Defines separate Conda environments
    for MultiQC, Seurat, and Scanpy.
-   **HPC execution:** Uses a SLURM profile with process-specific
    resource settings.

## Results showcase

### 1. Read-quality assessment

The raw-read QC review reported the expected 10x read layout:

-   **R1:** 28-bp barcode/UMI read
-   **R2:** approximately 91--92-bp transcript read

Per-base quality, per-sequence quality, N-content, sequence-length, and
adapter-content checks were broadly acceptable in the reviewed report.
High duplication in R2 was interpreted in the context of amplified
single-cell libraries, rather than treated as automatic evidence that
trimming or sample rejection was necessary.

<img width="1600" height="800" alt="fastqc_sequence_counts_plot" src="https://github.com/user-attachments/assets/83344fe8-4b8b-4488-8c9d-1f3399ad0ccc" />

### 2. Cell Ranger processing

Cell Ranger was run for each of the six samples to produce per-sample
metrics, web summaries, and filtered feature-barcode matrices.

<img width="929" height="369" alt="Screenshot 2026-10-19 145332" src="https://github.com/user-attachments/assets/9d878e54-0108-49c6-a1c2-6bff4377ec27" />

### 3. Exploratory Seurat and Scanpy Analysis

Seurat and Scanpy were used to perform exploratory single-cell RNA-seq
analysis independently for each biological sample. To keep this README
concise, BBM5 is presented as a representative example.

#### Seurat — BBM5

<img width="1650" height="1275" alt="BBM5_seurat_pages-to-jpg-0004" src="https://github.com/user-attachments/assets/e1f74075-2ea0-4ecf-974c-c6e948793d51" />
<img width="1650" height="1275" alt="BBM5_seurat_pages-to-jpg-0003" src="https://github.com/user-attachments/assets/90425b25-ec99-4295-abe2-17e8b3126b38" />
<img width="1650" height="1275" alt="BBM5_seurat_pages-to-jpg-0002" src="https://github.com/user-attachments/assets/81ff7a41-3e64-4d90-a2fa-aaa2444e9a58" />
<img width="1650" height="1275" alt="BBM5_seurat_pages-to-jpg-0001" src="https://github.com/user-attachments/assets/6f819b5c-cb21-4a9a-9e5d-816422a72973" />


#### Scanpy — BBM5

<img width="2100" height="600" alt="BBM5_scanpy_pages-to-jpg-0001" src="https://github.com/user-attachments/assets/77a3fe30-08d3-4b47-8493-f0d215d335e9" />
<img width="1800" height="750" alt="BBM5_scanpy_pages-to-jpg-0002" src="https://github.com/user-attachments/assets/9d8750b4-ed2a-4253-b2f5-6163a6cc4db6" />


**Interpretation note:** BBM5 is shown as a representative sample;
the figures do not represent a joint analysis of all six samples.
Clusters with the same numeric label across separate analyses should
not automatically be interpreted as the same biological cell type.
PCA and UMAP are exploratory visualizations and do not independently
establish biological differences or causality.

### Full analysis reports

Complete PDF reports are available for all six samples. Each PDF preserves
the full set of figures generated by its respective analysis.

| Sample | Seurat report | Scanpy report |
|---|---|---|
| BBM0 | [View PDF](results/BBM0_seurat.pdf) | [View PDF](results/BBM0_scanpy.pdf) |
| BBM5 | [View PDF](results/BBM5_seurat.pdf) | [View PDF](results/BBM5_scanpy.pdf) |
| BSM1 | [View PDF](results/BSM1_seurat.pdf) | [View PDF](results/BSM1_scanpy.pdf) |
| hMG1 | [View PDF](results/hMG1_seurat.pdf) | [View PDF](results/hMG1_scanpy.pdf) |
| hMG2 | [View PDF](results/hMG2_seurat.pdf) | [View PDF](results/hMG2_scanpy.pdf) |
| hMG3 | [View PDF](results/hMG3_seurat.pdf) | [View PDF](results/hMG3_scanpy.pdf) |

Additional reports and metrics:

- [Cell Ranger metrics table](results/cellranger_metrics.tsv)
- [Cell Ranger summary](results/cellranger_summary.html)
- [Final MultiQC report](results/multiqc_final.html)

### 4. Workflow outputs

The workflow is designed to generate:

-   FastQC output and a raw-read MultiQC report
-   Cell Ranger `metrics_summary.csv` and `web_summary.html` files for
    each sample
-   Filtered feature-barcode matrices for downstream analysis
-   Per-sample Seurat `.rds` objects and PDF figures
-   Per-sample Scanpy `.h5ad` objects and PDF figures
-   A final MultiQC report

Selected reports and summary metrics are included in `results/` for
review. Large raw inputs and intermediate analysis files are excluded
from the repository.

## Dataset

The case study uses data associated with:

-   **GEO:**
    [GSE281490](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE281490)
-   **BioProject:**
    [PRJNA1184144](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1184144)
-   **SRA study:**
    [SRP544221](https://www.ncbi.nlm.nih.gov/sra/?term=SRP544221)

### Sample-to-SRA Run Mapping

The six biological samples processed were `BBM0`, `BBM5`, `BSM1`,
`hMG1`, `hMG2`, and `hMG3`.

| Sample | SRA runs | Lane range |
|---|---|---|
| hMG3 | SRR31298425 | L001 |
| hMG2 | SRR31298426 | L001 |
| hMG1 | SRR31298427 | L001 |
| BBM5 | SRR31298428–SRR31298432 | L001–L005 |
| BSM1 | SRR31298433–SRR31298440 | L001–L008 |
| BBM0 | SRR31298441–SRR31298448 | L001–L008 |

The study was reported as private at the time of the original analysis.
Confirm current data availability and applicable data-use conditions
before redistributing data or metadata. Raw FASTQ files are not stored
in this repository.

## Software

| Tool | Role | Image |
|---|---|---|
| Nextflow | Workflow orchestration | <img width="613" height="123" alt="Nextflow_logo" src="https://github.com/user-attachments/assets/a0a10264-55ef-4fe3-b54e-9fef53bb431d" /> | 
| FastQC 0.11.5 | Read-level QC | <img width="1254" height="1254" alt="FastQC Interlocking Logo" src="https://github.com/user-attachments/assets/607d6c0c-ff4b-4801-afdc-d04b32833211" /> |
| MultiQC 1.35 | QC report aggregation | <img width="1318" height="250" alt="68747470733a2f2f6769746875622e636f6d2f4d756c746951432f4d756c746951432f7261772f6d61696e2f646f63732f696d616765732f4d756c746951435f6c6f676f2e706e67" src="https://github.com/user-attachments/assets/80d4e222-614f-41f0-ad32-fc2599d2b634" /> | 
| Cell Ranger 10.0.0 | 10x Genomics processing | <img width="402" height="402" alt="9890104" src="https://github.com/user-attachments/assets/9ff6f395-07ca-4ca5-8a5b-86c8cc46f590" /> |
| Seurat 5.3.0 | R-based exploratory analysis | <img width="4167" height="4164" alt="SeuratV5" src="https://github.com/user-attachments/assets/9efa68df-5718-4418-82cf-2d17ed14be54" /> | 
| Scanpy 1.11 | Python-based exploratory analysis | <img width="1755" height="774" alt="Scanpy_Logo_RGB" src="https://github.com/user-attachments/assets/a8aaa6d0-1c8b-4bc9-a03f-eafef2dc34b7" /> |


The repository contains Conda environment manifests in `envs/`. FastQC
and Cell Ranger are configured as HPC modules; module names and resource
settings may need adaptation for another cluster. Cell Ranger uses the
`refdata-gex-GRCh38-2024-A` reference, which must be installed
separately and accessible to the compute nodes.

## Running the workflow

The workflow expects pre-downloaded paired-end FASTQ files named using
this pattern:

``` text
SAMPLE_S1_L###_R1_001.fastq.gz
SAMPLE_S1_L###_R2_001.fastq.gz
```

Place them in `results/fastqs/`. The complete case-study dataset
contains 48 FASTQ files.

Check the input count:

``` bash
find results/fastqs -maxdepth 1 -type f -name '*.fastq.gz' | wc -l
```

Expected count for the complete dataset:

``` text
48
```

### Configure the environment

Configure Conda initialization, the Cell Ranger reference, and the
required HPC modules for your system. For example:

``` bash
export CONDA_SH=/path/to/conda.sh
export TRANSCRIPTOME=/path/to/refdata-gex-GRCh38-2024-A
```

Replace example paths with real locations accessible from the compute
nodes. Confirm that the Nextflow version, Java, FastQC, Cell Ranger, and
SLURM configuration are compatible with the target HPC.

### Run raw-read QC

From the repository root:

``` bash
sbatch scripts/submit_pipeline.sh
```

The submission script defaults to QC mode. Review:

``` text
results/multiqc/raw/multiqc_raw.html
results/fastqc/
```

### Run the complete workflow

After reviewing raw-read QC:

``` bash
PIPELINE_MODE=all sbatch scripts/submit_pipeline.sh
```

The workflow uses Nextflow resume behavior to reuse eligible completed
tasks. Actual HPC setup and resource requirements depend on the target
cluster.

## Output files

| Output | Location |
|---|---|
| FastQC reports | `results/fastqc/` |
| Raw-read MultiQC | `results/multiqc/raw/multiqc_raw.html` |
| Cell Ranger metrics | `results/cellranger/<sample>/metrics_summary.csv` |
| Cell Ranger web summary | `results/cellranger/<sample>/web_summary.html` |
| Filtered feature-barcode matrix | `results/cellranger/<sample>/filtered_feature_bc_matrix/` |
| Seurat PDF and object | `results/seurat/` |
| Scanpy PDF and object | `results/scanpy/` |
| Final MultiQC report | `results/multiqc/final/multiqc_final.html` |

The following paths describe the expected outputs of a full workflow run.
The repository's `results/` directory contains the selected reports
linked in the Results showcase section.

---

## Contributors and Affiliation

This project was conducted by **Mohamed Elsisi** and **Dr. Hesham Abdullah** and belongs to **Horizon Science Communication LLC**.

### Contributors

**Dr. Hesham Abdullah**
- Email: `ceo@horizon-sci-comm.us`
- GitHub: [GitHub Profile](https://github.com/HABIOTECH)
- LinkedIn: [Dr. Hesham Abdullah](https://www.linkedin.com/in/hesham-abdullah-40bab1190/)

  
**Mohamed Elsisi**
- Email: `M.Elsisi@horizon-sci-comm.us`
- GitHub: [MohamedElsisii](https://github.com/MohamedElsisii)
- LinkedIn: [Mohamed Elsisi](https://www.linkedin.com/in/mohamed-elsisii)

### Organization

**Horizon Science Communication LLC**

<img width="1448" height="1086" alt="Screenshot 2026-06-21 211033" src="https://github.com/user-attachments/assets/eff03aac-180b-46c2-baac-8433ebdafdad" />

[Website](https://www.horizon-sci-comm.us/) 

This work is affiliated with Horizon Science Communication LLC. The project repository documents the workflow, computational methods, and exploratory analysis outputs developed for this case study.

---

For questions about the workflow or collaboration, please contact the contributors using the contact information above.
