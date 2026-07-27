# Cross-Disease Oxidative Stress Biomarker Discovery in Neurodegenerative Diseases

## A reproducible transcriptomic meta-analysis workflow for identifying conserved oxidative stress biomarkers across Alzheimer's disease, Parkinson's disease, and Amyotrophic Lateral Sclerosis

---

## Overview

Neurodegenerative diseases (NDs) such as Alzheimer's disease (AD), Parkinson's disease (PD), and Amyotrophic Lateral Sclerosis (ALS) share several pathogenic mechanisms despite their distinct clinical presentations. Among these, oxidative stress is recognized as a major contributor to neuronal dysfunction, mitochondrial impairment, and progressive neurodegeneration. Identifying oxidative stress-related genes that are consistently dysregulated across multiple neurodegenerative diseases may improve our understanding of shared disease mechanisms and highlight candidate biomarkers for future diagnostic and therapeutic research.

This repository contains the complete R-based computational workflow developed for the transcriptomic meta-analysis accompanying the manuscript **"Cross-Disease Oxidative Stress Biomarker Discovery in Neurodegenerative Diseases."** The pipeline integrates publicly available transcriptomic datasets from the NCBI Gene Expression Omnibus (GEO), performs differential expression analyses within individual studies, combines evidence through meta-analysis, and identifies oxidative stress-related genes that are consistently altered across AD, PD, and ALS.

The repository emphasizes transparency and reproducibility by providing the complete analysis pipeline used to generate the processed results, publication-ready figures, and summary tables reported in the accompanying manuscript.

---

## Objectives

The primary objectives of this project are to:

* Identify oxidative stress-related differentially expressed genes in AD, PD, and ALS.
* Integrate transcriptomic evidence across multiple independent GEO studies.
* Perform cross-disease transcriptomic meta-analysis.
* Prioritize conserved oxidative stress biomarkers across neurodegenerative diseases.
* Generate publication-ready figures and summary tables.
* Provide a fully reproducible computational workflow.

---

## Repository Structure

```text
OxStress_MetaAnalysis_Repository/

├── Data/
│   ├── Raw/
│   ├── Processed/
│   └── GeneLists/
│
├── Functions/
│   ├── io_functions.R
│   ├── expression_functions.R
│   └── stats_functions.R
│
├── Scripts/
│   ├── 00_config.R
│   ├── 00_run_all.R
│   ├── 01_data_download.R
│   ├── 02_preprocessing.R
│   ├── 03_differential_expression.R
│   ├── 04_meta_analysis.R
│   ├── 05_cross_disease.R
│   └── 06_figures.R
│
├── Results/
├── Logs/
└── README.md
```

---

## Analysis Workflow

```text
GEO Datasets
      │
      ▼
Data Acquisition
      │
      ▼
Quality Control & Preprocessing
      │
      ▼
Differential Expression Analysis
      │
      ▼
Oxidative Stress Gene Extraction
      │
      ▼
Transcriptomic Meta-analysis
      │
      ▼
Cross-Disease Integration
      │
      ▼
Biomarker Prioritization
      │
      ▼
Figure & Table Generation
```

---

## Datasets

The analysis uses publicly available transcriptomic datasets obtained from the NCBI Gene Expression Omnibus (GEO).

| Disease                       | GEO Accession      |
| ----------------------------- | ------------------ |
| Alzheimer's disease           | GSE5281, GSE63061  |
| Parkinson's disease           | GSE6613, GSE8397   |
| Amyotrophic Lateral Sclerosis | GSE67196, GSE76220 |

---

## Software Requirements

The workflow was developed in **R (version 4.6 or later)** using Bioconductor and CRAN packages, including:

* GEOquery
* limma
* metafor
* meta
* tidyverse
* data.table
* ggplot2
* pheatmap

Package installation and configuration are handled through:

```r
source("Scripts/00_config.R")
```

---

## Running the Analysis

Clone the repository:

```bash
git clone https://github.com/gilwilliamsga-design/OxStress_MetaAnalysis_Repository.git
```

Open the project in R or RStudio and execute the complete workflow:

```r
source("Scripts/00_run_all.R")
```

---

## Outputs

Running the pipeline generates:

* Differential expression results
* Transcriptomic meta-analysis results
* Prioritized cross-disease oxidative stress biomarkers
* Heatmaps
* Volcano plots
* Forest plots
* Publication-ready figures
* Publication-ready summary tables

---

## Reproducibility

All analyses are performed using publicly available GEO datasets. Executing the analysis pipeline reproduces the processed results, figures, and tables presented in the accompanying manuscript.

---

## Citation

If you use or adapt this workflow, please cite the associated manuscript. Citation details will be updated following publication.

---

## Author

**Gil Williams Amana**


GitHub: https://github.com/gilwilliamsga-design
