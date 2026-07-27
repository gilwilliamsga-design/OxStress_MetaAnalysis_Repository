# Cross-Disease Oxidative Stress Biomarker Discovery in Neurodegenerative Diseases

## A reproducible transcriptomic meta-analysis workflow for identifying conserved oxidative stress biomarkers across Alzheimer's disease, Parkinson's disease, and Amyotrophic Lateral Sclerosis

---

## Overview

Oxidative stress is a central mechanism implicated in the pathogenesis of Alzheimer's disease (AD), Parkinson's disease (PD), and Amyotrophic Lateral Sclerosis (ALS). Although these disorders differ in their clinical manifestations, they share molecular alterations associated with redox imbalance, mitochondrial dysfunction, impaired antioxidant defenses, and progressive neuronal loss. Identifying oxidative stress-related genes that are consistently dysregulated across these diseases may improve our understanding of shared pathogenic mechanisms and facilitate the discovery of broadly applicable biomarkers.

This repository contains the complete R-based computational workflow developed for the transcriptomic meta-analysis accompanying the manuscript **"Cross-Disease Oxidative Stress Biomarker Discovery in Neurodegenerative Diseases."** The pipeline integrates publicly available transcriptomic datasets from the NCBI Gene Expression Omnibus (GEO), performs differential expression analyses within individual studies, combines evidence through meta-analysis, and identifies oxidative stress-related genes that are consistently altered across AD, PD, and ALS.

All analyses, figures, and tables reported in the manuscript can be reproduced using the scripts provided in this repository.

---

## Study Summary

| Feature                  | Description                                                                                         |
| ------------------------ | --------------------------------------------------------------------------------------------------- |
| **Diseases**             | Alzheimer's disease (AD), Parkinson's disease (PD), Amyotrophic Lateral Sclerosis (ALS)             |
| **Data source**          | NCBI Gene Expression Omnibus (GEO)                                                                  |
| **Datasets analyzed**    | 9 disease-specific transcriptomic datasets                                                          |
| **Programming language** | R (v4.6 or later)                                                                                   |
| **Primary analyses**     | Differential expression analysis, transcriptomic meta-analysis, cross-disease biomarker integration |
| **Primary outcome**      | Identification of conserved oxidative stress biomarkers across neurodegenerative diseases           |

---

## Repository Structure

```text
OxStress_MetaAnalysis_Repository/
│
├── Data/          # Raw, processed, and gene list files
├── Functions/     # Reusable R functions
├── Scripts/       # Complete analysis pipeline
├── Results/       # Generated figures, tables, and outputs
├── Logs/          # Analysis logs
├── .gitignore
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

The transcriptomic meta-analysis was performed using **nine disease-specific datasets** obtained from the NCBI Gene Expression Omnibus (GEO).

| Disease                       | GEO Dataset           |
| ----------------------------- | --------------------- |
| Alzheimer's disease           | GSE5281               |
| Alzheimer's disease           | GSE63061              |
| Alzheimer's disease           | GSE26927 (AD subset)  |
| Parkinson's disease           | GSE6613               |
| Parkinson's disease           | GSE8397               |
| Parkinson's disease           | GSE26927 (PD subset)  |
| Amyotrophic Lateral Sclerosis | GSE67196              |
| Amyotrophic Lateral Sclerosis | GSE76220              |
| Amyotrophic Lateral Sclerosis | GSE26927 (ALS subset) |

---

## Reproducing the Analysis

Clone the repository:

```bash
git clone https://github.com/gilwilliamsga-design/OxStress_MetaAnalysis_Repository.git
```

Open the project in R or RStudio and execute:

```r
source("Scripts/00_config.R")
source("Scripts/00_run_all.R")
```

Executing the pipeline reproduces the processed results, publication-ready figures, and summary tables generated for the accompanying manuscript.

---

## Outputs

The workflow generates:

* Differential expression results
* Transcriptomic meta-analysis results
* Prioritized cross-disease oxidative stress biomarkers
* Heatmaps
* Volcano plots
* Forest plots
* Publication-ready figures
* Publication-ready summary tables

---

## Citation

If you use this repository or adapt the analysis workflow, please cite the associated manuscript. Citation information will be updated following publication.

---

## Author

**Gil Williams Amana**


GitHub: https://github.com/gilwilliamsga-design
