# Cross-Disease Meta-analysis of Oxidative Stress Genes in Neurodegenerative Diseases

**Author:** Gil Williams Amana  
**Version:** 1.0.0

---

## Overview

This project performs a transcriptomic meta-analysis of oxidative stress gene expression across three major neurodegenerative diseases:

- **Alzheimer's Disease (AD)**
- **Amyotrophic Lateral Sclerosis (ALS)**
- **Parkinson's Disease (PD)**

Data from nine GEO datasets are integrated, differentially expressed oxidative stress genes are identified per disease, and cross-disease candidate biomarkers are prioritised.

---

## Repository Structure

```
OxStress_MetaAnalysis/
├── README.md
├── Scripts/
│   ├── 00_config.R                    # Global settings, paths, thresholds
│   ├── 01_data_download.R             # Download GEO datasets
│   ├── 02_preprocessing.R             # Probe annotation, gene extraction, normalisation
│   ├── 03_differential_expression.R   # limma DEA per dataset per disease
│   ├── 04_meta_analysis.R             # Stouffer meta-analysis per disease
│   ├── 05_cross_disease.R             # Cross-disease integration & biomarker ranking
│   └── 06_figures.R                   # All publication figures
├── Functions/
│   ├── io_functions.R                 # File I/O helpers
│   ├── expression_functions.R         # Extract, normalise, z-score
│   └── stats_functions.R             # meta_stouffer(), make_deg_table()
├── Data/
│   ├── Raw/                           # Downloaded GEO series matrix files
│   ├── Processed/                     # Per-dataset expression matrices (.RData)
│   └── GeneLists/
│       └── antioxidant_genes.csv      # Curated oxidative stress gene set
├── Results/
│   ├── Figures/                       # PNG, TIFF, PDF, SVG outputs
│   ├── Tables/                        # CSV result tables
│   ├── MetaAnalysis/                  # Meta-analysis RData objects
│   └── QC/                            # QC plots and logs
└── Logs/
    └── pipeline.log
```

---

## Datasets

| GEO ID   | Disease    | Samples |
|----------|------------|---------|
| GSE5281  | AD         | 161     |
| GSE63061 | AD         | 388     |
| GSE26927 | AD/ALS/PD  | 118     |
| GSE67196 | ALS        | 53      |
| GSE76220 | ALS        | 21      |
| GSE6613  | PD         | 105     |
| GSE8397  | PD         | 47      |

---

## How to Run

```r
source("Scripts/00_config.R")
source("Scripts/01_data_download.R")
source("Scripts/02_preprocessing.R")
source("Scripts/03_differential_expression.R")
source("Scripts/04_meta_analysis.R")
source("Scripts/05_cross_disease.R")
source("Scripts/06_figures.R")
```

---

## Dependencies

```r
install.packages(c("tidyverse", "pheatmap", "VennDiagram", "ggplot2", "grid", "BiocManager"))
BiocManager::install(c("GEOquery", "limma", "Biobase", "sva"))
```

---

## Key Findings

- **PRDX2** and **SOD1** were identified as shared oxidative stress biomarkers significantly dysregulated across all three neurodegenerative diseases (FDR < 0.05).
- Both genes were consistently downregulated across AD, ALS, and PD cohorts.
