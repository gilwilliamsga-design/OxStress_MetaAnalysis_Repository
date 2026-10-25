# Cross-Disease Oxidative Stress Biomarker Discovery in Neurodegenerative Diseases

## A reproducible transcriptomic meta-analysis workflow for identifying conserved oxidative stress biomarkers across Alzheimer's disease, Parkinson's disease, and Amyotrophic Lateral Sclerosis

---

## Overview

Oxidative stress is a central mechanism implicated in the pathogenesis of Alzheimer's disease (AD), Parkinson's disease (PD), and Amyotrophic Lateral Sclerosis (ALS). Although these disorders differ in their clinical manifestations, they share molecular alterations associated with redox imbalance, mitochondrial dysfunction, impaired antioxidant defenses, and progressive neuronal loss. Identifying oxidative stress-related genes that are consistently dysregulated across these diseases may improve our understanding of shared pathogenic mechanisms and facilitate the discovery of broadly applicable biomarkers.

This repository contains the complete R-based computational workflow developed for the transcriptomic meta-analysis accompanying the manuscript **"Cross-Disease Oxidative Stress Biomarker Discovery in Neurodegenerative Diseases."** The pipeline integrates publicly available transcriptomic datasets from the NCBI Gene Expression Omnibus (GEO), performs differential expression analyses within individual studies, combines evidence through meta-analysis, and identifies oxidative stress-related genes that are consistently altered across AD, PD, and ALS.

**Reproducibility status:** The full pipeline has completed end-to-end in the current R environment, and the regenerated outputs passed basic consistency checks. A clean-environment rerun and independent verification of all sample classifications and biological assumptions remain outstanding. A successful run alone does not establish that the comparisons or meta-analysis assumptions are valid.

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

Open the repository root in R or RStudio. Install the required packages if they are not already available:

```r
install.packages(c("dplyr", "tibble", "tidyr", "ggplot2", "VennDiagram"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("GEOquery", "Biobase", "limma", "edgeR"))
```

Then run the pipeline:

```r
source("Scripts/00_run_all.R")
```

**Validation note:** The pipeline has run end-to-end in the current environment, but has not yet been independently reproduced in a clean environment. Sample classifications, expression transformations, and cross-study comparability still require careful verification. The meta-analysis combines signed Z-scores using the unweighted Stouffer method; `Meta_logFC` is the arithmetic mean of study-level log-fold changes, not a pooled effect estimate. The scripts include checks intended to catch some invalid comparisons.
**GSE76220 sample-identity discrepancy:** During quality control, we identified a discrepancy between sample identifiers in the supplementary RPKM expression matrix and sample titles in the corresponding GEO metadata. The supplementary matrix contains a sample labelled `42c`, whereas the GEO metadata lists `85_sALS`, preventing unambiguous reconciliation of all sample identities. Because reliable sample-to-group assignment is essential for differential-expression analysis, GSE76220 is retained in the repository but excluded from the primary ALS meta-analysis pending verification. Its results should not be interpreted as part of the primary findings until this discrepancy is resolved.
---

## Outputs

The workflow generates:

* Differential expression results
* Transcriptomic meta-analysis results
* Prioritized cross-disease oxidative stress biomarkers
* Heatmaps
* Forest plots
* Publication-ready figures
* Publication-ready summary tables

## Preliminary Results

The current exploratory analysis combines two ALS studies, GSE67196 and the ALS subset of GSE26927, after excluding GSE76220 from the primary ALS meta-analysis because of an unresolved discrepancy between supplementary sample identifiers and GEO metadata.

Across the 18 oxidative stress-related genes included in the ALS analysis, **PRDX2 was the only gene with a false-discovery-rate-adjusted p-value below 0.05** (combined \(Z = -3.85\), \(p = 1.16 \times 10^{-4}\), FDR = 0.00209). Its estimated effects were negative in both contributing studies.

The current cross-disease integration also lists PRDX2 as significant in the AD, ALS, and PD meta-analysis outputs. This represents a preliminary, hypothesis-generating signal rather than validation of a diagnostic or prognostic biomarker. The analyses remain subject to dataset-specific sample classification, study comparability, and independent reproducibility checks.

These results may change as data provenance issues are resolved and the analysis workflow undergoes further validation.

---

## Citation

If you use this repository or adapt the analysis workflow, please cite the associated manuscript. Citation information will be updated following publication.

---

## Author

**Gil Williams Amana**


GitHub: https://github.com/gilwilliamsga-design
