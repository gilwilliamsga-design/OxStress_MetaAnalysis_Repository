###############################################################################
# 00_config.R
#
# Project : Cross-Disease Meta-analysis of Oxidative Stress Genes
#           in Alzheimer's Disease, ALS, and Parkinson's Disease
# Author  : Gil Williams Amana
# Version : 1.0.0
#
# Purpose : Central configuration file. Source this at the top of every
#           script. Defines all paths, thresholds, and global options so
#           nothing is hard-coded downstream.
###############################################################################

# ── Project information ───────────────────────────────────────────────────────
PROJECT_NAME <- "OxStress_MetaAnalysis"
AUTHOR       <- "Gil Williams Amana"
VERSION      <- "1.0.0"

# ── Directory structure ───────────────────────────────────────────────────────
ROOT_DIR      <- getwd()          # project root (where this repo lives)
# If you don't use {here}, replace with: ROOT_DIR <- getwd()

DATA_DIR      <- file.path(ROOT_DIR, "Data")
RAW_DIR       <- file.path(DATA_DIR,  "Raw")
PROCESSED_DIR <- file.path(DATA_DIR,  "Processed")
GENELIST_DIR  <- file.path(DATA_DIR,  "GeneLists")

RESULTS_DIR   <- file.path(ROOT_DIR, "Results")
FIGURE_DIR    <- file.path(RESULTS_DIR, "Figures")
TABLE_DIR     <- file.path(RESULTS_DIR, "Tables")
META_DIR      <- file.path(RESULTS_DIR, "MetaAnalysis")
QC_DIR        <- file.path(RESULTS_DIR, "QC")
LOG_DIR       <- file.path(ROOT_DIR, "Logs")

# Create all directories (safe to re-run)
for (d in c(RAW_DIR, PROCESSED_DIR, GENELIST_DIR,
            FIGURE_DIR, TABLE_DIR, META_DIR, QC_DIR, LOG_DIR)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ── GEO dataset registry ──────────────────────────────────────────────────────
GEO_DATASETS <- list(
  AD  = c("GSE5281",  "GSE63061", "GSE26927"),
  ALS = c("GSE67196", "GSE76220", "GSE26927"),
  PD  = c("GSE6613",  "GSE8397",  "GSE26927")
)

# ── Oxidative stress gene set ─────────────────────────────────────────────────
ANTIOXIDANT_GENES <- list(
  SOD                     = c("SOD1", "SOD2", "SOD3"),
  Catalase                = c("CAT"),
  Glutathione_Peroxidases = c("GPX1", "GPX2", "GPX3", "GPX4", "GPX5", "GPX6"),
  Glutathione_Related     = c("GSR",  "GSTM1", "GSTM2", "GSTM3",
                               "GSTM4", "GSTM5", "GSTP1", "GSTT1"),
  Peroxiredoxins          = c("PRDX1", "PRDX2", "PRDX3", "PRDX4", "PRDX5"),
  Thioredoxin             = c("TXN",  "TXNRD1"),
  Nrf2_Pathway            = c("NFE2L2", "KEAP1"),
  Mitochondrial_Biogenesis = c("NRF1", "TFAM")
)

ALL_ANTIOXIDANT_GENES <- unlist(ANTIOXIDANT_GENES, use.names = FALSE)

# ── Statistical thresholds ────────────────────────────────────────────────────
ALPHA          <- 0.05       # FDR cutoff for significance
LOGFC_THRESH   <- 0          # minimum |logFC| (0 = no filter beyond FDR)
META_METHOD    <- "Stouffer"  # meta-analysis method

# ── Figure export settings ────────────────────────────────────────────────────
FIG_WIDTH  <- 8.5   # inches
FIG_HEIGHT <- 11    # inches
FIG_DPI    <- 600   # raster resolution
FIG_FORMAT <- "tiff"

HEATMAP_COLORS <- colorRampPalette(c("navy", "white", "firebrick3"))(100)

# ── Reproducibility ───────────────────────────────────────────────────────────
set.seed(12345)
options(stringsAsFactors = FALSE, scipen = 999)

message("Config loaded. Project root: ", ROOT_DIR)
