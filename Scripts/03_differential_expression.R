###############################################################################
# 03_differential_expression.R
#
# Project : OxStress Meta-analysis
# Author  : Gil Williams Amana
#
# Purpose : Run limma differential expression for each dataset.
#           Disease vs. Control comparisons are defined here.
#           Results saved as RDS + CSV per dataset.
#
# Inputs  : Data/Processed/expr_*.rds
#           GEO pData (re-fetched from cache in Data/Raw/)
#
# Outputs : Results/Tables/<gse_id>_DEG.csv
#           Data/Processed/results_<gse_id>.rds
###############################################################################

source("Scripts/00_config.R")
source("Functions/io_functions.R")
source("Functions/expression_functions.R")
source("Functions/stats_functions.R")

# Check required packages before running the analysis
required_packages <- c("GEOquery", "Biobase", "limma", "edgeR",
                       "dplyr", "tibble")

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1),
          quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required packages: ",
    paste(missing_packages, collapse = ", "),
    ". Install them before running this script."
  )
}

library(GEOquery)
library(Biobase)
library(limma)
library(dplyr)
library(tibble)
write_log("Starting 03_differential_expression.R")

cat("\n=== DIFFERENTIAL EXPRESSION ANALYSIS ===\n\n")

common_genes <- load_rds(file.path(PROCESSED_DIR, "common_genes.rds"))

# Helper: extract common-gene rows from a processed matrix
subset_common <- function(mat) mat[common_genes, , drop = FALSE]

# ── DEA wrapper ───────────────────────────────────────────────────────────────
#' Run limma on a matrix with a single binary factor (Case vs Control).
#' @param expr      Gene × sample matrix (common genes, log scale).
#' @param group     Factor with levels c("Control", <case_label>).
#' @param covariate Optional named list of additional covariates (e.g. region).
#' @param coef      Coefficient to extract; defaults to the case level.
run_limma <- function(expr, group, covariate = NULL, coef = NULL) {
  expr <- as.matrix(expr)
  storage.mode(expr) <- "numeric"

  if (ncol(expr) != length(group)) {
    stop("Expression columns (", ncol(expr),
         ") do not match group labels (", length(group), ").")
  }
  if (nrow(expr) < 1L || ncol(expr) < 2L) {
    stop("Expression matrix must contain at least one gene and two samples.")
  }
  if (any(!is.finite(expr))) {
    stop("Expression matrix contains NA, NaN, or infinite values.")
  }

  group <- droplevels(factor(group))
  if (anyNA(group) || nlevels(group) != 2L) {
    stop("Differential expression requires exactly two observed groups; found: ",
         paste(levels(group), collapse = ", "))
  }
  if (any(table(group) < 2L)) {
    stop("Each group must contain at least two samples for this limma workflow.")
  }

  if (!is.null(covariate)) {
    if (is.null(names(covariate)) || any(!nzchar(names(covariate)))) {
      stop("Covariates must be supplied as a named list.")
    }
    if (any(vapply(covariate, length, integer(1)) != ncol(expr))) {
      stop("Each covariate must have one value per expression sample.")
    }
    if (any(vapply(covariate, function(x) anyNA(x), logical(1)))) {
      stop("Covariates contain missing values.")
    }
  }

  df_design           <- data.frame(Group = group)
  formula_str         <- "~ Group"

  if (!is.null(covariate)) {
    for (nm in names(covariate)) {
      df_design[[nm]]  <- covariate[[nm]]
      formula_str      <- paste(nm, "+", formula_str)
    }
    formula_str <- paste("~", formula_str)  # rebuild properly
    formula_str <- paste0("~ ", paste(names(covariate), collapse = " + "), " + Group")
  }

  design <- model.matrix(as.formula(formula_str), data = df_design)
  if (is.null(coef)) coef <- grep("Group", colnames(design), value = TRUE)

  fit <- lmFit(expr, design)
  fit <- eBayes(fit)
  topTable(fit, coef = coef, number = Inf, sort.by = "P")
}

# ═══════════════════════════════════════════════════════════════════════════════
# ALZHEIMER'S DISEASE
# ═══════════════════════════════════════════════════════════════════════════════

# ── GSE5281 ───────────────────────────────────────────────────────────────────
cat("GSE5281 (AD) ...\n")
gse5281     <- getGEO("GSE5281", GSEMatrix = TRUE, destdir = RAW_DIR)
meta_5281   <- pData(gse5281[[1]])
expr_5281   <- subset_common(load_rds(file.path(PROCESSED_DIR, "expr_5281.rds")))

meta_5281$Disease <- ifelse(grepl("Alzheimer", meta_5281[["Disease State:ch1"]],
                                  ignore.case = TRUE), "AD", "Control")

# Balance regions
keep_regions <- c("Entorhinal Cortex", "Primary Visual Cortex", "Superior Frontal Gyrus")
meta_5281$RegionName <- as.character(meta_5281[["Organ Region:ch1"]])
meta_5281$RegionName <- gsub("\u00A0", " ", meta_5281$RegionName, fixed = TRUE)
meta_5281$RegionName <- trimws(meta_5281$RegionName)
keep_idx <- meta_5281$RegionName %in% keep_regions

meta_5281_bal <- meta_5281[keep_idx, ]
expr_5281_bal <- expr_5281[, keep_idx, drop = FALSE]

group_5281 <- factor(
  meta_5281_bal$Disease,
  levels = c("Control", "AD")
)
region_5281 <- factor(meta_5281_bal$RegionName)

results_5281 <- run_limma(expr_5281_bal, group_5281,
                          covariate = list(Region = region_5281))
cat("  Significant genes (FDR<0.05):", sum(results_5281$adj.P.Val < ALPHA), "\n")

save_rds(results_5281, file.path(PROCESSED_DIR, "results_5281.rds"))
save_csv(make_deg_table(results_5281), file.path(TABLE_DIR, "GSE5281_DEG.csv"))

# ── GSE63061 ──────────────────────────────────────────────────────────────────
cat("GSE63061 (AD) ...\n")
gse63061    <- getGEO("GSE63061", GSEMatrix = TRUE, destdir = RAW_DIR)
meta_63061  <- pData(gse63061[[1]])
expr_63061  <- subset_common(load_rds(file.path(PROCESSED_DIR, "expr_63061.rds")))

# Keep case-control samples only
keep_63061 <- meta_63061[["included in case -control study:ch1"]] == "yes" &
              meta_63061[["status:ch1"]] %in% c("AD", "CTL")
meta_63061_bal <- meta_63061[keep_63061, ]
expr_63061_bal <- expr_63061[, keep_63061, drop = FALSE]

group_63061 <- factor(meta_63061_bal[["status:ch1"]], levels = c("CTL", "AD"))
results_63061 <- run_limma(expr_63061_bal, group_63061)
cat("  Significant genes (FDR<0.05):", sum(results_63061$adj.P.Val < ALPHA), "\n")

save_rds(results_63061, file.path(PROCESSED_DIR, "results_63061.rds"))
save_csv(make_deg_table(results_63061), file.path(TABLE_DIR, "GSE63061_DEG.csv"))

# ── GSE26927 – AD subset ──────────────────────────────────────────────────────
cat("GSE26927 (AD subset) ...\n")
gse26927    <- getGEO("GSE26927", GSEMatrix = TRUE, destdir = RAW_DIR)
meta_26927  <- pData(gse26927[[1]])
expr_26927  <- subset_common(load_rds(file.path(PROCESSED_DIR, "expr_26927.rds")))

meta_26927_AD <- meta_26927[meta_26927[["disease:ch1"]] == "Alzheimer's disease", ]
meta_26927_AD$Group <- ifelse(grepl("Disease", meta_26927_AD$title), "AD", "Control")

expr_26927_AD <- expr_26927[, rownames(meta_26927_AD), drop = FALSE]


group_26927_AD <- factor(meta_26927_AD$Group, levels = c("Control", "AD"))
results_26927_AD <- run_limma(expr_26927_AD, group_26927_AD)
cat("  Significant genes (FDR<0.05):", sum(results_26927_AD$adj.P.Val < ALPHA), "\n")

save_rds(results_26927_AD, file.path(PROCESSED_DIR, "results_26927_AD.rds"))
save_csv(make_deg_table(results_26927_AD),
         file.path(TABLE_DIR, "GSE26927_AD_DEG.csv"))

# ═══════════════════════════════════════════════════════════════════════════════
# AMYOTROPHIC LATERAL SCLEROSIS
# ═══════════════════════════════════════════════════════════════════════════════

# ── GSE67196 ──────────────────────────────────────────────────────────────────
cat("GSE67196 (ALS) — edgeR + limma-voom, frontal cortex only...\n")

# Load GEO metadata
gse67196 <- getGEO("GSE67196", GSEMatrix = TRUE, destdir = RAW_DIR)
meta_67196 <- pData(gse67196[[1]])

# Load original raw counts
count_file <- file.path(
  RAW_DIR, "GSE67196",
  "GSE67196_Petrucelli2015_ALS_genes.rawcount.txt.gz"
)
als_raw <- read.delim(
  count_file, check.names = FALSE, stringsAsFactors = FALSE
)

# Select frontal-cortex count columns
fcx_samples <- grep("_fcx$", names(als_raw), value = TRUE)
if (length(fcx_samples) == 0L) {
  stop("No frontal-cortex (_fcx) columns found.")
}

# Clean gene IDs and sum duplicate IDs
gene_ids <- trimws(as.character(als_raw$GeneID))
valid <- !is.na(gene_ids) & nzchar(gene_ids)

fcx_counts <- as.matrix(als_raw[valid, fcx_samples, drop = FALSE])
gene_ids <- gene_ids[valid]
fcx_counts <- rowsum(fcx_counts, group = gene_ids, reorder = FALSE)
storage.mode(fcx_counts) <- "integer"

# Select frontal-cortex metadata and match by numeric sample ID
fcx_meta <- meta_67196[
  grepl("_FCX$", meta_67196$title, ignore.case = TRUE),
  ,
  drop = FALSE
]

count_ids <- as.integer(sub("^ALS", "", sub("_fcx$", "", colnames(fcx_counts))))
geo_ids <- as.integer(sub("_FCX$", "", fcx_meta$title, ignore.case = TRUE))

match_idx <- match(count_ids, geo_ids)
if (anyNA(match_idx) || anyDuplicated(geo_ids)) {
  stop("GSE67196 sample matching failed; check count columns and GEO metadata.")
}

fcx_meta <- fcx_meta[match_idx, , drop = FALSE]
rownames(fcx_meta) <- colnames(fcx_counts)

# Define Control vs ALS; Control is the reference group
genotype <- trimws(as.character(fcx_meta[["genotype:ch1"]]))
allowed_genotypes <- c("Healthy", "sALS", "c9ALS")

if (anyNA(genotype) || any(!genotype %in% allowed_genotypes)) {
  stop("Unexpected or missing GSE67196 genotype labels: ",
       paste(unique(genotype[is.na(genotype) |
         !genotype %in% allowed_genotypes]), collapse = ", "))
}

fcx_meta$Group <- factor(
  ifelse(genotype == "Healthy", "Control", "ALS"),
  levels = c("Control", "ALS")
)

if (anyNA(fcx_meta$Group) || length(unique(fcx_meta$Group)) != 2L) {
  stop("GSE67196 must contain both Control and ALS samples.")
}
stopifnot(identical(colnames(fcx_counts), rownames(fcx_meta)))

# edgeR: filter raw counts and normalize library sizes
y <- edgeR::DGEList(counts = fcx_counts, group = fcx_meta$Group)
keep <- edgeR::filterByExpr(y, group = fcx_meta$Group)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- edgeR::calcNormFactors(y, method = "TMM")

# limma-voom differential expression
design <- stats::model.matrix(~ Group, data = fcx_meta)
v <- limma::voom(y, design, plot = FALSE)
fit <- limma::lmFit(v, design)
fit <- limma::eBayes(fit)

results_all <- limma::topTable(
  fit, coef = "GroupALS", number = Inf, sort.by = "none"
)
results_all$Gene <- rownames(results_all)

# Retain the 18 genes used by the meta-analysis
common_genes <- c(
  "CAT", "GPX1", "GPX3", "GPX4", "GSR", "GSTM1", "GSTM3", "GSTM4",
  "GSTP1", "GSTT1", "KEAP1", "NRF1", "PRDX1", "PRDX2", "PRDX3",
  "PRDX4", "SOD1", "TFAM"
)

results_67196 <- results_all[
  results_all$Gene %in% common_genes,
  c("Gene", "logFC", "P.Value", "adj.P.Val"),
  drop = FALSE
]

# Preserve common-gene ordering
results_67196 <- results_67196[
  match(common_genes, results_67196$Gene),
  ,
  drop = FALSE
]
results_67196 <- results_67196[!is.na(results_67196$Gene), , drop = FALSE]

cat("  Frontal-cortex samples:", ncol(fcx_counts), "\n")
cat("  Genes retained after filtering:", nrow(y), "\n")
cat("  Oxidative-stress genes found:", nrow(results_67196), "\n")
cat("  Significant genes (FDR<0.05):",
    sum(results_67196$adj.P.Val < ALPHA, na.rm = TRUE), "\n")

save_rds(results_67196, file.path(PROCESSED_DIR, "results_67196.rds"))
save_csv(make_deg_table(results_67196), file.path(TABLE_DIR, "GSE67196_DEG.csv"))

# ── GSE76220 ──────────────────────────────────────────────────────────────────
cat("GSE76220 (ALS) ...\n")
gse76220    <- getGEO("GSE76220", GSEMatrix = TRUE, destdir = RAW_DIR)
meta_76220  <- pData(gse76220[[1]])
expr_76220  <- subset_common(load_rds(file.path(PROCESSED_DIR, "expr_76220.rds")))

meta_76220$Group <- factor(
  ifelse(grepl("Control", meta_76220$title), "Control", "ALS"),
  levels = c("Control", "ALS"))

results_76220 <- run_limma(expr_76220, meta_76220$Group)
# Direction must be verified from the GEO sample annotations before interpreting logFC.
cat("  Significant genes (FDR<0.05):", sum(results_76220$adj.P.Val < ALPHA), "\n")

save_rds(results_76220, file.path(PROCESSED_DIR, "results_76220.rds"))
save_csv(make_deg_table(results_76220), file.path(TABLE_DIR, "GSE76220_DEG.csv"))

# ── GSE26927 – ALS subset ─────────────────────────────────────────────────────
cat("GSE26927 (ALS subset) ...\n")
meta_26927_ALS <- meta_26927[meta_26927[["disease:ch1"]] ==
                               "Amyotrophic lateral sclerosis", ]
meta_26927_ALS$Group <- ifelse(grepl("Disease", meta_26927_ALS$title),
                               "ALS", "Control")
expr_26927_ALS <- expr_26927[, rownames(meta_26927_ALS), drop = FALSE]


group_26927_ALS <- factor(meta_26927_ALS$Group, levels = c("Control", "ALS"))
results_26927_ALS <- run_limma(expr_26927_ALS, group_26927_ALS)
cat("  Significant genes (FDR<0.05):", sum(results_26927_ALS$adj.P.Val < ALPHA), "\n")

save_rds(results_26927_ALS, file.path(PROCESSED_DIR, "results_26927_ALS.rds"))
save_csv(make_deg_table(results_26927_ALS),
         file.path(TABLE_DIR, "GSE26927_ALS_DEG.csv"))

# ═══════════════════════════════════════════════════════════════════════════════
# PARKINSON'S DISEASE
# ═══════════════════════════════════════════════════════════════════════════════

# ── GSE6613 ───────────────────────────────────────────────────────────────────
cat("GSE6613 (PD) ...\n")
gse6613    <- getGEO("GSE6613", GSEMatrix = TRUE, destdir = RAW_DIR)
meta_6613  <- pData(gse6613[[1]])
expr_6613  <- subset_common(load_rds(file.path(PROCESSED_DIR, "expr_6613.rds")))

# Keep only PD patients and healthy controls
keep_6613 <- meta_6613$characteristics_ch1 %in%
             c("Parkinson's disease", "healthy control")
meta_6613_bal <- meta_6613[keep_6613, ]
expr_6613_bal <- expr_6613[, rownames(meta_6613_bal), drop = FALSE]
expr_6613_log <- expr_6613_bal  # already transformed by preprocessing; avoid double-log

group_6613 <- factor(
  ifelse(meta_6613_bal$characteristics_ch1 == "Parkinson's disease",
         "PD", "Control"),
  levels = c("Control", "PD"))

results_6613 <- run_limma(expr_6613_log, group_6613)
cat("  Significant genes (FDR<0.05):", sum(results_6613$adj.P.Val < ALPHA), "\n")

save_rds(results_6613, file.path(PROCESSED_DIR, "results_6613.rds"))
save_csv(make_deg_table(results_6613), file.path(TABLE_DIR, "GSE6613_DEG.csv"))

# ── GSE8397 ───────────────────────────────────────────────────────────────────
cat("GSE8397 (PD) ...\n")
gse8397    <- getGEO("GSE8397", GSEMatrix = TRUE, destdir = RAW_DIR)
meta_8397  <- pData(gse8397[[1]])
expr_8397  <- subset_common(load_rds(file.path(PROCESSED_DIR, "expr_8397.rds")))

group_8397 <- factor(
  ifelse(grepl("Parkinson", meta_8397$title, ignore.case = TRUE),
         "PD", "Control"),
  levels = c("Control", "PD"))

results_8397 <- run_limma(expr_8397, group_8397)
cat("  Significant genes (FDR<0.05):", sum(results_8397$adj.P.Val < ALPHA), "\n")

save_rds(results_8397, file.path(PROCESSED_DIR, "results_8397.rds"))
save_csv(make_deg_table(results_8397), file.path(TABLE_DIR, "GSE8397_DEG.csv"))

# ── GSE26927 – PD subset ──────────────────────────────────────────────────────
cat("GSE26927 (PD subset) ...\n")
meta_26927_PD <- meta_26927[meta_26927[["disease:ch1"]] == "Parkinson's disease", ]
meta_26927_PD$Group <- ifelse(grepl("Disease", meta_26927_PD$title),
                              "PD", "Control")
expr_26927_PD <- expr_26927[, rownames(meta_26927_PD), drop = FALSE]


group_26927_PD <- factor(meta_26927_PD$Group, levels = c("Control", "PD"))
results_26927_PD <- run_limma(expr_26927_PD, group_26927_PD)
cat("  Significant genes (FDR<0.05):", sum(results_26927_PD$adj.P.Val < ALPHA), "\n")

save_rds(results_26927_PD, file.path(PROCESSED_DIR, "results_26927_PD.rds"))
save_csv(make_deg_table(results_26927_PD),
         file.path(TABLE_DIR, "GSE26927_PD_DEG.csv"))

write_log("03_differential_expression.R complete")
cat("\nDifferential expression complete. Results in:", TABLE_DIR, "\n")
