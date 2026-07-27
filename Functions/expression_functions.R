###############################################################################
# expression_functions.R
#
# Project : OxStress Meta-analysis
# Author  : Gil Williams Amana
#
# Purpose : Functions for extracting antioxidant genes from GEO ExpressionSet
#           objects, normalising, and z-scoring expression matrices.
###############################################################################

library(dplyr)
library(Biobase)

# ── Gene-symbol column detection ──────────────────────────────────────────────
#' Try several common column names to find gene symbols in fData().
#' Returns a character vector the same length as nrow(fData).
get_gene_symbols <- function(gse_obj) {
  fd  <- Biobase::fData(gse_obj)
  candidates <- c("GENE_SYMBOL", "Gene Symbol", "Symbol", "Gene.symbol",
                  "gene_symbol", "SYMBOL")
  col <- Filter(function(x) x %in% colnames(fd), candidates)[1]
  if (is.na(col))
    stop("Could not detect gene-symbol column. Available columns:\n  ",
         paste(colnames(fd), collapse = ", "))
  toupper(trimws(fd[[col]]))
}

# ── Extract & average antioxidant genes from an ExpressionSet ────────────────
#' @param gse_obj  An ExpressionSet (one element from getGEO() list).
#' @param genes    Character vector of gene symbols to extract.
#' @return         Matrix: rows = unique genes, columns = samples.
extract_antioxidants <- function(gse_obj, genes) {
  expr         <- Biobase::exprs(gse_obj)
  gene_symbols <- get_gene_symbols(gse_obj)
  idx          <- gene_symbols %in% toupper(genes)

  if (sum(idx) == 0)
    warning("No antioxidant genes found in this dataset.")

  tmp       <- as.data.frame(expr[idx, , drop = FALSE])
  tmp$Gene  <- gene_symbols[idx]

  # Average probes mapping to the same gene symbol
  tmp <- tmp %>%
    group_by(Gene) %>%
    summarise(across(everything(), mean, na.rm = TRUE)) %>%
    as.data.frame()

  rownames(tmp) <- tmp$Gene
  tmp$Gene      <- NULL
  as.matrix(tmp)
}

# ── Log2 transformation (if data appears to be on raw scale) ─────────────────
#' Applies log2(x + 1) only when max(x) > 100 (raw counts / intensities).
maybe_log2 <- function(mat) {
  if (max(mat, na.rm = TRUE) > 100) {
    message("  Data appears un-logged (max = ", round(max(mat, na.rm = TRUE)),
            "). Applying log2(x + 1).")
    mat <- log2(mat + 1)
  } else {
    message("  Data appears already log-scaled. No transformation applied.")
  }
  mat
}

# ── Row-wise z-score ──────────────────────────────────────────────────────────
#' Z-score each gene (row) across samples (columns).
zscore_matrix <- function(mat) {
  t(scale(t(mat)))
}

# ── Quantile normalisation across arrays ─────────────────────────────────────
#' Thin wrapper around limma::normalizeBetweenArrays.
quantile_normalise <- function(mat) {
  limma::normalizeBetweenArrays(mat, method = "quantile")
}
