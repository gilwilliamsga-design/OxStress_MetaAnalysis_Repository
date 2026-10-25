###############################################################################
# stats_functions.R
#
# Project : OxStress Meta-analysis
# Author  : Gil Williams Amana
#
# Purpose : Statistical helper functions:
#             - Stouffer meta-analysis
#             - DEG table formatting
#             - Biomarker scoring
###############################################################################

# ── Unweighted Stouffer meta-analysis ───────────────────────────────────────
#' Combine log-fold-changes and p-values from multiple studies using the
#' sign-aware Stouffer Z-score method.
#'
#' @param logfcs  Numeric vector of log2 fold-changes (one per study).
#' @param pvals   Numeric vector of raw p-values (same length as logfcs).
#' @return Named numeric vector: Meta_logFC, Meta_Z, Meta_P.
meta_stouffer <- function(logfcs, pvals) {
  if (length(logfcs) != length(pvals))
    stop("logfcs and pvals must be the same length.")
  if (length(logfcs) == 0L)
    stop("At least one study is required.")
  if (any(!is.finite(logfcs)) || any(!is.finite(pvals)))
    stop("logFC and p-values must be finite; missing study results are not silently dropped.")
  if (any(pvals < 0 | pvals > 1))
    stop("P-values must be between 0 and 1.")

  # Avoid infinite Z-scores when software reports p = 0 due to numerical
  # underflow. This is a numerical safeguard, not additional evidence.
  p_safe <- pmax(pvals, .Machine$double.xmin)
  signs  <- sign(logfcs)
  z      <- qnorm(p_safe / 2, lower.tail = FALSE) * signs
  z_meta <- sum(z) / sqrt(length(z))
  p_meta <- 2 * pnorm(abs(z_meta), lower.tail = FALSE)
  # Descriptive summary: arithmetic mean of study-level logFCs (not a pooled effect).
  lfc    <- mean(logfcs)

  c(Meta_logFC = lfc, Meta_Z = z_meta, Meta_P = p_meta)
}

# ── Apply meta-analysis across all genes in a merged data frame ──────────────
#' @param meta_df   Data frame with columns Gene + pairs of logFC_* / P_*.
#' @param lfc_cols  Character vector of logFC column names (matching pval_cols).
#' @param pval_cols Character vector of P-value column names.
#' @return Data frame with Gene, Meta_logFC, Meta_Z, Meta_P, FDR.
run_meta_analysis <- function(meta_df, lfc_cols, pval_cols) {
  results <- t(apply(meta_df, 1, function(row) {
    lfc <- as.numeric(row[lfc_cols])
    pv  <- as.numeric(row[pval_cols])
    meta_stouffer(lfc, pv)
  }))

  out <- data.frame(Gene = meta_df$Gene, results)
  out$FDR <- p.adjust(out$Meta_P, method = "BH")
  out[order(out$FDR), ]
}

# ── Tidy DEG result table ─────────────────────────────────────────────────────
#' Ensure a limma topTable result has a Gene column and is sorted by FDR.
make_deg_table <- function(x) {
  if (!"Gene" %in% colnames(x))
    x <- tibble::rownames_to_column(x, "Gene")
  dplyr::arrange(x, adj.P.Val)
}

# ── Biomarker scoring ─────────────────────────────────────────────────────────
#' Add Disease_Count, Mean_logFC, and Score columns to a cross-disease summary.
#' @param summary_df  Data frame with columns *_FDR and *_logFC for each disease.
#' @param fdr_cols    Column names holding per-disease FDR values.
#' @param lfc_cols    Column names holding per-disease meta logFC values.
add_biomarker_score <- function(summary_df, fdr_cols, lfc_cols) {
  fdr_mat <- summary_df[, fdr_cols, drop = FALSE]
  lfc_mat <- summary_df[, lfc_cols, drop = FALSE]

  summary_df$Disease_Count <- rowSums(fdr_mat < 0.05, na.rm = TRUE)
  summary_df$Mean_logFC    <- rowMeans(abs(lfc_mat), na.rm = TRUE)
  summary_df$Score         <- summary_df$Disease_Count * summary_df$Mean_logFC

  dplyr::arrange(summary_df, dplyr::desc(Score))
}
