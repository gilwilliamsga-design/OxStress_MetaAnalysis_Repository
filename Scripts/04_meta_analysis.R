###############################################################################
# 04_meta_analysis.R
#
# Project : OxStress Meta-analysis
# Author  : Gil Williams Amana
#
# Purpose : Combine per-dataset limma results into disease-level meta-analyses
#           using the sign-aware Stouffer method.
#
# Inputs  : Data/Processed/results_*.rds
#
# Outputs : Results/MetaAnalysis/AD_meta_results.csv
#           Results/MetaAnalysis/ALS_meta_results.csv
#           Results/MetaAnalysis/PD_meta_results.csv
#           Data/Processed/{AD,ALS,PD}_meta.rds  (merged logFC/P data frames)
###############################################################################

source("Scripts/00_config.R")
source("Functions/io_functions.R")
source("Functions/stats_functions.R")

library(dplyr)
library(tibble)

write_log("Starting 04_meta_analysis.R")
cat("\n=== META-ANALYSIS (Stouffer method) ===\n\n")

# ── Load all DEA results ──────────────────────────────────────────────────────
load_results <- function(id) {
  r <- load_rds(file.path(PROCESSED_DIR, paste0("results_", id, ".rds")))
  if (!"Gene" %in% colnames(r)) r <- rownames_to_column(r, "Gene")
  r %>% select(Gene, logFC, P.Value)
}

results_5281       <- load_results("5281")
results_63061      <- load_results("63061")
results_26927_AD   <- load_results("26927_AD")
results_67196      <- load_results("67196")
results_76220      <- load_results("76220")
results_26927_ALS  <- load_results("26927_ALS")
results_6613       <- load_results("6613")
results_8397       <- load_results("8397")
results_26927_PD   <- load_results("26927_PD")

# ── Helper: merge three studies on Gene, rename columns ──────────────────────
merge_studies <- function(r1, r2, r3, sfx1, sfx2, sfx3) {
  Reduce(function(a, b) merge(a, b, by = "Gene"), list(r1, r2, r3)) %>%
    setNames(c("Gene",
               paste0("logFC_", sfx1), paste0("P_", sfx1),
               paste0("logFC_", sfx2), paste0("P_", sfx2),
               paste0("logFC_", sfx3), paste0("P_", sfx3)))
}

# ── AD meta ──────────────────────────────────────────────────────────────────
AD_meta <- merge_studies(results_5281, results_63061, results_26927_AD,
                         "5281", "63061", "26927")
save_rds(AD_meta, file.path(PROCESSED_DIR, "AD_meta.rds"))

AD_results <- run_meta_analysis(AD_meta,
                                lfc_cols  = c("logFC_5281", "logFC_63061", "logFC_26927"),
                                pval_cols = c("P_5281",     "P_63061",     "P_26927"))
cat("AD: ", nrow(AD_results[AD_results$FDR < ALPHA, ]),
    "significant gene(s) at FDR <", ALPHA, "\n")
print(AD_results[AD_results$FDR < ALPHA, ])

save_rds(AD_results,  file.path(PROCESSED_DIR,  "AD_results.rds"))
save_csv(AD_results,  file.path(META_DIR, "AD_meta_results.csv"))

# ── ALS meta ─────────────────────────────────────────────────────────────────
ALS_meta <- merge_studies(results_67196, results_76220, results_26927_ALS,
                          "67196", "76220", "26927")
save_rds(ALS_meta, file.path(PROCESSED_DIR, "ALS_meta.rds"))

ALS_results <- run_meta_analysis(ALS_meta,
                                 lfc_cols  = c("logFC_67196", "logFC_76220", "logFC_26927"),
                                 pval_cols = c("P_67196",     "P_76220",     "P_26927"))
cat("ALS:", nrow(ALS_results[ALS_results$FDR < ALPHA, ]),
    "significant gene(s) at FDR <", ALPHA, "\n")
print(ALS_results[ALS_results$FDR < ALPHA, ])

save_rds(ALS_results, file.path(PROCESSED_DIR,  "ALS_results.rds"))
save_csv(ALS_results, file.path(META_DIR, "ALS_meta_results.csv"))

# ── PD meta ──────────────────────────────────────────────────────────────────
PD_meta <- merge_studies(results_6613, results_8397, results_26927_PD,
                         "6613", "8397", "26927")
save_rds(PD_meta, file.path(PROCESSED_DIR, "PD_meta.rds"))

PD_results <- run_meta_analysis(PD_meta,
                                lfc_cols  = c("logFC_6613", "logFC_8397", "logFC_26927"),
                                pval_cols = c("P_6613",     "P_8397",     "P_26927"))
cat("PD: ", nrow(PD_results[PD_results$FDR < ALPHA, ]),
    "significant gene(s) at FDR <", ALPHA, "\n")
print(PD_results[PD_results$FDR < ALPHA, ])

save_rds(PD_results, file.path(PROCESSED_DIR,  "PD_results.rds"))
save_csv(PD_results, file.path(META_DIR, "PD_meta_results.csv"))

write_log("04_meta_analysis.R complete")
cat("\nMeta-analysis complete. Results in:", META_DIR, "\n")
