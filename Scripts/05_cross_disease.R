###############################################################################
# 05_cross_disease.R
#
# Project : OxStress Meta-analysis
# Author  : Gil Williams Amana
#
# Purpose : Integrate meta-analysis results across AD, ALS and PD.
#           Rank candidate biomarkers by cross-disease significance and
#           effect size. Identify genes significant in all three diseases.
#
# Inputs  : Data/Processed/{AD,ALS,PD}_results.rds
#
# Outputs : Results/Tables/CrossDisease_Summary.csv
#           Results/Tables/Candidate_Biomarker_Ranking.csv
#           Results/Tables/Common_Biomarkers_All3.csv
###############################################################################

source("Scripts/00_config.R")
source("Functions/io_functions.R")
source("Functions/stats_functions.R")

library(dplyr)
library(tibble)

write_log("Starting 05_cross_disease.R")
cat("\n=== CROSS-DISEASE INTEGRATION ===\n\n")

# ── Load meta-analysis results ────────────────────────────────────────────────
AD_results  <- load_rds(file.path(PROCESSED_DIR, "AD_results.rds"))
ALS_results <- load_rds(file.path(PROCESSED_DIR, "ALS_results.rds"))
PD_results  <- load_rds(file.path(PROCESSED_DIR, "PD_results.rds"))


# ── Build cross-disease summary table ────────────────────────────────────────
cross_disease <- AD_results %>%
  dplyr::select(Gene, AD_logFC = Meta_logFC, AD_FDR = FDR) %>%
  full_join(
    ALS_results %>%
      dplyr::select(Gene, ALS_logFC = Meta_logFC, ALS_FDR = FDR),
    by = "Gene"
  ) %>%
  full_join(
    PD_results %>%
      dplyr::select(Gene, PD_logFC = Meta_logFC, PD_FDR = FDR),
    by = "Gene"
  ) %>%
  mutate(
    AD_sig  = ifelse(!is.na(AD_FDR)  & AD_FDR  < ALPHA, "Yes", "No"),
    ALS_sig = ifelse(!is.na(ALS_FDR) & ALS_FDR < ALPHA, "Yes", "No"),
    PD_sig  = ifelse(!is.na(PD_FDR)  & PD_FDR  < ALPHA, "Yes", "No")
  )

save_csv(cross_disease, file.path(TABLE_DIR, "CrossDisease_Summary.csv"))

# ── Biomarker scoring and ranking ─────────────────────────────────────────────
ranking <- add_biomarker_score(
  cross_disease,
  fdr_cols = c("AD_FDR", "ALS_FDR", "PD_FDR"),
  lfc_cols = c("AD_logFC", "ALS_logFC", "PD_logFC")
)

save_csv(ranking, file.path(TABLE_DIR, "Candidate_Biomarker_Ranking.csv"))

# ── Genes significant in all three diseases ───────────────────────────────────
common_biomarkers <- ranking %>%
  filter(Disease_Count == 3)

cat("=== SHARED BIOMARKERS (All 3 Diseases) ===\n")
print(common_biomarkers[, c("Gene", "AD_logFC", "ALS_logFC", "PD_logFC",
                             "AD_FDR", "ALS_FDR", "PD_FDR", "Score")])

save_csv(common_biomarkers, file.path(TABLE_DIR, "Common_Biomarkers_All3.csv"))

# ── Summary of direction ──────────────────────────────────────────────────────
direction_summary <- cross_disease %>%
  mutate(
    Dir_AD  = case_when(AD_sig == "Yes" & AD_logFC > 0 ~ "Up",
                        AD_sig == "Yes" & AD_logFC < 0 ~ "Down",
                        TRUE ~ "NS"),
    Dir_ALS = case_when(ALS_sig == "Yes" & ALS_logFC > 0 ~ "Up",
                        ALS_sig == "Yes" & ALS_logFC < 0 ~ "Down",
                        TRUE ~ "NS"),
    Dir_PD  = case_when(PD_sig == "Yes" & PD_logFC > 0 ~ "Up",
                        PD_sig == "Yes" & PD_logFC < 0 ~ "Down",
                        TRUE ~ "NS")
  ) %>%
  select(Gene, Dir_AD, Dir_ALS, Dir_PD)

cat("\n=== DIRECTIONAL CONSISTENCY ===\n")
print(direction_summary %>% filter(Dir_AD != "NS" | Dir_ALS != "NS" | Dir_PD != "NS"))

save_csv(direction_summary, file.path(TABLE_DIR, "Direction_Summary.csv"))

write_log("05_cross_disease.R complete")
cat("\nCross-disease integration complete. Tables in:", TABLE_DIR, "\n")
