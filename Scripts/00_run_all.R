###############################################################################
# 00_run_all.R
#
# Project : OxStress Meta-analysis
# Author  : Gil Williams Amana
#
# Purpose : Master runner — executes the full pipeline in order.
#           Safe to re-run; steps that have already saved their outputs
#           will still re-execute (change to if(!file.exists()) blocks
#           to enable caching).
###############################################################################

cat("╔══════════════════════════════════════════════════════╗\n")
cat("║  OxStress Meta-analysis Pipeline                    ║\n")
cat("║  Author: Gil Williams Amana                         ║\n")
cat("╚══════════════════════════════════════════════════════╝\n\n")

steps <- list(
  list(script = "Scripts/00_config.R",                 label = "Configuration"),
  list(script = "Scripts/01_data_download.R",          label = "Data Download"),
  list(script = "Scripts/02_preprocessing.R",          label = "Preprocessing"),
  list(script = "Scripts/03_differential_expression.R",label = "Differential Expression"),
  list(script = "Scripts/04_meta_analysis.R",          label = "Meta-analysis"),
  list(script = "Scripts/05_cross_disease.R",          label = "Cross-disease Integration"),
  list(script = "Scripts/06_figures.R",                label = "Figure Generation")
)

for (step in steps) {
  cat(sprintf("\n[%s]\n", step$label))
  cat(rep("-", 50), "\n", sep = "")
  tryCatch(
    source(step$script),
    error = function(e) {
      cat("ERROR in", step$script, ":\n  ", conditionMessage(e), "\n")
      cat("Stopping pipeline.\n")
      stop(conditionMessage(e))
    }
  )
  cat(sprintf("[%s] DONE\n", step$label))
}

cat("\n╔══════════════════════════════════════════════════════╗\n")
cat("║  Pipeline complete. See Results/ for outputs.       ║\n")
cat("╚══════════════════════════════════════════════════════╝\n")
