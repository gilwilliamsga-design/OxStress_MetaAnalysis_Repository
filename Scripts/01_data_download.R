###############################################################################
# 01_data_download.R
#
# Project : OxStress Meta-analysis
# Author  : Gil Williams Amana
#
# Purpose : Download all GEO datasets and supplementary files.
#           Skips datasets already present on disk.
#           Also writes Data/GeneLists/antioxidant_genes.csv.
#
# Run after : 00_config.R
###############################################################################

source("Scripts/00_config.R")
source("Functions/io_functions.R")

library(GEOquery)

write_log("Starting 01_data_download.R")

# ── Download GEO series matrix files ─────────────────────────────────────────
all_gse_ids <- unique(unlist(GEO_DATASETS))

cat("\n=== DOWNLOADING GEO DATASETS ===\n\n")

for (gse_id in all_gse_ids) {
  dest <- file.path(RAW_DIR, gse_id)

  # Check whether already downloaded
  cached <- list.files(dest, pattern = "series_matrix", full.names = TRUE)
  if (length(cached) > 0) {
    cat(gse_id, ": already cached – skipping\n")
    next
  }

  cat(gse_id, ": downloading ... ")
  tryCatch({
    getGEO(gse_id, GSEMatrix = TRUE, destdir = RAW_DIR)
    cat("done\n")
  }, error = function(e) {
    cat("FAILED:", conditionMessage(e), "\n")
    write_log(paste("Download failed:", gse_id, conditionMessage(e)))
  })
}

# ── Download supplementary files for RNA-seq datasets ────────────────────────
# GSE67196 – raw count table (ALS)
cat("\nDownloading supplementary file for GSE67196 ...\n")
tryCatch({
  getGEOSuppFiles("GSE67196", baseDir = RAW_DIR)
  cat("done\n")
}, error = function(e) cat("FAILED:", conditionMessage(e), "\n"))

# GSE76220 – RPKM table (ALS)
cat("Downloading supplementary file for GSE76220 ...\n")
tryCatch({
  getGEOSuppFiles("GSE76220", baseDir = RAW_DIR)
  cat("done\n")
}, error = function(e) cat("FAILED:", conditionMessage(e), "\n"))

# ── Write antioxidant gene list ───────────────────────────────────────────────
cat("\nWriting antioxidant gene list ...\n")

gene_df <- data.frame(
  Gene     = ALL_ANTIOXIDANT_GENES,
  Category = rep(names(ANTIOXIDANT_GENES), lengths(ANTIOXIDANT_GENES))
)

save_csv(gene_df, file.path(GENELIST_DIR, "antioxidant_genes.csv"))

cat("\nDownload complete. Datasets in:", RAW_DIR, "\n")
write_log("01_data_download.R complete")
