###############################################################################
# 02_preprocessing.R
#
# Project : OxStress Meta-analysis
# Author  : Gil Williams Amana
#
# Purpose : Load each GEO dataset, extract antioxidant genes, apply log2
#           transformation where necessary, and save per-dataset expression
#           matrices as RDS files in Data/Processed/.
#
# Inputs  : GEO series matrix files (Data/Raw/)
#           GSE67196 raw count .txt.gz
#           GSE76220 RPKM .txt.gz
#           Data/GeneLists/antioxidant_genes.csv
#
# Outputs : Data/Processed/expr_<gse_id>.rds  (one per dataset)
###############################################################################

source("Scripts/00_config.R")
source("Functions/io_functions.R")
source("Functions/expression_functions.R")

library(GEOquery)
library(Biobase)
library(dplyr)

write_log("Starting 02_preprocessing.R")

cat("\n=== PREPROCESSING ===\n\n")

# ── Load gene list ────────────────────────────────────────────────────────────
gene_df <- load_csv(file.path(GENELIST_DIR, "antioxidant_genes.csv"))
genes   <- gene_df$Gene

# ── Microarray datasets (ExpressionSet) ──────────────────────────────────────
microarray_ids <- c("GSE5281", "GSE63061", "GSE6613", "GSE8397",
                    "GSE26927", "GSE67196_array")

process_microarray <- function(gse_id) {
  cat(gse_id, ": loading ... ")
  gse  <- getGEO(gse_id, GSEMatrix = TRUE, destdir = RAW_DIR,\n                 returnType = "ExpressionSet")
  obj  <- gse[[1]]
  expr <- extract_antioxidants(obj, genes)
  expr <- maybe_log2(expr)
  cat(nrow(expr), "genes x", ncol(expr), "samples\n")
  expr
}

# GSE5281 – Alzheimer's Disease
expr_5281 <- process_microarray("GSE5281")
save_rds(expr_5281, file.path(PROCESSED_DIR, "expr_5281.rds"))

# GSE63061 – Alzheimer's Disease
expr_63061 <- process_microarray("GSE63061")
save_rds(expr_63061, file.path(PROCESSED_DIR, "expr_63061.rds"))

# GSE26927 – Multi-disease (AD / ALS / PD subsets used in DEA)
expr_26927 <- process_microarray("GSE26927")
save_rds(expr_26927, file.path(PROCESSED_DIR, "expr_26927.rds"))

# GSE6613 – Parkinson's Disease (blood microarray)
expr_6613 <- process_microarray("GSE6613")
save_rds(expr_6613, file.path(PROCESSED_DIR, "expr_6613.rds"))

# GSE8397 – Parkinson's Disease (substantia nigra)
expr_8397 <- process_microarray("GSE8397")
save_rds(expr_8397, file.path(PROCESSED_DIR, "expr_8397.rds"))

# ── GSE67196 – ALS raw count table ───────────────────────────────────────────
cat("GSE67196 (RNA-seq counts): loading supplementary file ... ")

counts_file <- file.path(RAW_DIR, "GSE67196",
                         "GSE67196_Petrucelli2015_ALS_genes.rawcount.txt.gz")
check_file(counts_file)

als_counts <- read.delim(gzfile(counts_file), check.names = FALSE)

expr_67196 <- als_counts %>%
  filter(GeneID %in% genes) %>%
  { mat <- as.data.frame(.); rownames(mat) <- mat$GeneID; mat[, -(1:2)] } %>%
  as.matrix()

mode(expr_67196) <- "numeric"
expr_67196 <- log2(expr_67196 + 1)   # raw counts always need log2

cat(nrow(expr_67196), "genes x", ncol(expr_67196), "samples\n")
save_rds(expr_67196, file.path(PROCESSED_DIR, "expr_67196.rds"))

# ── GSE76220 – ALS RPKM table ────────────────────────────────────────────────
cat("GSE76220 (RNA-seq RPKM): loading supplementary file ... ")

rpkm_file <- file.path(RAW_DIR, "GSE76220",
                       "GSE76220_ALS_LCM_RPKM.txt.gz")
check_file(rpkm_file)

als_rpkm <- read.delim(gzfile(rpkm_file), check.names = FALSE)

expr_76220 <- als_rpkm %>%
  filter(Gene_symbol %in% genes) %>%
  { mat <- as.data.frame(.)
    rownames(mat) <- mat$Gene_symbol
    mat %>% select(-Gene_symbol, -NCBI_mRNA, -NCBI_Protein, -Gene_Name) } %>%
  as.matrix()

mode(expr_76220) <- "numeric"

cat(nrow(expr_76220), "genes x", ncol(expr_76220), "samples\n")
save_rds(expr_76220, file.path(PROCESSED_DIR, "expr_76220.rds"))

# ── Identify common genes across all datasets ─────────────────────────────────
cat("\nIdentifying common genes across all datasets ...\n")

all_expr <- list(
  GSE5281  = expr_5281,
  GSE63061 = expr_63061,
  GSE26927 = expr_26927,
  GSE6613  = expr_6613,
  GSE8397  = expr_8397,
  GSE67196 = expr_67196,
  GSE76220 = expr_76220
)

common_genes <- Reduce(intersect, lapply(all_expr, rownames))
cat(length(common_genes), "genes shared across all datasets:\n")
cat(paste(sort(common_genes), collapse = ", "), "\n\n")

save_rds(common_genes, file.path(PROCESSED_DIR, "common_genes.rds"))
save_csv(data.frame(Gene = common_genes),
         file.path(PROCESSED_DIR, "common_genes.csv"))

write_log(paste("02_preprocessing.R complete.", length(common_genes), "common genes."))
cat("Preprocessing complete.\n")
