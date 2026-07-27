###############################################################################
# 06_figures.R
#
# Project : OxStress Meta-analysis
# Author  : Gil Williams Amana
#
# Purpose : Generate all publication-quality figures.
#           Each figure is exported as PNG (600 dpi), TIFF (600 dpi, LZW),
#           and PDF.
#
# Figures:
#   Figure 1 – Study Design and Analytical Workflow
#   Figure 2 – Heatmap: meta-analysis logFC across AD, ALS, PD
#   Figure 3 – Venn Diagram: shared significant genes
#   Figure 4A – Forest Plot: SOD1 across all datasets
#   Figure 4B – Forest Plot: PRDX2 across all datasets
#   Figure 5  – Candidate Biomarker Ranking Bar Chart
###############################################################################

source("Scripts/00_config.R")
source("Functions/io_functions.R")
source("Functions/stats_functions.R")

library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(grid)
library(VennDiagram)

write_log("Starting 06_figures.R")
cat("\n=== GENERATING FIGURES ===\n\n")

# ── Load inputs ───────────────────────────────────────────────────────────────
AD_results   <- load_rds(file.path(PROCESSED_DIR, "AD_results.rds"))
ALS_results  <- load_rds(file.path(PROCESSED_DIR, "ALS_results.rds"))
PD_results   <- load_rds(file.path(PROCESSED_DIR, "PD_results.rds"))

AD_meta  <- load_rds(file.path(PROCESSED_DIR, "AD_meta.rds"))
ALS_meta <- load_rds(file.path(PROCESSED_DIR, "ALS_meta.rds"))
PD_meta  <- load_rds(file.path(PROCESSED_DIR, "PD_meta.rds"))

ranking  <- read.csv(file.path(TABLE_DIR, "Candidate_Biomarker_Ranking.csv"),
                     stringsAsFactors = FALSE)

# ── Build heatmap data frame ──────────────────────────────────────────────────
heatmap_df <- AD_results %>%
  select(Gene, AD = Meta_logFC) %>%
  left_join(ALS_results %>% select(Gene, ALS = Meta_logFC), by = "Gene") %>%
  left_join(PD_results  %>% select(Gene, PD  = Meta_logFC), by = "Gene")

# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 1 – Study Design and Analytical Workflow
# ═══════════════════════════════════════════════════════════════════════════════
cat("Figure 1: Study workflow ...\n")

draw_figure1 <- function() {

  # ---- colour palette ----
  COL <- list(
    ad = "#4F81BD", als = "#70AD47", pd = "#C0504D",
    qc = "#E8F1FB", deg = "#EAF7EA", ox = "#FDE9D9",
    int = "#F2F2F2", rank = "#FFF4CC", vis = "#EFE5FA",
    border = "#4A4A4A", text = "#222222", arrow = "#555555"
  )

  draw_box <- function(x, y, w, h, fill = "white", border = COL$border, lwd = 1.3) {
    grid.roundrect(x = unit(x, "npc"), y = unit(y, "npc"),
                   width = unit(w, "npc"), height = unit(h, "npc"),
                   r = unit(0.02, "snpc"), gp = gpar(fill = fill, col = border, lwd = lwd))
  }

  draw_label <- function(x, y, label, size = 10, face = "plain", col = COL$text) {
    grid.text(label, x = unit(x, "npc"), y = unit(y, "npc"),
              gp = gpar(fontsize = size, fontface = face, col = col))
  }

  draw_arrow <- function(x0, y0, x1, y1) {
    grid.lines(x = unit(c(x0, x1), "npc"), y = unit(c(y0, y1), "npc"),
               arrow = arrow(type = "closed", length = unit(0.12, "in")),
               gp = gpar(lwd = 1.8, col = COL$arrow))
  }

  grid.newpage()
  grid.rect(gp = gpar(fill = "white", col = NA))

  draw_label(0.50, 0.965, "Study Design and Analytical Workflow",
             size = 16, face = "bold")

  # Disease header
  draw_box(0.50, 0.915, 0.70, 0.048, fill = "#1F4E79", border = "#1F4E79")
  draw_label(0.50, 0.915, "GEO Transcriptomic Datasets",
             size = 13, face = "bold", col = "white")

  # Disease boxes
  diseases <- list(
    list(cx = 0.18, fill = COL$als,
         title = "Amyotrophic\nLateral Sclerosis",
         body  = "GSE67196 | GSE76220\nGSE26927 (ALS subset)"),
    list(cx = 0.50, fill = COL$ad,
         title = "Alzheimer's\nDisease",
         body  = "GSE5281 | GSE63061\nGSE26927 (AD subset)"),
    list(cx = 0.82, fill = COL$pd,
         title = "Parkinson's\nDisease",
         body  = "GSE6613 | GSE8397\nGSE26927 (PD subset)")
  )
  for (d in diseases) {
    draw_box(d$cx, 0.84, 0.23, 0.10, fill = "white", border = d$fill, lwd = 1.8)
    grid.roundrect(x = unit(d$cx, "npc"), y = unit(0.895, "npc"),
                   width = unit(0.23, "npc"), height = unit(0.025, "npc"),
                   r = unit(0.01, "snpc"), gp = gpar(fill = d$fill, col = d$fill))
    draw_label(d$cx, 0.895, d$title, size = 9.5, face = "bold", col = "white")
    draw_label(d$cx, 0.840, d$body, size = 8.5)
    draw_arrow(d$cx, 0.79, 0.50, 0.752)
  }

  # Horizontal merge rail
  grid.lines(x = unit(c(0.18, 0.82), "npc"), y = unit(c(0.752, 0.752), "npc"),
             gp = gpar(col = COL$arrow, lwd = 1.5))

  # Workflow steps
  steps <- list(
    list(y = 0.72, lab = "Quality Control\nNormalization & Probe Annotation",
         fill = COL$qc),
    list(y = 0.63, lab = "Differential Expression Analysis\n(limma)",
         fill = COL$deg),
    list(y = 0.54, lab = "Oxidative Stress Gene Extraction",
         fill = COL$ox),
    list(y = 0.45, lab = "Cross-Disease Integration",
         fill = COL$int),
    list(y = 0.36, lab = "Candidate Biomarker Prioritization",
         fill = COL$rank),
    list(y = 0.27, lab = "Visualization\nHeatmaps • Forest Plots • Venn Diagrams",
         fill = COL$vis)
  )

  for (i in seq_along(steps)) {
    s <- steps[[i]]
    draw_box(0.50, s$y, 0.42, 0.065, fill = s$fill, border = "#888888", lwd = 1.0)
    draw_label(0.50, s$y, s$lab, size = 10, face = "bold")
    if (i < length(steps))
      draw_arrow(0.50, s$y - 0.033, 0.50, steps[[i + 1]]$y + 0.033)
  }

  # Output box
  draw_box(0.50, 0.10, 0.55, 0.08, fill = "#FFF9E6", border = "#B8860B", lwd = 1.8)
  draw_label(0.50, 0.125, "OUTPUTS", size = 12, face = "bold", col = "#B8860B")
  draw_label(0.50, 0.085,
             "Shared & Disease-specific Oxidative Stress Biomarkers\nHeatmaps  •  Forest Plots  •  Venn Diagram  •  Ranking Tables",
             size = 8.5)
  draw_arrow(0.50, 0.242, 0.50, 0.140)
}

# Export
for (fmt in c("png", "tiff", "pdf")) {
  stem <- file.path(FIGURE_DIR, "Figure1_Study_Workflow")
  if (fmt == "png")
    png(paste0(stem, ".png"), width = FIG_WIDTH, height = FIG_HEIGHT,
        units = "in", res = FIG_DPI, bg = "white")
  else if (fmt == "tiff")
    tiff(paste0(stem, ".tiff"), width = FIG_WIDTH, height = FIG_HEIGHT,
         units = "in", res = FIG_DPI, compression = "lzw", bg = "white")
  else
    pdf(paste0(stem, ".pdf"), width = FIG_WIDTH, height = FIG_HEIGHT, bg = "white")

  draw_figure1()
  dev.off()
}
cat("  Saved Figure 1\n")

# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 2 – Heatmap
# ═══════════════════════════════════════════════════════════════════════════════
cat("Figure 2: Heatmap ...\n")

heatmap_long <- heatmap_df %>%
  pivot_longer(cols = c(AD, ALS, PD), names_to = "Disease", values_to = "Meta_logFC") %>%
  mutate(Disease = factor(Disease, levels = c("AD", "ALS", "PD")))

fig2 <- ggplot(heatmap_long, aes(Disease, Gene, fill = Meta_logFC)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0, name = "Meta logFC"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x  = element_text(face = "bold", size = 13),
    axis.text.y  = element_text(face = "italic", size = 10),
    plot.title   = element_text(face = "bold", size = 15, hjust = 0.5),
    legend.title = element_text(face = "bold")
  ) +
  labs(
    title = "Meta-analysis of Oxidative Stress Genes\nAcross Alzheimer's Disease, ALS and Parkinson's Disease",
    x = "", y = ""
  )

save_publication_plot(fig2, file.path(FIGURE_DIR, "Figure2_MetaAnalysis_Heatmap"),
                      width = 7, height = 9)
cat("  Saved Figure 2\n")

# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 3 – Venn Diagram
# ═══════════════════════════════════════════════════════════════════════════════
cat("Figure 3: Venn diagram ...\n")

AD_genes  <- AD_results$Gene[AD_results$FDR < ALPHA]
ALS_genes <- ALS_results$Gene[ALS_results$FDR < ALPHA]
PD_genes  <- PD_results$Gene[PD_results$FDR < ALPHA]

venn_plot <- venn.diagram(
  x = list("Alzheimer's Disease" = AD_genes,
           "ALS"                  = ALS_genes,
           "Parkinson's Disease"  = PD_genes),
  filename  = NULL,
  fill      = c("#4C72B0", "#55A868", "#C44E52"),
  alpha     = 0.55,
  cex       = 1.8, cat.cex = 1.5, cat.fontface = "bold",
  lwd       = 2.0, margin = 0.12
)

for (fmt in c("png", "tiff")) {
  stem <- file.path(FIGURE_DIR, "Figure3_Venn_SharedBiomarkers")
  if (fmt == "png")
    png(paste0(stem, ".png"), width = 3200, height = 3200, res = 400)
  else
    tiff(paste0(stem, ".tiff"), width = 3200, height = 3200, res = 400,
         compression = "lzw")
  grid.draw(venn_plot)
  dev.off()
}
cat("  Saved Figure 3\n")

# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 4 – Forest Plots
# ═══════════════════════════════════════════════════════════════════════════════
cat("Figure 4: Forest plots ...\n")

build_forest_df <- function(gene) {
  data.frame(
    Dataset = c("GSE5281","GSE63061","GSE26927_AD",
                "GSE67196","GSE76220","GSE26927_ALS",
                "GSE6613","GSE8397","GSE26927_PD",
                "Meta AD","Meta ALS","Meta PD"),
    Disease = c("AD","AD","AD","ALS","ALS","ALS","PD","PD","PD","AD","ALS","PD"),
    logFC = c(
      AD_meta[AD_meta$Gene == gene, "logFC_5281"],
      AD_meta[AD_meta$Gene == gene, "logFC_63061"],
      AD_meta[AD_meta$Gene == gene, "logFC_26927"],
      ALS_meta[ALS_meta$Gene == gene, "logFC_67196"],
      ALS_meta[ALS_meta$Gene == gene, "logFC_76220"],
      ALS_meta[ALS_meta$Gene == gene, "logFC_26927"],
      PD_meta[PD_meta$Gene == gene, "logFC_6613"],
      PD_meta[PD_meta$Gene == gene, "logFC_8397"],
      PD_meta[PD_meta$Gene == gene, "logFC_26927"],
      AD_results[AD_results$Gene == gene, "Meta_logFC"],
      ALS_results[ALS_results$Gene == gene, "Meta_logFC"],
      PD_results[PD_results$Gene == gene, "Meta_logFC"]
    ),
    IsMeta = c(rep(FALSE, 9), rep(TRUE, 3))
  )
}

forest_plot <- function(gene) {
  df <- build_forest_df(gene)
  ggplot(df, aes(x = logFC, y = reorder(Dataset, logFC),
                 color = Disease, shape = IsMeta, size = IsMeta)) +
    geom_vline(xintercept = 0, linetype = 2, color = "grey40") +
    geom_point() +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 18), guide = "none") +
    scale_size_manual(values  = c(`FALSE` = 3,  `TRUE` = 5),  guide = "none") +
    scale_color_manual(values = c(AD = "#4C72B0", ALS = "#55A868", PD = "#C44E52")) +
    theme_classic(base_size = 14) +
    labs(title = paste(gene, "Expression Across Neurodegenerative Diseases"),
         x = "log2 Fold Change", y = "")
}

for (gene in c("SOD1", "PRDX2")) {
  fig_id <- ifelse(gene == "SOD1", "Figure4A", "Figure4B")
  p      <- forest_plot(gene)
  save_publication_plot(p,
                        file.path(FIGURE_DIR, paste0(fig_id, "_", gene, "_ForestPlot")),
                        width = 8, height = 6)
}
cat("  Saved Figure 4A & 4B\n")

# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 5 – Candidate Biomarker Ranking
# ═══════════════════════════════════════════════════════════════════════════════
cat("Figure 5: Biomarker ranking ...\n")

top10 <- ranking %>% slice(1:10)

fig5 <- ggplot(top10, aes(x = reorder(Gene, Score), y = Score, fill = Disease_Count)) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_fill_gradient(low = "#A6CEE3", high = "#1F78B4",
                      name = "Diseases\nSignificant") +
  theme_classic(base_size = 15) +
  labs(
    title = "Top Candidate Oxidative Stress Biomarkers",
    subtitle = paste0("Scored by (# diseases significant) × (mean |logFC|)"),
    x = "", y = "Biomarker Score"
  ) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

save_publication_plot(fig5,
                      file.path(FIGURE_DIR, "Figure5_Biomarker_Ranking"),
                      width = 8, height = 6)
cat("  Saved Figure 5\n")

# ── Session information ───────────────────────────────────────────────────────
save_session_info()

write_log("06_figures.R complete")
cat("\nAll figures complete. Saved in:", FIGURE_DIR, "\n")
