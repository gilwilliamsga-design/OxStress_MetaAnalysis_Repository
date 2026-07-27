###############################################################################
# io_functions.R
#
# Project : OxStress Meta-analysis
# Author  : Gil Williams Amana
#
# Purpose : Standardised I/O helpers used throughout the pipeline.
#           Source after 00_config.R so LOG_DIR and RESULTS_DIR are defined.
###############################################################################

ensure_directory <- function(path) {
  if (!dir.exists(path))
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

# ── CSV ───────────────────────────────────────────────────────────────────────
save_csv <- function(df, filename) {
  ensure_directory(dirname(filename))
  write.csv(df, filename, row.names = FALSE)
  message("Saved: ", filename)
}

load_csv <- function(filename) {
  if (!file.exists(filename)) stop("Missing file: ", filename)
  read.csv(filename, stringsAsFactors = FALSE)
}

# ── RDS ───────────────────────────────────────────────────────────────────────
save_rds <- function(object, filename) {
  ensure_directory(dirname(filename))
  saveRDS(object, filename)
  message("Saved: ", filename)
}

load_rds <- function(filename) {
  if (!file.exists(filename)) stop("Missing file: ", filename)
  readRDS(filename)
}

# ── Publication figure export (ggplot) ────────────────────────────────────────
save_publication_plot <- function(plot, stem,
                                  width  = FIG_WIDTH,
                                  height = FIG_HEIGHT,
                                  dpi    = FIG_DPI) {
  ensure_directory(dirname(stem))
  ggplot2::ggsave(paste0(stem, ".png"),  plot, width = width, height = height, dpi = dpi)
  ggplot2::ggsave(paste0(stem, ".tiff"), plot, width = width, height = height, dpi = dpi,
                  compression = "lzw")
  ggplot2::ggsave(paste0(stem, ".pdf"),  plot, width = width, height = height)
  message("Figures saved: ", basename(stem), ".{png,tiff,pdf}")
}

# ── Pipeline log ──────────────────────────────────────────────────────────────
write_log <- function(text, logfile = file.path(LOG_DIR, "pipeline.log")) {
  ensure_directory(dirname(logfile))
  cat(format(Sys.time()), ":", text, "\n", file = logfile, append = TRUE)
}

# ── Session info ──────────────────────────────────────────────────────────────
save_session_info <- function(outfile = file.path(RESULTS_DIR, "SessionInfo.txt")) {
  capture.output(sessionInfo(), file = outfile)
  message("Session info written: ", outfile)
}

check_file <- function(filename) {
  if (!file.exists(filename))
    stop("Required file not found:\n  ", filename)
}
