#!/usr/bin/env Rscript
#
# Combined ADT QC + filtering pipeline -- consolidates scripts 07 and 08.
#
#   1. QC report figure (script 08) -- computed FIRST, on the raw
#      (pre-filter) object, since it documents the threshold-justification
#      lines for the filters about to be applied below. Doesn't chain
#      sequentially with the filter step -- it's a diagnostic on raw data.
#   2. Antibody filter + cell filter (script 07): antibody detected in
#      >=100 cells; cell has nCount_ADT>=20 AND nFeature_ADT>=20.
#
# Runs on top of the RNA-side quantile-trimmed checkpoint (05). Note: the
# resulting checkpoint here is -- same as before combining -- not read by
# any later step in this rebuild (script 09 branches from 05, script 18
# branches from 16, both skipping past this one). Kept for its QC
# figure/filter record, not because downstream steps depend on it.
#
# Usage:
#   conda activate citeseq-pipeline
#   Rscript 02_adt_qc_filter.R

suppressMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
})

# Auto-detect the project root from this script's own location (scripts/
# is always one level below it), so paths work regardless of where the
# repo is cloned.
script_args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", script_args[grep("^--file=", script_args)])
project_dir <- dirname(dirname(normalizePath(script_path)))
results_dir <- file.path(project_dir, "results")
figures_dir <- results_dir
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

min_adt_feature        <- 20
min_adt_count          <- 20
min_cells_per_antibody <- 100

log_msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), sprintf(...), "\n", sep = "")
}

## =============================================================================
## Step 1 (script 08): ADT QC report figure, on the RAW object
## =============================================================================

log_msg("Step 1/2: ADT QC report figure (raw, pre-filter)...")
obj_raw <- readRDS(file.path(results_dir, "01_rna_seurat_object_raw.rds"))
log_msg("Raw: %d cells, %d antibodies (ADT)", ncol(obj_raw), nrow(obj_raw[["ADT"]]))

isotype_controls <- c("Rat-IgG1-1", "Rat-IgG1-2", "Rat-IgG2b", "Rag-IgG2c")
obj_raw[["percent.isotype"]] <- PercentageFeatureSet(obj_raw, assay = "ADT", features = isotype_controls)

qc_df <- FetchData(obj_raw, vars = c("nCount_ADT", "nFeature_ADT", "percent.isotype", "pool"))

adt_counts_raw <- GetAssayData(obj_raw, assay = "ADT", layer = "counts")
cells_per_antibody_raw <- rowSums(adt_counts_raw > 0)
rm(adt_counts_raw, obj_raw); gc()

pool_colors <- c(L_pool = "#4C9A6A", E2_pool = "#E45756")

p_a <- ggplot(qc_df, aes(x = nCount_ADT, y = nFeature_ADT)) +
  geom_hex(bins = 60) +
  scale_fill_viridis_c(option = "plasma", trans = "log10", name = "cells") +
  labs(title = "B. Antibodies vs. Total Counts", x = "Total ADT UMI count", y = "Number of antibodies") +
  theme_minimal(base_size = 12) +
  theme(panel.border = element_rect(color = "grey40", fill = NA), panel.grid = element_blank())

p_b <- ggplot(qc_df, aes(x = "", y = nFeature_ADT)) +
  geom_violin(fill = "#B279A2", trim = TRUE) +
  geom_hline(yintercept = min_adt_feature, color = "black", linetype = "dashed", linewidth = 0.6) +
  labs(title = "C. Antibody Counts per Cell", x = NULL, y = "Number of antibodies per cell") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.border = element_rect(color = "grey40", fill = NA), panel.grid = element_blank())

antibody_df <- data.frame(cells_per_antibody = cells_per_antibody_raw)
p_c <- ggplot(antibody_df, aes(x = "", y = cells_per_antibody)) +
  geom_violin(fill = "#B279A2", trim = TRUE) +
  geom_hline(yintercept = min_cells_per_antibody, color = "black", linetype = "dashed", linewidth = 0.6) +
  scale_y_log10(labels = scales::comma) +
  labs(title = "D. Cell Counts per Antibody", x = NULL, y = "Number of cells per antibody (log10)") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.border = element_rect(color = "grey40", fill = NA), panel.grid = element_blank())

p_d <- ggplot(qc_df, aes(x = "", y = nCount_ADT)) +
  geom_violin(fill = "#F8766D", trim = TRUE) +
  labs(title = "A. Total Count", x = NULL, y = "Total ADT UMI count") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.border = element_rect(color = "grey40", fill = NA), panel.grid = element_blank())

p_e <- ggplot(qc_df, aes(x = "", y = percent.isotype)) +
  geom_violin(fill = "#B279A2", trim = TRUE) +
  labs(title = "E. Isotype Control Fraction", x = NULL, y = "% Isotype control reads") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.border = element_rect(color = "grey40", fill = NA), panel.grid = element_blank())

wt <- wilcox.test(nCount_ADT ~ pool, data = qc_df)
sig_label <- if (wt$p.value < 0.0001) "****" else if (wt$p.value < 0.001) "***" else
             if (wt$p.value < 0.01) "**" else if (wt$p.value < 0.05) "*" else "ns"
p_label <- if (wt$p.value < 2.2e-16) "p < 2.2e-16" else sprintf("p = %.2g", wt$p.value)
bracket_y  <- max(qc_df$nCount_ADT) * 1.25
text_y     <- max(qc_df$nCount_ADT) * 1.32

p_f <- ggplot(qc_df, aes(x = pool, y = nCount_ADT, fill = pool)) +
  geom_violin(trim = TRUE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  annotate("segment", x = 1, xend = 2, y = bracket_y, yend = bracket_y) +
  annotate("segment", x = 1, xend = 1, y = bracket_y, yend = bracket_y * 0.94) +
  annotate("segment", x = 2, xend = 2, y = bracket_y, yend = bracket_y * 0.94) +
  annotate("text", x = 1.5, y = text_y,
           label = sprintf("%s (Wilcoxon, %s)", sig_label, p_label), size = 3.5) +
  scale_fill_manual(values = pool_colors, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  labs(title = "F. Batch Effect Between Pools", x = NULL, y = "Total ADT UMI count") +
  theme_minimal(base_size = 12) +
  theme(panel.border = element_rect(color = "grey40", fill = NA), panel.grid = element_blank())

p_combined <- wrap_plots(p_d, p_a, p_b, p_c, p_e, p_f, ncol = 2) +
  plot_annotation(title = "ADT", theme = theme(plot.title = element_text(size = 15, face = "bold")))
ggsave(file.path(figures_dir, "02_adt_qc_report_summary.png"), p_combined, width = 12, height = 13, dpi = 150)
log_msg("Saved results/02_adt_qc_report_summary.png")
log_msg("Wilcoxon test nCount_ADT ~ pool: p = %.4g", wt$p.value)

rm(qc_df, cells_per_antibody_raw, antibody_df, p_a, p_b, p_c, p_d, p_e, p_f, p_combined); gc()

## =============================================================================
## Step 2 (script 07): ADT antibody filter + cell filter
## =============================================================================

log_msg("Step 2/2: ADT filtering (antibody + cell)...")
obj <- readRDS(file.path(results_dir, "01_rna_seurat_object_quantile_trim.rds"))
n_cells_before <- ncol(obj)
n_ab_before <- nrow(obj[["ADT"]])
log_msg("Before ADT filtering: %d cells, %d antibodies", n_cells_before, n_ab_before)

adt_counts <- GetAssayData(obj, assay = "ADT", layer = "counts")
cells_per_antibody <- Matrix::rowSums(adt_counts > 0)
antibodies_keep <- names(cells_per_antibody)[cells_per_antibody >= min_cells_per_antibody]
n_ab_removed <- n_ab_before - length(antibodies_keep)
log_msg("Keeping %d / %d antibodies detected in >= %d cells (%d removed)",
        length(antibodies_keep), n_ab_before, min_cells_per_antibody, n_ab_removed)

obj[["ADT"]] <- subset(obj[["ADT"]], features = antibodies_keep)

adt_counts <- GetAssayData(obj, assay = "ADT", layer = "counts")
obj$nCount_ADT   <- Matrix::colSums(adt_counts)
obj$nFeature_ADT <- Matrix::colSums(adt_counts > 0)
rm(adt_counts, cells_per_antibody); gc()

log_msg("Recomputed nCount_ADT min=%d, nFeature_ADT min=%d",
        min(obj$nCount_ADT), min(obj$nFeature_ADT))
n_fail <- sum(obj$nCount_ADT < min_adt_count | obj$nFeature_ADT < min_adt_feature)
log_msg("Cells failing nCount_ADT>=%d or nFeature_ADT>=%d: %d", min_adt_count, min_adt_feature, n_fail)

obj <- subset(obj, subset = nCount_ADT >= min_adt_count & nFeature_ADT >= min_adt_feature)
n_cells_after <- ncol(obj)
log_msg("Cells after ADT cell filter: %d (removed %d)", n_cells_after, n_cells_before - n_cells_after)

qc_filter <- data.frame(
  step = c("antibody_filter_min_cells_100", "cell_filter_nCount_nFeature_ge_20"),
  cells_before = c(n_cells_before, n_cells_before),
  cells_after  = c(n_cells_before, n_cells_after),
  antibodies_before = c(n_ab_before, length(antibodies_keep)),
  antibodies_after  = c(length(antibodies_keep), length(antibodies_keep))
)
write.csv(qc_filter, file.path(results_dir, "02_adt_qc_filter.csv"), row.names = FALSE)
log_msg("Saved results/02_adt_qc_filter.csv")

saveRDS(obj, file.path(results_dir, "02_adt_seurat_object_filtered.rds"))
log_msg("Saved results/02_adt_seurat_object_filtered.rds")
log_msg("Done: %d cells, %d antibodies (ADT).", ncol(obj), nrow(obj[["ADT"]]))
