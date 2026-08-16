#!/usr/bin/env Rscript
#
# Combined Seurat object build + RNA QC + filtering + doublet removal
# pipeline -- consolidates scripts 01, 02, 03, 04, 05, 06, 09 into one
# script, covering everything from the raw 10x-format source data to a
# fully filtered, doublet-free RNA population.
#
#   1. Build Seurat object, RNA + ADT assays (script 01)
#   2. QC report figure (script 06) -- computed on the raw (pre-filter)
#      object, since it documents the threshold-justification lines for
#      the filters about to be applied below. Doesn't chain sequentially
#      with the rest -- it's a diagnostic on raw data, not a step that
#      consumes/produces a filtered checkpoint.
#   3. Cell filter (script 02): nFeature_RNA >= 200
#   4. Gene filter (script 03): detected in >= 100 cells
#   5. Mito filter (script 04): percent.mt < 20%
#   6. Quantile trim (script 05): 2-98th percentile, per pool
#   7. Doublet removal (script 09): DoubletFinder + HTO, per-lane, union
#
# Two intermediate checkpoints are still saved, not just the final
# filtered object -- because 02_adt_qc_filter.R reads both of them
# directly (its own QC-report step reads the raw object; its filtering
# step reads the quantile-trimmed-but-pre-doublet-removal population):
#   - 01_rna_seurat_object_raw.rds
#   - 01_rna_seurat_object_quantile_trim.rds
# The other old per-step checkpoints (02/03/04) only existed for
# resuming between separate script invocations, which doesn't apply
# once combined -- see README's GitHub section on results/*.rds size.
# Each step's filter CSV is still written individually either way,
# since those are cheap and valuable on their own.
#
# Usage:
#   conda activate citeseq-pipeline
#   Rscript 01_rna_qc_filter.R

suppressMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
  library(DoubletFinder)
})

set.seed(42)

# Auto-detect the project root from this script's own location (scripts/
# is always one level below it), so paths work regardless of where the
# repo is cloned.
script_args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", script_args[grep("^--file=", script_args)])
project_dir <- dirname(dirname(normalizePath(script_path)))
data_dir    <- file.path(project_dir, "data")
results_dir <- file.path(project_dir, "results")
figures_dir <- results_dir
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

hto_source   <- file.path(project_dir, "reference", "hto_doublet_calls.csv")
min_genes      <- 200
min_cells_gene <- 100
mito_max       <- 20
quantile_lo    <- 0.02
quantile_hi    <- 0.98
doublet_rate   <- 0.075

log_msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), sprintf(...), "\n", sep = "")
}

## =============================================================================
## Step 1/7 (script 01): build Seurat object (RNA + ADT) from raw 10x data
## =============================================================================

log_msg("Step 1/7: Building Seurat object from raw 10x-format data...")

log_msg("Reading RNA matrix...")
rna_counts <- ReadMtx(
  mtx      = file.path(data_dir, "GSM5008737_RNA_3P-matrix.mtx.gz"),
  cells    = file.path(data_dir, "GSM5008737_RNA_3P-barcodes.tsv.gz"),
  features = file.path(data_dir, "GSM5008737_RNA_3P-features.tsv.gz"),
  feature.column = 2
)
log_msg("RNA matrix: %d genes x %d cells", nrow(rna_counts), ncol(rna_counts))

log_msg("Reading ADT matrix...")
adt_counts <- ReadMtx(
  mtx      = file.path(data_dir, "GSM5008738_ADT_3P-matrix.mtx.gz"),
  cells    = file.path(data_dir, "GSM5008738_ADT_3P-barcodes.tsv.gz"),
  features = file.path(data_dir, "GSM5008738_ADT_3P-features.tsv.gz"),
  feature.column = 2
)
log_msg("ADT matrix: %d antibodies x %d cells", nrow(adt_counts), ncol(adt_counts))

stopifnot(identical(colnames(rna_counts), colnames(adt_counts)))

log_msg("Assembling Seurat object...")
obj <- CreateSeuratObject(counts = rna_counts, assay = "RNA", project = "citeseq_pipeline")
obj[["ADT"]] <- CreateAssay5Object(counts = adt_counts)
rm(rna_counts, adt_counts); gc()

lane <- sub("_[ACGT]+$", "", colnames(obj))
obj$lane <- lane
obj$pool <- ifelse(grepl("^E2", lane), "E2_pool", "L_pool")

log_msg("Seurat object built: %d cells, assays: %s", ncol(obj), paste(Assays(obj), collapse = ", "))
log_msg("Pool breakdown: %s", paste(capture.output(print(table(obj$pool))), collapse = " | "))

saveRDS(obj, file.path(results_dir, "01_rna_seurat_object_raw.rds"))
log_msg("Saved results/01_rna_seurat_object_raw.rds")
rm(obj); gc()

## =============================================================================
## Step 2/7 (script 06): RNA QC report figure, on the RAW object
## =============================================================================

log_msg("Step 2/7: RNA QC report figure (raw, pre-filter)...")

# Re-read from disk rather than reusing the in-memory object from step 1 --
# R's copy-on-modify semantics would otherwise deep-copy the full RNA+ADT
# object the moment a slot is touched below, transiently doubling memory
# and risking an OOM kill on a 161k-cell object (this crashed once before
# this re-read-from-disk fix was added).
obj_raw <- readRDS(file.path(results_dir, "01_rna_seurat_object_raw.rds"))
obj_raw[["percent.mt"]] <- PercentageFeatureSet(obj_raw, assay = "RNA", pattern = "^MT-")
ribo_lines <- readLines(file.path(project_dir, "reference", "KEGG_RIBOSOME.txt"))
ribo_genes <- ribo_lines[nzchar(ribo_lines) & !grepl("^>", ribo_lines) & ribo_lines != "KEGG_RIBOSOME"]
ribo_genes_present <- intersect(ribo_genes, rownames(obj_raw[["RNA"]]))
obj_raw[["percent.ribo"]] <- PercentageFeatureSet(obj_raw, assay = "RNA", features = ribo_genes_present)

qc_df <- FetchData(obj_raw, vars = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo", "pool"))

rna_counts_raw <- GetAssayData(obj_raw, assay = "RNA", layer = "counts")
cells_per_gene <- rowSums(rna_counts_raw > 0)
rm(rna_counts_raw, obj_raw); gc()

pool_colors <- c(L_pool = "#4C9A6A", E2_pool = "#E45756")

p_total <- ggplot(qc_df, aes(x = "", y = nCount_RNA)) +
  geom_violin(fill = "#F8766D", trim = TRUE) +
  labs(title = "A. Total Count", x = NULL, y = "Total UMI count") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.border = element_rect(color = "grey40", fill = NA), panel.grid = element_blank())

p_a <- ggplot(qc_df, aes(x = "", y = nFeature_RNA)) +
  geom_violin(fill = "#B279A2", trim = TRUE) +
  geom_hline(yintercept = min_genes, color = "black", linetype = "dashed", linewidth = 0.6) +
  labs(title = "C. Gene Counts per Cell", x = NULL, y = "Number of genes per cell") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.border = element_rect(color = "grey40", fill = NA), panel.grid = element_blank())

p_b <- ggplot(qc_df, aes(x = "", y = percent.mt)) +
  geom_violin(fill = "#B279A2", trim = TRUE) +
  geom_hline(yintercept = mito_max, color = "black", linetype = "dashed", linewidth = 0.6) +
  labs(title = "E. Mitochondrial Fraction", x = NULL, y = "% Mitochondrial reads") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.border = element_rect(color = "grey40", fill = NA), panel.grid = element_blank())

gene_df <- data.frame(cells_per_gene = cells_per_gene)
p_c <- ggplot(gene_df, aes(x = "", y = cells_per_gene)) +
  geom_violin(fill = "#B279A2", trim = TRUE) +
  geom_hline(yintercept = min_cells_gene, color = "black", linetype = "dashed", linewidth = 0.6) +
  scale_y_log10(labels = scales::comma) +
  labs(title = "D. Cell Counts per Gene", x = NULL, y = "Number of cells per gene (log10)") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.border = element_rect(color = "grey40", fill = NA), panel.grid = element_blank())

p_d <- ggplot(qc_df, aes(x = "", y = percent.ribo)) +
  geom_violin(fill = "#B279A2", trim = TRUE) +
  labs(title = "F. Ribosomal RNA Fraction", x = NULL, y = "% Ribosomal reads") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.border = element_rect(color = "grey40", fill = NA), panel.grid = element_blank())

p_e <- ggplot(qc_df, aes(x = nCount_RNA, y = nFeature_RNA)) +
  geom_hex(bins = 60) +
  scale_fill_viridis_c(option = "plasma", trans = "log10", name = "cells") +
  labs(title = "B. Genes vs. Total Counts", x = "Total UMI count", y = "Number of genes") +
  theme_minimal(base_size = 12) +
  theme(panel.border = element_rect(color = "grey40", fill = NA), panel.grid = element_blank())

wt <- wilcox.test(nCount_RNA ~ pool, data = qc_df)
sig_label <- if (wt$p.value < 0.0001) "****" else if (wt$p.value < 0.001) "***" else
             if (wt$p.value < 0.01) "**" else if (wt$p.value < 0.05) "*" else "ns"
p_label <- if (wt$p.value < 2.2e-16) "p < 2.2e-16" else sprintf("p = %.2g", wt$p.value)
bracket_y  <- max(qc_df$nCount_RNA) * 1.25
text_y     <- max(qc_df$nCount_RNA) * 1.32

p_f <- ggplot(qc_df, aes(x = pool, y = nCount_RNA, fill = pool)) +
  geom_violin(trim = TRUE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  annotate("segment", x = 1, xend = 2, y = bracket_y, yend = bracket_y) +
  annotate("segment", x = 1, xend = 1, y = bracket_y, yend = bracket_y * 0.94) +
  annotate("segment", x = 2, xend = 2, y = bracket_y, yend = bracket_y * 0.94) +
  annotate("text", x = 1.5, y = text_y,
           label = sprintf("%s (Wilcoxon, %s)", sig_label, p_label), size = 3.5) +
  scale_fill_manual(values = pool_colors, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  labs(title = "G. Batch Effect Between Pools", x = NULL, y = "Total UMI count") +
  theme_minimal(base_size = 12) +
  theme(panel.border = element_rect(color = "grey40", fill = NA), panel.grid = element_blank())

p_combined <- wrap_plots(p_total, p_e, p_a, p_c, p_b, p_d, p_f, ncol = 2) +
  plot_annotation(title = "RNA", theme = theme(plot.title = element_text(size = 15, face = "bold")))
ggsave(file.path(figures_dir, "01_rna_qc_report_summary.png"), p_combined, width = 12, height = 17, dpi = 150)
log_msg("Saved results/01_rna_qc_report_summary.png")

rm(qc_df, cells_per_gene, gene_df, p_total, p_a, p_b, p_c, p_d, p_e, p_f, p_combined); gc()

## =============================================================================
## Step 3/7 (script 02): cell filter, nFeature_RNA >= 200
## =============================================================================

log_msg("Step 3/7: Cell filter (nFeature_RNA >= %d)...", min_genes)
obj <- readRDS(file.path(results_dir, "01_rna_seurat_object_raw.rds"))
n_before <- ncol(obj)

n_fail <- sum(obj$nFeature_RNA < min_genes)
log_msg("Cells with nFeature_RNA < %d: %d", min_genes, n_fail)
obj <- subset(obj, subset = nFeature_RNA >= min_genes)
n_after <- ncol(obj)
log_msg("Cells after cell filter: %d (removed %d)", n_after, n_before - n_after)

write.csv(data.frame(step = "cell_filter_min_genes_200", cells_before = n_before,
                      cells_after = n_after, cells_removed = n_before - n_after),
          file.path(results_dir, "01_rna_qc_filter_cell.csv"), row.names = FALSE)
log_msg("Saved results/01_rna_qc_filter_cell.csv")

## =============================================================================
## Step 4/7 (script 03): gene filter, detected in >= 100 cells
## =============================================================================

log_msg("Step 4/7: Gene filter (>= %d cells/gene)...", min_cells_gene)
n_genes_before <- nrow(obj[["RNA"]])
n_cells <- ncol(obj)

rna_counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
gene_cell_counts <- Matrix::rowSums(rna_counts > 0)
keep_genes <- names(gene_cell_counts)[gene_cell_counts >= min_cells_gene]
n_genes_removed <- n_genes_before - length(keep_genes)
log_msg("Genes with cell count < %d: %d", min_cells_gene, n_genes_removed)

obj[["RNA"]] <- subset(obj[["RNA"]], features = keep_genes)
n_genes_after <- nrow(obj[["RNA"]])
log_msg("Genes after gene filter: %d (removed %d)", n_genes_after, n_genes_removed)
rm(rna_counts, gene_cell_counts); gc()

write.csv(data.frame(step = "gene_filter_min_cells_100", genes_before = n_genes_before,
                      genes_after = n_genes_after, genes_removed = n_genes_removed, cells = n_cells),
          file.path(results_dir, "01_rna_qc_filter_gene.csv"), row.names = FALSE)
log_msg("Saved results/01_rna_qc_filter_gene.csv")

## =============================================================================
## Step 5/7 (script 04): mito filter, percent.mt < 20%
## =============================================================================

log_msg("Step 5/7: Mito filter (percent.mt < %d%%)...", mito_max)
n_before <- ncol(obj)
obj$percent.mt <- PercentageFeatureSet(obj, pattern = "^MT-", assay = "RNA")
n_fail <- sum(obj$percent.mt >= mito_max)
log_msg("Cells with percent.mt >= %d%%: %d", mito_max, n_fail)

obj <- subset(obj, subset = percent.mt < mito_max)
n_after <- ncol(obj)
log_msg("Cells after mito filter: %d (removed %d)", n_after, n_before - n_after)

write.csv(data.frame(step = "mito_filter_lt_20pct", cells_before = n_before,
                      cells_after = n_after, cells_removed = n_before - n_after),
          file.path(results_dir, "01_rna_qc_filter_mito.csv"), row.names = FALSE)
log_msg("Saved results/01_rna_qc_filter_mito.csv")

## =============================================================================
## Step 6/7 (script 05): quantile trim, 2-98th percentile, per pool
## =============================================================================

log_msg("Step 6/7: Quantile trim (%.0f-%.0fth percentile, per pool)...", quantile_lo * 100, quantile_hi * 100)
n_before <- ncol(obj)
pools <- sort(unique(obj$pool))

keep <- rep(FALSE, ncol(obj))
names(keep) <- colnames(obj)

for (pl in pools) {
  cells_in_pool <- colnames(obj)[obj$pool == pl]
  count_vals <- obj$nCount_RNA[cells_in_pool]
  feat_vals  <- obj$nFeature_RNA[cells_in_pool]

  q_count <- quantile(count_vals, c(quantile_lo, quantile_hi))
  q_feat  <- quantile(feat_vals, c(quantile_lo, quantile_hi))
  log_msg("Pool %s: nCount_RNA [%.0f%%,%.0f%%]=[%.0f,%.0f]  nFeature_RNA [%.0f%%,%.0f%%]=[%.0f,%.0f]",
          pl, quantile_lo * 100, quantile_hi * 100, q_count[1], q_count[2],
          quantile_lo * 100, quantile_hi * 100, q_feat[1], q_feat[2])

  pool_keep <- count_vals >= q_count[1] & count_vals <= q_count[2] &
               feat_vals  >= q_feat[1]  & feat_vals  <= q_feat[2]
  keep[cells_in_pool] <- pool_keep
}

obj <- subset(obj, cells = colnames(obj)[keep])
n_after <- ncol(obj)
log_msg("Cells after quantile trim: %d (removed %d)", n_after, n_before - n_after)

write.csv(data.frame(step = "quantile_trim_2_98pct_per_pool", cells_before = n_before,
                      cells_after = n_after, cells_removed = n_before - n_after),
          file.path(results_dir, "01_rna_qc_filter_quantile.csv"), row.names = FALSE)
log_msg("Saved results/01_rna_qc_filter_quantile.csv")

# Saved here, not just as a final-object-only checkpoint, because
# 02_adt_qc_filter.R depends on this exact quantile-trimmed-but-
# pre-doublet-removal population (ADT filtering is deliberately run
# before RNA-side doublet removal) -- without this file, that script
# has no valid input.
saveRDS(obj, file.path(results_dir, "01_rna_seurat_object_quantile_trim.rds"))
log_msg("Saved results/01_rna_seurat_object_quantile_trim.rds")

## =============================================================================
## Step 7/7 (script 09): doublet removal, DoubletFinder + HTO, per-lane, union
## =============================================================================

log_msg("Step 7/7: Doublet removal (DoubletFinder + HTO, per-lane, union)...")
n_before <- ncol(obj)
lanes <- sort(unique(obj$lane))
log_msg("Detected %d lanes: %s", length(lanes), paste(lanes, collapse = ", "))

classifications <- setNames(vector("list", length(lanes)), lanes)

for (ln in lanes) {
  cells_in_lane <- colnames(obj)[obj$lane == ln]
  log_msg("Lane %s: %d cells -- preprocessing...", ln, length(cells_in_lane))
  obj_lane <- subset(obj, cells = cells_in_lane)

  obj_lane <- NormalizeData(obj_lane, verbose = FALSE)
  obj_lane <- FindVariableFeatures(obj_lane, nfeatures = 2000, verbose = FALSE)
  obj_lane <- ScaleData(obj_lane, verbose = FALSE)
  obj_lane <- RunPCA(obj_lane, npcs = 10, seed.use = 42, verbose = FALSE)

  optimal_pK <- 0.01
  nExp_lane <- round(doublet_rate * ncol(obj_lane))
  log_msg("Lane %s: running doubletFinder (pK=%s, nExp=%d)...", ln, optimal_pK, nExp_lane)
  # doubletFinder()'s only randomness is two plain sample() calls internally
  # (confirmed by inspecting its source) -- it takes no seed argument itself
  # and otherwise just consumes whatever position the single top-of-script
  # set.seed(42) RNG stream happens to be at by the time this lane's turn
  # comes up, which depends on exactly how many prior random draws happened
  # in every earlier lane. Reseed fresh right here so each lane's doublet
  # calls are reproducible on their own, independent of loop position.
  set.seed(42)
  obj_lane <- doubletFinder(obj_lane, PCs = 1:10, pN = 0.25, pK = optimal_pK,
                             nExp = nExp_lane, sct = FALSE)

  doublet_col <- grep("^DF.classifications", colnames(obj_lane@meta.data), value = TRUE)
  lane_calls <- setNames(obj_lane@meta.data[[doublet_col]], colnames(obj_lane))
  classifications[[ln]] <- lane_calls
  log_msg("Lane %s classifications: %s", ln,
          paste(names(table(lane_calls)), table(lane_calls), sep = "=", collapse = ", "))

  rm(obj_lane, lane_calls); gc()
}

all_calls <- unlist(classifications, use.names = TRUE)
names(all_calls) <- sub("^[^.]+\\.", "", names(all_calls))
obj$doublet_finder_call <- all_calls[colnames(obj)]
log_msg("DoubletFinder calls across all lanes:")
print(table(obj$doublet_finder_call))

log_msg("Reading HTO-based doublet calls (%s)...", hto_source)
hto_df <- read.csv(hto_source, stringsAsFactors = FALSE)
hto_lookup <- setNames(hto_df$hto_doublet_call, hto_df$barcode)
obj$hto_doublet_call <- unname(hto_lookup[colnames(obj)])
log_msg("HTO-based calls:")
print(table(obj$hto_doublet_call, useNA = "ifany"))

xt <- table(DoubletFinder = obj$doublet_finder_call, HTO = obj$hto_doublet_call, useNA = "ifany")
log_msg("Cross-tab: DoubletFinder vs. HTO-based")
print(xt)

n_both <- sum(obj$doublet_finder_call == "Doublet" & obj$hto_doublet_call == "Doublet", na.rm = TRUE)
n_df_only <- sum(obj$doublet_finder_call == "Doublet" & obj$hto_doublet_call == "Singlet", na.rm = TRUE)
n_hto_only <- sum(obj$doublet_finder_call == "Singlet" & obj$hto_doublet_call == "Doublet", na.rm = TRUE)
log_msg("Flagged by both: %d | DoubletFinder only: %d | HTO only: %d", n_both, n_df_only, n_hto_only)

obj$is_doublet <- (obj$doublet_finder_call == "Doublet") |
                   (!is.na(obj$hto_doublet_call) & obj$hto_doublet_call == "Doublet")
n_removed <- sum(obj$is_doublet)
obj <- subset(obj, subset = !is_doublet)
n_after <- ncol(obj)
log_msg("Cells after doublet removal: %d (removed %d)", n_after, n_before - n_after)

write.csv(data.frame(step = c("doubletfinder_only", "hto_only", "both", "union_removed"),
                      n_cells = c(n_df_only, n_hto_only, n_both, n_removed)),
          file.path(results_dir, "01_rna_doublet_summary.csv"), row.names = FALSE)
log_msg("Saved results/01_rna_doublet_summary.csv")

## =============================================================================
## Final checkpoint -- replaces what was previously 01_rna_seurat_object_doublet_filtered.rds
## =============================================================================

saveRDS(obj, file.path(results_dir, "01_rna_seurat_object_doublet_filtered.rds"))
log_msg("Saved results/01_rna_seurat_object_doublet_filtered.rds")
log_msg("Done: %d cells, %d genes (RNA) -- fully QC'd, filtered, doublet-free.", ncol(obj), nrow(obj[["RNA"]]))
