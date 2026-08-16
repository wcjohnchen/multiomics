#!/usr/bin/env Rscript
#
# Regenerates report.html from report_template.html, using this run's
# actual results/ files -- the template keeps the exact page formatting
# fixed (same design, layout, prose, methodology explanations); only the
# 6 figures and 15 numeric/statistical values are pulled fresh each run.
#
# Every value is read from a results/ CSV or a results/ PNG -- none of it
# re-derives anything by loading the large .rds checkpoints, so this runs
# in seconds, not minutes.
#
# Usage:
#   conda activate citeseq-pipeline
#   Rscript 06_generate_report.R

suppressMessages({
  library(base64enc)
})

script_args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", script_args[grep("^--file=", script_args)])
project_dir <- dirname(dirname(normalizePath(script_path)))
results_dir <- file.path(project_dir, "results")
template_path <- file.path(project_dir, "scripts", "report_template.html")
output_path <- file.path(project_dir, "report.html")

log_msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), sprintf(...), "\n", sep = "")
}

fmt_int <- function(x) format(round(x), big.mark = ",", scientific = FALSE, trim = TRUE)

fmt_p <- function(p) {
  if (p < 2.2e-16) "p &lt; 2.2e-16" else sprintf("p = %.2g", p)
}

## =============================================================================
## Step 1/3: read every results/ file this template needs
## =============================================================================

log_msg("Step 1/3: Reading results/ files...")

rna_cell_filter   <- read.csv(file.path(results_dir, "01_rna_qc_filter_cell.csv"))
rna_gene_filter   <- read.csv(file.path(results_dir, "01_rna_qc_filter_gene.csv"))
rna_doublets      <- read.csv(file.path(results_dir, "01_rna_doublet_summary.csv"))
rna_qc_stats      <- read.csv(file.path(results_dir, "01_rna_qc_stats.csv"))
adt_qc_filter     <- read.csv(file.path(results_dir, "02_adt_qc_filter.csv"))
adt_qc_stats      <- read.csv(file.path(results_dir, "02_adt_qc_stats.csv"))
rna_broad         <- read.csv(file.path(results_dir, "03_rna_broad_annotation.csv"))
adt_broad         <- read.csv(file.path(results_dir, "04_adt_broad_annotation.csv"))
wnn_detailed      <- read.csv(file.path(results_dir, "05_wnn_detailed_annotation.csv"))
wnn_broad         <- read.csv(file.path(results_dir, "05_wnn_broad_annotation.csv"))

## =============================================================================
## Step 2/3: compute every template token
## =============================================================================

log_msg("Step 2/3: Computing template values...")

raw_cells      <- rna_cell_filter$cells_before[1]
rna_cells_qc   <- adt_qc_filter$cells_before[1]   # doublet-filtered checkpoint's cell count
genes_retained <- rna_gene_filter$genes_after[1]
qc_pass_pct    <- sprintf("%.1f", 100 * rna_cells_qc / raw_cells)

rna_wilcoxon_p    <- fmt_p(rna_qc_stats$value[rna_qc_stats$metric == "wilcoxon_p_value"])
rna_batch_ratio   <- sprintf("%.2f", rna_qc_stats$value[rna_qc_stats$metric == "median_ratio_E2_pool_over_L_pool"])
adt_wilcoxon_p    <- fmt_p(adt_qc_stats$value[adt_qc_stats$metric == "wilcoxon_p_value"])

adt_ab_total    <- adt_qc_filter$antibodies_before[1]
adt_ab_retained <- adt_qc_filter$antibodies_after[nrow(adt_qc_filter)]

doublets_removed   <- rna_doublets$n_cells[rna_doublets$step == "union_removed"]
doublet_both_count <- rna_doublets$n_cells[rna_doublets$step == "both"]

wnn_clusters      <- nrow(wnn_detailed)
detailed_lineages <- length(unique(wnn_detailed$broad_lineage))
broad_lineages    <- length(unique(wnn_broad$l1_lineage))

rna_broad_categories <- length(unique(rna_broad$broad_lineage))
adt_broad_categories <- length(unique(adt_broad$broad_lineage))

tokens <- list(
  RNA_CELLS_QC          = fmt_int(rna_cells_qc),
  RAW_CELLS             = fmt_int(raw_cells),
  GENES_RETAINED        = fmt_int(genes_retained),
  QC_PASS_PCT           = qc_pass_pct,
  RNA_BATCH_RATIO       = rna_batch_ratio,
  RNA_WILCOXON_P        = rna_wilcoxon_p,
  ADT_AB_RETAINED       = fmt_int(adt_ab_retained),
  ADT_AB_TOTAL          = fmt_int(adt_ab_total),
  ADT_WILCOXON_P        = adt_wilcoxon_p,
  DOUBLETS_REMOVED      = fmt_int(doublets_removed),
  DOUBLET_BOTH_COUNT    = fmt_int(doublet_both_count),
  WNN_CLUSTERS          = fmt_int(wnn_clusters),
  DETAILED_LINEAGES     = fmt_int(detailed_lineages),
  BROAD_LINEAGES        = fmt_int(broad_lineages),
  RNA_BROAD_CATEGORIES  = fmt_int(rna_broad_categories),
  ADT_BROAD_CATEGORIES  = fmt_int(adt_broad_categories)
)

for (name in names(tokens)) {
  log_msg("  %s = %s", name, tokens[[name]])
}

image_files <- c(
  IMG_RNA_QC       = "01_rna_qc_report_summary.png",
  IMG_ADT_QC       = "02_adt_qc_report_summary.png",
  IMG_RNA_UMAP     = "03_rna_umap_broad_labels.png",
  IMG_ADT_UMAP     = "04_adt_umap_broad_labels.png",
  IMG_WNN_DETAILED = "05_wnn_umap_detailed_labels.png",
  IMG_WNN_BROAD    = "05_wnn_umap_broad_labels.png"
)

for (name in names(image_files)) {
  path <- file.path(results_dir, image_files[[name]])
  stopifnot(file.exists(path))
  tokens[[name]] <- base64enc::base64encode(path)
}
log_msg("Encoded %d figures", length(image_files))

## =============================================================================
## Step 3/3: fill the template and write report.html
## =============================================================================

log_msg("Step 3/3: Filling template and writing report.html...")

html <- paste(readLines(template_path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")

for (name in names(tokens)) {
  html <- gsub(paste0("{{", name, "}}"), tokens[[name]], html, fixed = TRUE)
}

remaining <- regmatches(html, gregexpr("\\{\\{[A-Z_]+\\}\\}", html))[[1]]
if (length(remaining) > 0) {
  stop("Unfilled template token(s) remain: ", paste(unique(remaining), collapse = ", "))
}

writeLines(html, output_path, useBytes = TRUE)
log_msg("Saved report.html")
