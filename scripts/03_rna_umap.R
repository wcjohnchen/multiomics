#!/usr/bin/env Rscript
#
# Combined RNA normalization -> UMAP -> clustering -> annotation pipeline
# -- consolidates scripts 10, 11, 12, 13, 14, 15, 16, 17. Covers
# everything from the doublet-filtered RNA population to the final
# labeled RNA-only UMAP plot.
#
#   1. LogNormalize (script 10)
#   2. HVG selection, top 2000 (script 11)
#   3. Scale + PCA, 30 PCs (script 12)
#   4. UMAP on pca (script 13)
#   5. FindNeighbors + FindClusters, 24 clusters (script 14)
#   6. FindAllMarkers, RNA restricted to HVGs (script 15)
#   7. Manual broad-lineage annotation (script 16), including the two
#      corrections made after visually cross-checking against the UMAP
#      embedding (cluster 15 -> Monocyte, cluster 22 -> ILC -- see README)
#   8. Labeled RNA UMAP plot (script 17)
#
# Only the FINAL object is saved as a checkpoint here (the old per-step
# checkpoints from 10-16 only existed for resuming between separate
# script invocations -- see README's GitHub section on results/*.rds
# size). Each step's own CSV/figure output is still written individually.
#
# Usage:
#   conda activate citeseq-pipeline
#   Rscript 03_rna_umap.R

suppressMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
})

set.seed(42)

# Auto-detect the project root from this script's own location (scripts/
# is always one level below it), so paths work regardless of where the
# repo is cloned.
script_args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", script_args[grep("^--file=", script_args)])
project_dir <- dirname(dirname(normalizePath(script_path)))
results_dir <- file.path(project_dir, "results")
figures_dir <- results_dir
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

rna_npcs <- 30
cluster_resolution <- 0.5

log_msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), sprintf(...), "\n", sep = "")
}

log_msg("Loading doublet-filtered Seurat object...")
obj <- readRDS(file.path(results_dir, "01_rna_seurat_object_doublet_filtered.rds"))
log_msg("%d cells, %d genes (RNA)", ncol(obj), nrow(obj[["RNA"]]))

## =============================================================================
## Step 1/8 (script 10): LogNormalize
## =============================================================================

log_msg("Step 1/8: Normalizing RNA (LogNormalize, scale.factor=10000)...")
DefaultAssay(obj) <- "RNA"
obj <- NormalizeData(obj, assay = "RNA", normalization.method = "LogNormalize",
                      scale.factor = 10000, verbose = FALSE)

raw_rna  <- GetAssayData(obj, assay = "RNA", layer = "counts")[, 1]
norm_rna <- GetAssayData(obj, assay = "RNA", layer = "data")[, 1]
example_gene <- names(raw_rna)[raw_rna > 0][1]
log_msg("Sanity check (cell 1, gene %s): raw count=%d, LogNormalize value=%.4f",
        example_gene, raw_rna[example_gene], norm_rna[example_gene])
rm(raw_rna, norm_rna); gc()

## =============================================================================
## Step 2/8 (script 11): HVG selection, top 2000
## =============================================================================

log_msg("Step 2/8: Finding variable features (RNA, vst, top 2000)...")
obj <- FindVariableFeatures(obj, assay = "RNA", selection.method = "vst",
                             nfeatures = 2000, verbose = FALSE)
hvgs <- VariableFeatures(obj, assay = "RNA")
log_msg("Selected %d variable genes. Top 10: %s", length(hvgs), paste(head(hvgs, 10), collapse = ", "))

## =============================================================================
## Step 3/8 (script 12): Scale + PCA, 30 PCs
## =============================================================================

log_msg("Step 3/8: Scaling data (RNA, HVGs only) + PCA (seed=42)...")
obj <- ScaleData(obj, assay = "RNA", verbose = FALSE)
obj <- RunPCA(obj, assay = "RNA", npcs = rna_npcs, reduction.name = "pca",
              seed.use = 42, verbose = FALSE)
obj[["RNA"]]$scale.data <- NULL
gc()

pca_var <- obj[["pca"]]@stdev^2
pct_var <- round(100 * pca_var / sum(pca_var), 2)
log_msg("PC1-5 %% variance explained: %s", paste(pct_var[1:5], collapse = ", "))

## =============================================================================
## Step 4/8 (script 13): UMAP on pca
## =============================================================================

log_msg("Step 4/8: Running UMAP (RNA pca, dims 1:30, seed=42)...")
obj <- RunUMAP(obj, reduction = "pca", dims = 1:30, reduction.name = "umap.rna",
               reduction.key = "rnaUMAP_", seed.use = 42, verbose = FALSE)

## =============================================================================
## Step 5/8 (script 14): FindNeighbors + FindClusters
## =============================================================================

log_msg("Step 5/8: FindNeighbors + FindClusters (algorithm=3/SLM, resolution=%.2f, seed=42)...",
        cluster_resolution)
obj <- FindNeighbors(obj, reduction = "pca", dims = 1:30, verbose = FALSE)
obj <- FindClusters(obj, algorithm = 3, resolution = cluster_resolution,
                     random.seed = 42, verbose = FALSE)
n_clusters <- length(unique(obj$seurat_clusters))
log_msg("Found %d clusters", n_clusters)
print(table(obj$seurat_clusters))

## =============================================================================
## Step 6/8 (script 15): FindAllMarkers, RNA restricted to HVGs
## =============================================================================

log_msg("Step 6/8: Running FindAllMarkers (RNA, restricted to 2000 HVGs)...")
markers <- FindAllMarkers(obj, assay = "RNA", features = VariableFeatures(obj, assay = "RNA"),
                           only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,
                           max.cells.per.ident = 500, verbose = FALSE)
log_msg("Found %d marker rows across %d clusters", nrow(markers), length(unique(markers$cluster)))
write.csv(markers, file.path(results_dir, "03_rna_cluster_markers.csv"), row.names = FALSE)
log_msg("Saved results/03_rna_cluster_markers.csv")

## =============================================================================
## Step 7/8 (script 16): manual broad-lineage annotation
## =============================================================================

log_msg("Step 7/8: Manual broad-lineage annotation...")

# cluster -> broad_lineage, based on top FindAllMarkers hits (step 6 above).
# Includes the two corrections made after visually cross-checking against
# the UMAP embedding (see step 8 plot / README): cluster 15's markers
# (TNF/CCL3/CCL3L1/CCL4L2/NFKBIA) also read as inflammatory monocyte/
# macrophage, and it embeds inside the Monocyte island, not near any T
# cell cluster -> Monocyte, not "T cell (activated)". Cluster 22's markers
# (KIT/GATA3/IL1R1/SPINK2) plus its position within the CD8 T/NK lymphoid
# branch (far from the true HSPC cluster 21, which sits isolated) fit ILC
# better -> ILC, not "HSPC/Basophil".
broad_lineage_map <- c(
  "0"="CD4 T", "1"="Monocyte", "2"="NK", "3"="CD4 T", "4"="CD8 T",
  "5"="CD4 T", "6"="Monocyte", "7"="B cell", "8"="Monocyte", "9"="CD8 T",
  "10"="Monocyte", "11"="B cell", "12"="Monocyte", "13"="CD8 T", "14"="CD8 T",
  "15"="Monocyte", "16"="gdT/NK", "17"="MAIT", "18"="Monocyte", "19"="DC",
  "20"="Proliferating", "21"="HSPC", "22"="ILC", "23"="pDC"
)

confidence_map <- c(
  "0"="high","1"="high","2"="high","3"="high","4"="high","5"="high",
  "6"="high","7"="high","8"="medium","9"="high","10"="high","11"="high",
  "12"="high","13"="high","14"="high","15"="high","16"="low","17"="high",
  "18"="low","19"="high","20"="high","21"="high","22"="medium","23"="high"
)

key_markers_map <- c(
  "0"="TNFRSF4, CD40LG, IL7R, MAF, GATA3", "1"="S100A8/9/12, LYZ, RETN",
  "2"="NCR1, KLRF1, KLRC1, CD160", "3"="LEF1, CCR7, TCF7, MAL",
  "4"="CD8A, CD8B, GZMH, ZNF683", "5"="CCR7, LEF1, TRABD2A",
  "6"="IL6, IL1A, IL1B, CCL20", "7"="IGHD, TCL1A, FCER2, IGHM",
  "8"="HLA-DQA2, EGR2/3, IFI30", "9"="CD8B, LEF1, CCR7, THEMIS",
  "10"="SIGLEC1, MARCO, IFI44L", "11"="TNFRSF13B, IGHA2, IGHG1-3",
  "12"="CDKN1C, C1QA, MS4A4A", "13"="CD8B, CCR7, LEF1, THEMIS",
  "14"="GZMK, CD8A, CD8B, CCL5", "15"="TNF, CCL3/4, NFKBIA",
  "16"="TRDC, KLRC2, KLRF1, NCAM1", "17"="SLC4A10, KLRB1, TRGC1",
  "18"="TMEM176B, VSTM1, SIGLEC14", "19"="FCER1A, CD1C, CLEC10A, FLT3",
  "20"="MKI67, TOP2A, TYMS, BIRC5", "21"="CD34, PROM1, AVP",
  "22"="KIT, SPINK2, IL1R1, GATA3", "23"="CLEC4C, LILRA4, DNASE1L3"
)

clusters_present <- as.character(sort(unique(as.integer(as.character(obj$seurat_clusters)))))

# This hand-typed map is keyed by cluster *number* from one specific
# reference run (see the big comment above) -- it isn't guaranteed to cover
# every cluster number a differently-clustered run produces, even with the
# same seed and renv-pinned package versions (confirmed: a real rerun once
# produced one extra cluster). Rather than halting the whole pipeline on a
# mismatch, label any unmapped cluster "Unidentified" and continue -- this
# extends the map itself so every downstream lookup below just works
# unchanged.
unmapped_clusters <- setdiff(clusters_present, names(broad_lineage_map))
if (length(unmapped_clusters) > 0) {
  log_msg("WARNING: cluster(s) %s not in the hand-typed broad_lineage_map (likely",
          paste(unmapped_clusters, collapse = ", "))
  log_msg("  non-deterministic clustering vs. the reference run this map was built for)")
  log_msg("  -- labeling as 'Unidentified' rather than halting.")
  for (cl in unmapped_clusters) {
    broad_lineage_map[[cl]] <- "Unidentified"
    confidence_map[[cl]] <- "none"
    key_markers_map[[cl]] <- NA_character_
  }
}

# Opposite direction of the same problem: a cluster *number* the map has an
# entry for might not exist at all in this run's clustering (e.g. the map
# covers 0-23 but only 23 clusters formed). Left in, `table(obj$seurat_
# clusters)[names(broad_lineage_map)]` below would silently produce an NA
# n_cells row for it. Drop those stale entries so the output only ever
# reflects clusters that are actually present.
stale_clusters <- setdiff(names(broad_lineage_map), clusters_present)
if (length(stale_clusters) > 0) {
  log_msg("NOTE: cluster(s) %s from broad_lineage_map don't exist in this run's",
          paste(stale_clusters, collapse = ", "))
  log_msg("  clustering -- dropping them from the output rather than leaving stale rows.")
  broad_lineage_map <- broad_lineage_map[setdiff(names(broad_lineage_map), stale_clusters)]
  confidence_map <- confidence_map[setdiff(names(confidence_map), stale_clusters)]
  key_markers_map <- key_markers_map[setdiff(names(key_markers_map), stale_clusters)]
}

obj$broad_lineage <- unname(broad_lineage_map[as.character(obj$seurat_clusters)])
obj$annotation_confidence <- unname(confidence_map[as.character(obj$seurat_clusters)])

log_msg("Broad lineage breakdown:")
print(sort(table(obj$broad_lineage), decreasing = TRUE))

annotation_table <- data.frame(
  cluster = names(broad_lineage_map),
  n_cells = as.integer(table(obj$seurat_clusters)[names(broad_lineage_map)]),
  broad_lineage = broad_lineage_map,
  confidence = confidence_map,
  key_markers = key_markers_map,
  row.names = NULL
)
annotation_table <- annotation_table[order(as.integer(annotation_table$cluster)), ]
write.csv(annotation_table, file.path(results_dir, "03_rna_broad_annotation.csv"), row.names = FALSE)
log_msg("Saved results/03_rna_broad_annotation.csv")

## =============================================================================
## Step 8/8 (script 17): labeled RNA UMAP plot
## =============================================================================

log_msg("Step 8/8: Building labeled RNA UMAP plot...")

embed <- Embeddings(obj, "umap.rna")
df <- data.frame(UMAP_1 = embed[, 1], UMAP_2 = embed[, 2], broad_lineage = obj$broad_lineage)

centroids <- df %>%
  group_by(broad_lineage) %>%
  summarize(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), .groups = "drop")

lineage_colors <- c(
  "CD4 T" = "#4C78A8", "CD8 T" = "#72B7B2", "NK" = "#E45756",
  "B cell" = "#F58518", "Monocyte" = "#EECA3B", "DC" = "#BAB0AC",
  "MAIT" = "#B279A2", "HSPC" = "#000000", "Proliferating" = "#17BECF",
  "gdT/NK" = "#54A24B", "ILC" = "#9C755F", "pDC" = "#8C564B",
  "Unidentified" = "#7F7F7F"
)

missing_colors <- setdiff(unique(df$broad_lineage), names(lineage_colors))
if (length(missing_colors) > 0) {
  stop(sprintf("No color defined for broad_lineage categories: %s",
               paste(missing_colors, collapse = ", ")))
}

p <- ggplot(df, aes(UMAP_1, UMAP_2, color = broad_lineage)) +
  geom_point(size = 0.15, alpha = 0.35) +
  geom_text_repel(data = centroids, aes(x = UMAP_1, y = UMAP_2, label = broad_lineage),
                   inherit.aes = FALSE, color = "black", size = 5, fontface = "bold",
                   max.overlaps = 100, bg.color = "white", bg.r = 0.15, seed = 1) +
  scale_color_manual(values = lineage_colors, name = "Broad Lineage") +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  labs(title = "RNA UMAP: Broad Level Cell Annotation", x = "UMAP 1", y = "UMAP 2") +
  theme_minimal(base_size = 13) +
  theme(panel.border = element_rect(color = "grey40", fill = NA), panel.grid = element_blank())

ggsave(file.path(figures_dir, "03_rna_umap_broad_labels.png"), p, width = 11, height = 8.5, dpi = 150)
log_msg("Saved results/03_rna_umap_broad_labels.png")

## =============================================================================
## Final checkpoint -- replaces what was previously 03_rna_seurat_object_annotated.rds
## =============================================================================

saveRDS(obj, file.path(results_dir, "03_rna_seurat_object_annotated.rds"))
log_msg("Saved results/03_rna_seurat_object_annotated.rds")
log_msg("Done: %d cells, %d broad lineage categories.", ncol(obj), length(unique(obj$broad_lineage)))
