#!/usr/bin/env Rscript
#
# Combined ADT normalization -> UMAP -> clustering -> annotation pipeline
# -- consolidates scripts 18, 19, 20, 21, 22, 23, 24. Covers everything
# from the RNA-annotated object (which already carries the shared 228-
# antibody ADT counts) to the final labeled ADT-only UMAP plot.
#
#   1. CLR normalization, margin=2 (script 18)
#   2. Scale + PCA, all 228 antibodies, 30 PCs (script 19) -- no HVG-style
#      subselection; unlike RNA's 17,808 genes, all 228 antibodies were
#      deliberately curated to be informative, so all are used directly.
#   3. UMAP on apca (script 20)
#   4. FindNeighbors + FindClusters, 29 clusters, using distinct
#      `adt_clusters`/`adt_snn` names so as not to collide with the RNA
#      clustering already on this object (script 21)
#   5. FindAllMarkers, all 228 antibodies (script 22)
#   6. Manual broad-lineage annotation, using ADT surface-protein markers
#      -- several near-definitive single-protein calls (TCR-Va7.2 for
#      MAIT, TCR-Vg9/Vd2 for gdT, CD34/CD117/CD133 for HSPC) (script 23)
#   7. Labeled ADT UMAP plot (script 24)
#
# Only the FINAL object is saved as a checkpoint here (see README's
# GitHub section on results/*.rds size). Each step's own CSV/figure
# output is still written individually.
#
# Usage:
#   conda activate citeseq-pipeline
#   Rscript 04_adt_umap.R

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

adt_npcs <- 30
cluster_resolution <- 0.5

log_msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), sprintf(...), "\n", sep = "")
}

log_msg("Loading RNA-annotated Seurat object...")
obj <- readRDS(file.path(results_dir, "03_rna_seurat_object_annotated.rds"))
log_msg("%d cells, %d antibodies (ADT)", ncol(obj), nrow(obj[["ADT"]]))

## =============================================================================
## Step 1/7 (script 18): CLR normalization
## =============================================================================

log_msg("Step 1/7: Normalizing ADT (CLR, margin=2)...")
obj <- NormalizeData(obj, assay = "ADT", normalization.method = "CLR",
                      margin = 2, verbose = FALSE)

raw_adt  <- GetAssayData(obj, assay = "ADT", layer = "counts")[, 1]
norm_adt <- GetAssayData(obj, assay = "ADT", layer = "data")[, 1]
example_ab <- names(raw_adt)[raw_adt > 0][1]
log_msg("Sanity check (cell 1, antibody %s): raw count=%d, CLR value=%.4f",
        example_ab, raw_adt[example_ab], norm_adt[example_ab])
rm(raw_adt, norm_adt); gc()

## =============================================================================
## Step 2/7 (script 19): Scale + PCA, all 228 antibodies
## =============================================================================

log_msg("Step 2/7: Scaling data (ADT, all 228 antibodies) + PCA (seed=42)...")
adt_features <- rownames(obj[["ADT"]])
VariableFeatures(obj, assay = "ADT") <- adt_features

obj <- ScaleData(obj, assay = "ADT", features = adt_features, verbose = FALSE)
obj <- RunPCA(obj, assay = "ADT", features = adt_features, npcs = adt_npcs,
              reduction.name = "apca", seed.use = 42, verbose = FALSE)
obj[["ADT"]]$scale.data <- NULL
gc()

pca_var <- obj[["apca"]]@stdev^2
pct_var <- round(100 * pca_var / sum(pca_var), 2)
log_msg("PC1-5 %% variance explained: %s", paste(pct_var[1:5], collapse = ", "))

## =============================================================================
## Step 3/7 (script 20): UMAP on apca
## =============================================================================

log_msg("Step 3/7: Running UMAP (ADT apca, dims 1:30, seed=42)...")
obj <- RunUMAP(obj, reduction = "apca", dims = 1:30, reduction.name = "umap.adt",
               reduction.key = "adtUMAP_", seed.use = 42, verbose = FALSE)

## =============================================================================
## Step 4/7 (script 21): FindNeighbors + FindClusters
## =============================================================================

log_msg("Step 4/7: FindNeighbors + FindClusters (algorithm=3/SLM, resolution=%.2f, seed=42)...",
        cluster_resolution)
obj <- FindNeighbors(obj, reduction = "apca", dims = 1:30,
                      graph.name = c("adt_nn", "adt_snn"), verbose = FALSE)
obj <- FindClusters(obj, graph.name = "adt_snn", algorithm = 3,
                     resolution = cluster_resolution, random.seed = 42,
                     cluster.name = "adt_clusters", verbose = FALSE)
n_clusters <- length(unique(obj$adt_clusters))
log_msg("Found %d clusters", n_clusters)
print(table(obj$adt_clusters))

## =============================================================================
## Step 5/7 (script 22): FindAllMarkers, all 228 antibodies
## =============================================================================

log_msg("Step 5/7: Running FindAllMarkers (ADT, all 228 antibodies)...")
Idents(obj) <- obj$adt_clusters
markers <- FindAllMarkers(obj, assay = "ADT", only.pos = TRUE, min.pct = 0.25,
                           logfc.threshold = 0.25, max.cells.per.ident = 500,
                           verbose = FALSE)
log_msg("Found %d marker rows across %d clusters", nrow(markers), length(unique(markers$cluster)))
write.csv(markers, file.path(results_dir, "04_adt_cluster_markers.csv"), row.names = FALSE)
log_msg("Saved results/04_adt_cluster_markers.csv")

## =============================================================================
## Step 6/7 (script 23): manual broad-lineage annotation
## =============================================================================

log_msg("Step 6/7: Manual broad-lineage annotation...")

broad_lineage_map <- c(
  "0"="Monocyte","1"="CD4 T","2"="NK","3"="Monocyte","4"="CD4 T",
  "5"="CD8 T","6"="CD4 T","7"="NK","8"="CD8 T","9"="B cell",
  "10"="CD8 T","11"="CD8 T","12"="CD4 T","13"="B cell","14"="B cell",
  "15"="Monocyte","16"="Monocyte","17"="DC","18"="Monocyte","19"="MAIT",
  "20"="CD8 T","21"="MAIT","22"="gdT","23"="CD4 T","24"="CD4 T",
  "25"="gdT","26"="CD8 T","27"="CD8 T","28"="HSPC"
)

confidence_map <- c(
  "0"="high","1"="high","2"="high","3"="medium","4"="high","5"="high",
  "6"="high","7"="high","8"="high","9"="high","10"="high","11"="high",
  "12"="high","13"="high","14"="high","15"="high","16"="medium","17"="high",
  "18"="high","19"="high","20"="high","21"="high","22"="high","23"="low",
  "24"="high","25"="high","26"="high","27"="high","28"="high"
)

key_markers_map <- c(
  "0"="CD36, CD155, CLEC12A, CD64, CD86, CD11b",
  "1"="CD4, CD185(CXCR5), CD278(ICOS), CD3",
  "2"="CD335(NKp46), CD56, CD158b, CD16, CD337(NKp30)",
  "3"="CD155, CLEC12A, CD1d, CD86, Folate receptor",
  "4"="CD4, CD106, CD209, CD45RB, CD27, CD3",
  "5"="CD8, CD8a, TIGIT, CD314(NKG2D), CD3",
  "6"="CD4, CD106, CD209, CD45RO, CD127",
  "7"="CD335, CD158, CD56, CD337, CD16, CD161",
  "8"="CD8, CD8a, CD45RB, CD27, CD45RA (naive)",
  "9"="IgD, CD21, CD20, CD72, CD268, CD22, IgM",
  "10"="CD8, CD8a, CD185, CD275, CD314",
  "11"="CD8, CD8a, CD314, CD275, TIGIT",
  "12"="CD4, CD278(ICOS), CD25, CD127, CD28",
  "13"="CD20, CD21, CD268, CD22, CD19, CD1c",
  "14"="CD20, CD21, IgD, CD72, CD22, IgM, CD19",
  "15"="CD155, CD36, CD86, CD1d, CD169, CD64",
  "16"="CD123, CD86, CD16, CLEC12A, CD13",
  "17"="CD123, CD86, CD16, CD141(cDC1), CD172a",
  "18"="CD169, CD64, CD155, CLEC12A, CD36",
  "19"="TCR-Va7.2 (MAIT), CD161, CD117, CD195",
  "20"="CD8, CD314, CD8a, CD45RB, CD45RA",
  "21"="TCR-Va7.2 (MAIT), CD161, CD26, CD195",
  "22"="TCR-Vd2, TCR-Vg9 (gdT), CD161, CD3",
  "23"="CD271, CD158e1, CD4, CD307e (mixed/unclear)",
  "24"="CD195, CD4, CD275, CD185, CD271",
  "25"="TCR-Vd2, TCR-Vg9 (gdT), CD161, CD3",
  "26"="CD103, CD8, CD96, Integrin-7 (tissue-resident)",
  "27"="CD103, CD8, Integrin-7, CD96, CD73",
  "28"="CD117, CD133, CD34, CD71 (HSPC)"
)

clusters_present <- as.character(sort(unique(as.integer(as.character(obj$adt_clusters)))))
stopifnot(all(clusters_present %in% names(broad_lineage_map)))

obj$adt_broad_lineage <- unname(broad_lineage_map[as.character(obj$adt_clusters)])
obj$adt_annotation_confidence <- unname(confidence_map[as.character(obj$adt_clusters)])

log_msg("ADT broad lineage breakdown:")
print(sort(table(obj$adt_broad_lineage), decreasing = TRUE))

annotation_table <- data.frame(
  cluster = names(broad_lineage_map),
  n_cells = as.integer(table(obj$adt_clusters)[names(broad_lineage_map)]),
  broad_lineage = broad_lineage_map,
  confidence = confidence_map,
  key_markers = key_markers_map,
  row.names = NULL
)
write.csv(annotation_table, file.path(results_dir, "04_adt_broad_annotation.csv"), row.names = FALSE)
log_msg("Saved results/04_adt_broad_annotation.csv")

## =============================================================================
## Step 7/7 (script 24): labeled ADT UMAP plot
## =============================================================================

log_msg("Step 7/7: Building labeled ADT UMAP plot...")

embed <- Embeddings(obj, "umap.adt")
df <- data.frame(UMAP_1 = embed[, 1], UMAP_2 = embed[, 2], broad_lineage = obj$adt_broad_lineage)

centroids <- df %>%
  group_by(broad_lineage) %>%
  summarize(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), .groups = "drop")

lineage_colors <- c(
  "CD4 T" = "#4C78A8", "CD8 T" = "#72B7B2", "gdT" = "#54A24B",
  "MAIT" = "#B279A2", "NK" = "#E45756", "B cell" = "#F58518",
  "Monocyte" = "#EECA3B", "DC" = "#BAB0AC", "HSPC" = "#000000"
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
  labs(title = "ADT UMAP: Broad Level Cell Annotation", x = "UMAP 1", y = "UMAP 2") +
  theme_minimal(base_size = 13) +
  theme(panel.border = element_rect(color = "grey40", fill = NA), panel.grid = element_blank())

ggsave(file.path(figures_dir, "04_adt_umap_broad_labels.png"), p, width = 11, height = 8.5, dpi = 150)
log_msg("Saved results/04_adt_umap_broad_labels.png")

## =============================================================================
## Final checkpoint -- replaces what was previously 04_adt_seurat_object_annotated.rds
## =============================================================================

saveRDS(obj, file.path(results_dir, "04_adt_seurat_object_annotated.rds"))
log_msg("Saved results/04_adt_seurat_object_annotated.rds")
log_msg("Done: %d cells, %d broad lineage categories.", ncol(obj), length(unique(obj$adt_broad_lineage)))
