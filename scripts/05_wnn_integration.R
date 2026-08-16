#!/usr/bin/env Rscript
#
# WNN (weighted nearest neighbor) integration -> UMAP -> clustering -> annotation
#
#   1. FindMultiModalNeighbors, reductions pca+apca.
#   2. UMAP on the WNN graph.
#   3. FindClusters.
#   4. FindAllMarkers.
#   5. Detailed-lineage annotation.
#   6. Detailed-level labeled UMAP plot.
#   7. Broad-lineage annotation + labeled UMAP plot.
#
# Usage:
#   conda activate citeseq-pipeline
#   Rscript 05_wnn_integration.R

suppressMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
})

set.seed(42)

script_args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", script_args[grep("^--file=", script_args)])
project_dir <- dirname(dirname(normalizePath(script_path)))
results_dir <- file.path(project_dir, "results")
figures_dir <- results_dir
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

wnn_dims_rna <- 1:30
wnn_dims_adt <- 1:30
cluster_resolution <- 1.2

log_msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), sprintf(...), "\n", sep = "")
}

log_msg("Loading RNA branch (has pca) and ADT branch (has apca)...")
obj_rna <- readRDS(file.path(results_dir, "03_rna_seurat_object_annotated.rds"))
obj_adt <- readRDS(file.path(results_dir, "04_adt_seurat_object_annotated.rds"))
log_msg("RNA branch: %d cells, reductions: %s", ncol(obj_rna), paste(Reductions(obj_rna), collapse = ", "))
log_msg("ADT branch: %d cells, reductions: %s", ncol(obj_adt), paste(Reductions(obj_adt), collapse = ", "))


stopifnot(setequal(colnames(obj_rna), colnames(obj_adt)))
obj_adt <- obj_adt[, colnames(obj_rna)]


obj <- obj_rna
obj[["ADT"]]      <- obj_adt[["ADT"]]
obj[["apca"]]     <- obj_adt[["apca"]]
obj$nCount_ADT     <- obj_adt$nCount_ADT
obj$nFeature_ADT   <- obj_adt$nFeature_ADT
obj$adt_clusters   <- obj_adt$adt_clusters
obj$adt_broad_lineage        <- obj_adt$adt_broad_lineage
obj$adt_annotation_confidence <- obj_adt$adt_annotation_confidence
rm(obj_rna, obj_adt); gc()

log_msg("Merged: %d cells, reductions: %s", ncol(obj), paste(Reductions(obj), collapse = ", "))

## =============================================================================
## Step 1/7: FindMultiModalNeighbors (WNN)
## =============================================================================

log_msg("Step 1/7: Running FindMultiModalNeighbors (WNN)...")
obj <- FindMultiModalNeighbors(
  obj,
  reduction.list = list("pca", "apca"),
  dims.list = list(wnn_dims_rna, wnn_dims_adt),
  verbose = FALSE
)
log_msg("WNN graphs present: %s", paste(Graphs(obj), collapse = ", "))
log_msg("Modality weight summary (RNA):")
print(summary(obj$RNA.weight))

## =============================================================================
## Step 2/7: UMAP on the WNN graph
## =============================================================================

log_msg("Step 2/7: Running UMAP (WNN graph, seed=42)...")
obj <- RunUMAP(obj, nn.name = "weighted.nn", reduction.name = "wnn.umap",
               reduction.key = "wnnUMAP_", seed.use = 42, verbose = FALSE)

## =============================================================================
## Step 3/7: FindClusters on wsnn
## =============================================================================

log_msg("Step 3/7: FindClusters (wsnn, algorithm=3/SLM, resolution=%.2f, seed=42)...",
        cluster_resolution)
obj <- FindClusters(obj, graph.name = "wsnn", algorithm = 3,
                     resolution = cluster_resolution, random.seed = 42,
                     cluster.name = "wnn_clusters", verbose = FALSE)
n_clusters <- length(unique(obj$wnn_clusters))
log_msg("Found %d clusters", n_clusters)
print(table(obj$wnn_clusters))

## =============================================================================
## Step 4/7: FindAllMarkers, RNA + ADT
## =============================================================================

log_msg("Step 4/7: Running FindAllMarkers (RNA + ADT)...")
Idents(obj) <- obj$wnn_clusters

rna_markers <- FindAllMarkers(obj, assay = "RNA", features = VariableFeatures(obj, assay = "RNA"),
                               only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,
                               max.cells.per.ident = 500, verbose = FALSE)
log_msg("RNA: %d marker rows across %d clusters", nrow(rna_markers), length(unique(rna_markers$cluster)))
write.csv(rna_markers, file.path(results_dir, "05_wnn_cluster_markers_RNA.csv"), row.names = FALSE)
log_msg("Saved results/05_wnn_cluster_markers_RNA.csv")

adt_markers <- FindAllMarkers(obj, assay = "ADT", only.pos = TRUE, min.pct = 0.25,
                               logfc.threshold = 0.25, max.cells.per.ident = 500,
                               verbose = FALSE)
log_msg("ADT: %d marker rows across %d clusters", nrow(adt_markers), length(unique(adt_markers$cluster)))
write.csv(adt_markers, file.path(results_dir, "05_wnn_cluster_markers_ADT.csv"), row.names = FALSE)
log_msg("Saved results/05_wnn_cluster_markers_ADT.csv")

## =============================================================================
## Step 5/7: detailed-lineage annotation (Hao et al. 2021)
## =============================================================================

log_msg("Step 5/7: Detailed-lineage annotation (Hao et al. 2021, celltype.l2)...")

celltype_ref <- read.csv(file.path(project_dir, "reference", "celltype_annotations.csv"),
                          stringsAsFactors = FALSE)
match_idx <- match(colnames(obj), celltype_ref$barcode)
stopifnot(!anyNA(match_idx))  # every cell in our data must exist in the reference

obj$ref_celltype_l2 <- celltype_ref$celltype.l2[match_idx]

majority_label <- function(x) {
  tab <- sort(table(x), decreasing = TRUE)
  names(tab)[1]
}
pct_agreement <- function(x) {
  tab <- sort(table(x), decreasing = TRUE)
  unname(round(100 * tab[1] / sum(tab), 1))
}

broad_lineage_map <- sapply(split(obj$ref_celltype_l2, obj$wnn_clusters), majority_label)
broad_pct_map <- sapply(split(obj$ref_celltype_l2, obj$wnn_clusters), pct_agreement)

obj$wnn_broad_lineage <- unname(broad_lineage_map[as.character(obj$wnn_clusters)])

log_msg("WNN detailed lineage breakdown (from paper reference):")
print(sort(table(obj$wnn_broad_lineage), decreasing = TRUE))

annotation_table <- data.frame(
  cluster = names(broad_lineage_map),
  n_cells = as.integer(table(obj$wnn_clusters)[names(broad_lineage_map)]),
  broad_lineage = unname(broad_lineage_map),
  pct_agreement = unname(broad_pct_map[names(broad_lineage_map)]),
  row.names = NULL
)
write.csv(annotation_table, file.path(results_dir, "05_wnn_detailed_annotation.csv"), row.names = FALSE)
log_msg("Saved results/05_wnn_detailed_annotation.csv")

## =============================================================================
## Step 6/7: detailed-level labeled UMAP plot
## =============================================================================

log_msg("Step 6/7: detailed-level labeled UMAP plot...")

broad_categories <- sort(unique(obj$wnn_broad_lineage))
lineage_colors <- setNames(scales::hue_pal()(length(broad_categories)), broad_categories)

embed <- Embeddings(obj, "wnn.umap")
df_broad <- data.frame(UMAP_1 = embed[, 1], UMAP_2 = embed[, 2], broad_lineage = obj$wnn_broad_lineage)

centroids_broad <- df_broad %>%
  group_by(broad_lineage) %>%
  summarize(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), .groups = "drop")

p_broad <- ggplot(df_broad, aes(UMAP_1, UMAP_2, color = broad_lineage)) +
  geom_point(size = 0.15, alpha = 0.35) +
  geom_text_repel(data = centroids_broad, aes(x = UMAP_1, y = UMAP_2, label = broad_lineage),
                   inherit.aes = FALSE, color = "black", size = 5, fontface = "bold",
                   max.overlaps = 100, bg.color = "white", bg.r = 0.15, seed = 1) +
  scale_color_manual(values = lineage_colors, name = "Detailed Lineage") +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  labs(title = "WNN: Detailed Level Cell Annotation", x = "UMAP 1", y = "UMAP 2") +
  theme_minimal(base_size = 13) +
  theme(panel.border = element_rect(color = "grey40", fill = NA),
        panel.grid = element_blank())

ggsave(file.path(figures_dir, "05_wnn_umap_detailed_labels.png"), p_broad, width = 11, height = 8.5, dpi = 150)
log_msg("Saved results/05_wnn_umap_detailed_labels.png")

## =============================================================================
## Step 7/7: broad-lineage annotation (Hao et al. 2021)
## =============================================================================

log_msg("Step 7/7: Broad-lineage annotation (Hao et al. 2021, celltype.l1)...")

obj$ref_celltype_l1 <- celltype_ref$celltype.l1[match_idx]

l1_map <- sapply(split(obj$ref_celltype_l1, obj$wnn_clusters), majority_label)
l1_pct_map <- sapply(split(obj$ref_celltype_l1, obj$wnn_clusters), pct_agreement)

obj$wnn_l1_lineage <- unname(l1_map[as.character(obj$wnn_clusters)])

log_msg("WNN celltype.l1 breakdown:")
print(sort(table(obj$wnn_l1_lineage), decreasing = TRUE))

l1_table <- data.frame(
  cluster = names(l1_map),
  n_cells = as.integer(table(obj$wnn_clusters)[names(l1_map)]),
  l1_lineage = unname(l1_map),
  pct_agreement = unname(l1_pct_map[names(l1_map)]),
  row.names = NULL
)
write.csv(l1_table, file.path(results_dir, "05_wnn_broad_annotation.csv"), row.names = FALSE)
log_msg("Saved results/05_wnn_broad_annotation.csv")


l1_display_map <- c("other T" = "dnT/gdT/MAIT", "other" = "Other")
l1_display <- ifelse(obj$wnn_l1_lineage %in% names(l1_display_map),
                      unname(l1_display_map[obj$wnn_l1_lineage]),
                      obj$wnn_l1_lineage)

l1_categories <- sort(unique(l1_display))
l1_colors <- setNames(scales::hue_pal()(length(l1_categories)), l1_categories)

df_l1 <- data.frame(UMAP_1 = embed[, 1], UMAP_2 = embed[, 2], l1_lineage = l1_display)

centroids_l1 <- df_l1 %>%
  group_by(l1_lineage) %>%
  summarize(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), .groups = "drop")

p_l1 <- ggplot(df_l1, aes(UMAP_1, UMAP_2, color = l1_lineage)) +
  geom_point(size = 0.15, alpha = 0.35) +
  geom_text_repel(data = centroids_l1, aes(x = UMAP_1, y = UMAP_2, label = l1_lineage),
                   inherit.aes = FALSE, color = "black", size = 5, fontface = "bold",
                   max.overlaps = 100, bg.color = "white", bg.r = 0.15, seed = 1) +
  scale_color_manual(values = l1_colors, name = "Broad Lineage") +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  labs(title = "WNN: Broad Level Cell Annotation", x = "UMAP 1", y = "UMAP 2") +
  theme_minimal(base_size = 13) +
  theme(panel.border = element_rect(color = "grey40", fill = NA),
        panel.grid = element_blank())

ggsave(file.path(figures_dir, "05_wnn_umap_broad_labels.png"), p_l1, width = 11, height = 8.5, dpi = 150)
log_msg("Saved results/05_wnn_umap_broad_labels.png")


saveRDS(obj, file.path(results_dir, "05_wnn_seurat_object_detailed.rds"))
log_msg("Saved results/05_wnn_seurat_object_detailed.rds")
log_msg("Done: %d cells, %d WNN clusters, %d detailed (l2) categories, %d broad (l1) categories.",
        ncol(obj), n_clusters, length(unique(obj$wnn_broad_lineage)), length(unique(obj$wnn_l1_lineage)))
