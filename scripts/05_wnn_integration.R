#!/usr/bin/env Rscript
#
# Combined WNN (weighted nearest neighbor) integration -> UMAP -> clustering
# -> annotation pipeline -- consolidates scripts 25, 26, 27, 28, 29, 30, 31,
# 32. This is the actual "combine RNA and ADT" phase: each modality was
# normalized/scaled/PCA'd separately (03_rna_umap.R, 04_adt_umap.R); WNN
# is where they're integrated.
#
#   1. FindMultiModalNeighbors, reductions pca+apca, dims 1:30 each (25)
#   2. UMAP on the WNN graph (26)
#   3. FindClusters on wsnn, 49 clusters, `wnn_clusters` to avoid
#      colliding with the existing RNA (`seurat_clusters`) and ADT
#      (`adt_clusters`) columns already on this object (27)
#   4. FindAllMarkers, both RNA (2000 HVGs) and ADT (228 antibodies) (28)
#   5. Detailed-lineage annotation from Hao et al. 2021's own celltype.l2
#      reference labels (majority vote per cluster) -- exact barcode
#      match confirmed, reference/celltype_annotations.csv (29)
#   6. Detailed-level labeled UMAP plot (30)
#   7. Broad-lineage (celltype.l1, 8 categories) annotation + labeled
#      UMAP plot -- same reference, same majority-vote method as step 5,
#      shown alongside celltype.l2 for direct comparison of resolution
#
# No finer-grained (celltype.l3) annotation step -- see the note
# above the final checkpoint for why celltype.l2 is the level that
# actually suits this dataset's clustering, not celltype.l3.
#
# Only the FINAL object is saved as a checkpoint here (see README's
# GitHub section on results/*.rds size). Each step's own CSV/figure
# output is still written individually...
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

# Auto-detect the project root from this script's own location (scripts/
# is always one level below it), so paths work regardless of where the
# repo is cloned.
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

# The two branches were filtered independently (RNA QC vs. ADT antibody/
# cell filter, both applied on top of the same doublet-filtered
# checkpoint) -- confirm they still agree on the exact same cell set
# before merging, rather than assuming it.
stopifnot(setequal(colnames(obj_rna), colnames(obj_adt)))
obj_adt <- obj_adt[, colnames(obj_rna)]

# Merge: start from the RNA branch (has pca + RNA annotations), then
# bring in the ADT branch's own filtered/normalized ADT assay, its apca
# reduction, and its ADT-specific annotations -- replacing the stale,
# unfiltered ADT assay obj_rna still carries from before the two
# branches split (03_rna_umap.R never touches ADT).
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
## Step 1/7 (script 25): FindMultiModalNeighbors (WNN)
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
## Step 2/7 (script 26): UMAP on the WNN graph
## =============================================================================

log_msg("Step 2/7: Running UMAP (WNN graph, seed=42)...")
obj <- RunUMAP(obj, nn.name = "weighted.nn", reduction.name = "wnn.umap",
               reduction.key = "wnnUMAP_", seed.use = 42, verbose = FALSE)

## =============================================================================
## Step 3/7 (script 27): FindClusters on wsnn
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
## Step 4/7 (script 28): FindAllMarkers, RNA (2000 HVGs) + ADT (228 antibodies)
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
## Step 5/7 (script 29): detailed-lineage annotation, from Hao et al. 2021's
## own celltype.l2 reference labels (majority vote per WNN cluster) --
## replaces the earlier manual marker-based lookup. See README
## "Methodology notes" for why: this is a real, citable ground truth for
## this exact dataset (confirmed by exact barcode match), not an AI/
## human guess -- unlike the RNA-only and ADT-only annotation steps,
## which stay manual on purpose (see README).
## =============================================================================

log_msg("Step 5/7: Detailed-lineage annotation (Hao et al. 2021 celltype.l2, majority vote per cluster)...")

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
  # unname() is required: tab[1] carries its own category name (e.g.
  # "CD14 Mono"), and if left in, sapply() below would concatenate it
  # into the result's own name (e.g. "0.CD14 Mono" instead of "0"),
  # breaking the later broad_pct_map[names(broad_lineage_map)] lookup
  # -- confirmed by hitting exactly that bug (silent all-NA column).
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
write.csv(annotation_table, file.path(results_dir, "05_wnn_broad_annotation.csv"), row.names = FALSE)
log_msg("Saved results/05_wnn_broad_annotation.csv")

## =============================================================================
## Step 6/7 (script 30): detailed-level labeled UMAP plot
## =============================================================================

log_msg("Step 6/7: Building detailed-level labeled UMAP plot...")

# Auto-generated to exactly match whatever categories the paper's
# reference actually produces (no longer a fixed hand-picked set, since
# celltype.l2 has 31 possible values -- only some appear as a cluster
# majority here).
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

ggsave(file.path(figures_dir, "05_wnn_umap_broad_labels.png"), p_broad, width = 11, height = 8.5, dpi = 150)
log_msg("Saved results/05_wnn_umap_broad_labels.png")

## =============================================================================
## Step 7/7: broad-lineage annotation from celltype.l1 (8 categories) --
## same reference, same majority-vote method as Step 5/7, kept alongside
## celltype.l2 (not instead of it) for direct comparison. celltype.l1 was
## tested and rejected as the *primary* annotation earlier (it merges
## MAIT+gdT into "other T" and HSPC+ILC+Eryth into "other", losing real,
## spatially-distinct populations visible on the UMAP) -- this plot
## exists specifically to make that resolution loss visible, not as an
## equally-valid alternative to Step 5/7.
## =============================================================================

log_msg("Step 7/7: Broad-lineage annotation (Hao et al. 2021 celltype.l1, majority vote per cluster)...")

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
write.csv(l1_table, file.path(results_dir, "05_wnn_l1_annotation.csv"), row.names = FALSE)
log_msg("Saved results/05_wnn_l1_annotation.csv")

# Display-only relabel for the plot -- these are the paper's real
# celltype.l1 values (kept as-is in results/05_wnn_l1_annotation.csv
# above, unchanged): "other T" is a vague catch-all covering three
# distinct T-cell subtypes (gdT/MAIT/dnT, confirmed by cross-referencing
# the reference's own celltype.l2 column); "other" is just capitalized
# for consistency with the other title-case legend entries (CD4 T, CD8
# T, etc.). Relabeled here for the plot only.
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

ggsave(file.path(figures_dir, "05_wnn_umap_l1_labels.png"), p_l1, width = 11, height = 8.5, dpi = 150)
log_msg("Saved results/05_wnn_umap_l1_labels.png")

## =============================================================================
## Final checkpoint -- replaces what was previously 05_wnn_seurat_object_detailed.rds
##
## No separate detailed-level (celltype.l3) annotation step: l3 has 58
## possible categories but only ~23 actually appear as a majority in our
## 49 clusters, and its granularity is uneven relative to our own
## clustering (e.g. only 2 monocyte categories total, so 10 of our 49
## clusters all collapse to the same "CD14 Mono" label) -- broad
## (celltype.l2) is the level that actually suits this clustering.
## =============================================================================

saveRDS(obj, file.path(results_dir, "05_wnn_seurat_object_detailed.rds"))
log_msg("Saved results/05_wnn_seurat_object_detailed.rds")
log_msg("Done: %d cells, %d WNN clusters, %d detailed (l2) categories, %d broad (l1) categories.",
        ncol(obj), n_clusters, length(unique(obj$wnn_broad_lineage)), length(unique(obj$wnn_l1_lineage)))
