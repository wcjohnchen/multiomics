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
#   5. Manual broad-lineage annotation using combined RNA+ADT marker
#      evidence, 13 categories (29)
#   6. Broad-level labeled UMAP plot (30)
#   7. Manual detailed-level (per-cluster) annotation, 49 descriptive
#      labels (31)
#   8. Detailed-level UMAP plot: cluster numbers on the plot, colored by
#      broad lineage, grouped legend below (32)
#
# Only the FINAL object is saved as a checkpoint here (see README's
# GitHub section on results/*.rds size). Each step's own CSV/figure
# output is still written individually.
#
# Usage:
#   conda activate citeseq-pipeline
#   Rscript 05_wnn_integration.R

suppressMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
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

log_msg("Loading ADT-annotated object (has both pca and apca)...")
obj <- readRDS(file.path(results_dir, "04_adt_seurat_object_annotated.rds"))
log_msg("%d cells, reductions: %s", ncol(obj), paste(Reductions(obj), collapse = ", "))

## =============================================================================
## Step 1/8 (script 25): FindMultiModalNeighbors (WNN)
## =============================================================================

log_msg("Step 1/8: Running FindMultiModalNeighbors (WNN)...")
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
## Step 2/8 (script 26): UMAP on the WNN graph
## =============================================================================

log_msg("Step 2/8: Running UMAP (WNN graph, seed=42)...")
obj <- RunUMAP(obj, nn.name = "weighted.nn", reduction.name = "wnn.umap",
               reduction.key = "wnnUMAP_", seed.use = 42, verbose = FALSE)

## =============================================================================
## Step 3/8 (script 27): FindClusters on wsnn
## =============================================================================

log_msg("Step 3/8: FindClusters (wsnn, algorithm=3/SLM, resolution=%.2f, seed=42)...",
        cluster_resolution)
obj <- FindClusters(obj, graph.name = "wsnn", algorithm = 3,
                     resolution = cluster_resolution, random.seed = 42,
                     cluster.name = "wnn_clusters", verbose = FALSE)
n_clusters <- length(unique(obj$wnn_clusters))
log_msg("Found %d clusters", n_clusters)
print(table(obj$wnn_clusters))

## =============================================================================
## Step 4/8 (script 28): FindAllMarkers, RNA (2000 HVGs) + ADT (228 antibodies)
## =============================================================================

log_msg("Step 4/8: Running FindAllMarkers (RNA + ADT)...")
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
## Step 5/8 (script 29): manual broad-lineage annotation
## =============================================================================

log_msg("Step 5/8: Manual broad-lineage annotation...")

broad_lineage_map <- c(
  "0"="Monocyte","1"="CD4 T","2"="CD4 T","3"="CD4 T","4"="CD4 T",
  "5"="NK","6"="Monocyte","7"="CD8 T","8"="B cell","9"="CD8 T",
  "10"="Monocyte","11"="NK","12"="Monocyte","13"="B cell","14"="Monocyte",
  "15"="CD8 T","16"="Monocyte","17"="CD8 T","18"="CD8 T","19"="Monocyte",
  "20"="CD8 T","21"="NK","22"="Monocyte","23"="B cell","24"="Monocyte",
  "25"="CD8 T","26"="CD4 T","27"="Monocyte","28"="MAIT","29"="Monocyte",
  "30"="B cell","31"="CD8 T","32"="DC","33"="NK","34"="MAIT",
  "35"="NK","36"="CD4 T","37"="gdT","38"="Platelet","39"="CD4 T",
  "40"="gdT","41"="Monocyte","42"="Proliferating","43"="CD8 T","44"="NK",
  "45"="HSPC","46"="ILC","47"="CD4 T","48"="Erythrocyte"
)

confidence_map <- c(
  "0"="high","1"="high","2"="high","3"="high","4"="high",
  "5"="high","6"="high","7"="high","8"="high","9"="high",
  "10"="high","11"="high","12"="high","13"="high","14"="high",
  "15"="high","16"="medium","17"="high","18"="medium","19"="high",
  "20"="high","21"="high","22"="high","23"="high","24"="high",
  "25"="high","26"="high","27"="high","28"="high","29"="high",
  "30"="high","31"="high","32"="high","33"="high","34"="high",
  "35"="high","36"="low","37"="high","38"="high","39"="medium",
  "40"="high","41"="high","42"="high","43"="high","44"="low",
  "45"="high","46"="low","47"="high","48"="high"
)

key_markers_map <- c(
  "0"="RNA: VSTM1,RETN / ADT: CD155,CD86,CD14",
  "1"="RNA: CD40LG,IL7R / ADT: CD4,CD25,CD127,CD278",
  "2"="RNA: CCR7,LEF1 / ADT: CD4,CD45RB,CD27",
  "3"="RNA: CCR6,TNFRSF4 / ADT: CD4,CD25,CD45RO",
  "4"="RNA: LEF1,CCR7 / ADT: CD4,CD3,CD27",
  "5"="RNA: KLRC1,NCR1 / ADT: CD335,CD56,CD16",
  "6"="RNA: SIGLEC1,MARCO / ADT: CD169,CD64,CD11b",
  "7"="RNA: CD8A/B,LEF1,CCR7 / ADT: CD8,CD8a,CD27",
  "8"="RNA: TCL1A,IGHD / ADT: IgD,CD21,CD20,IgM",
  "9"="RNA: CD8A/B,CCR7 / ADT: CD8,CD45RA,CD27",
  "10"="RNA: IL1A,CCL20,IL6 / ADT: CD36,CD64,CD11b",
  "11"="RNA: KLRF1,NCR1,CD160 / ADT: CD335,CD158,CD56",
  "12"="RNA: TNF,NFKBIA,CCL3 / ADT: CD36,CD64,CD86,CD11b",
  "13"="RNA: TNFRSF13B,IGHG / ADT: CD20,CD21,CD19",
  "14"="RNA: HLA-DQA2,IFI30 / ADT: CD155,CD64,CD192",
  "15"="RNA: GZMH,EOMES,CD8A/B / ADT: CD8,TIGIT,CD8a",
  "16"="RNA: CDKN1C,C1QA / ADT: CD86,CD123,CD16",
  "17"="RNA: GZMK,CD8A/B / ADT: CD8,CD314,CD8a",
  "18"="RNA: KLRC2,ZNF683,GZMH / ADT: CD8,CD158,CD56",
  "19"="RNA: THBS1,PTGES,TNFAIP6 / ADT: CD36,CD64,CD11b",
  "20"="RNA: GZMK,CD8A/B / ADT: CD8,CD103,CD314",
  "21"="RNA: KLRC2,KIR2DL3,TRDC / ADT: CD158b,CD335,CD56",
  "22"="RNA: IFI27,EGR2/3,SIGLEC1 / ADT: CD169,CD64,CD86",
  "23"="RNA: TCL1A,IGHD,IGHM / ADT: CD20,IgD,CD19",
  "24"="RNA: VMO1,CDKN1C,C1QA / ADT: CD123,CD16,CD86",
  "25"="RNA: ZNF683,GZMH,GNLY / ADT: CD8,CD56,CD314",
  "26"="RNA: LEF1,CCR7 / ADT: CD4,CD278,CD27",
  "27"="RNA: TMEM144,CYP1B1 / ADT: CD169,CD36,CD64",
  "28"="RNA: SLC4A10,GZMK,KLRB1 / ADT: TCR-Va7.2(MAIT),CD161",
  "29"="RNA: RETN,S100A8/12 / ADT: CD155,CD1d,CD86",
  "30"="RNA: IGHA2,IGHG,JCHAIN / ADT: CD20,IgM,CD19 (plasma-like)",
  "31"="RNA: CD8A/B,NELL2 / ADT: CD8,CD314,CD45RA",
  "32"="RNA: FCER1A,CD1C,CLEC10A,FLT3 / ADT: CD1c,HLA-DR (cDC2)",
  "33"="RNA: XCL1,XCL2,NCAM1,KLRC1 / ADT: CD335,CD56,CD117 (CD56bright)",
  "34"="RNA: SLC4A10,GZMK,CCR6 / ADT: TCR-Va7.2(MAIT),CD161",
  "35"="RNA: KIR2DL1,NCAM1 / ADT: CD158b,CD158e1,CD56,CD57",
  "36"="RNA: ZNF683,GZMH (mixed) / ADT: CD271,CD4,CD307e",
  "37"="RNA: TRDV2,TRGV9,TRDC (gdT genes) / ADT: TCR-Vd2,TCR-Vg9",
  "38"="RNA: GP9,ITGA2B,PF4,PPBP / ADT: CD42b,CD61,CD9 (platelet)",
  "39"="RNA: GZMH,MAF,CD40LG / ADT: CD4,CD195,CD3",
  "40"="RNA: TRGC1,TRDC,KLRC1 (gdT genes) / ADT: TCR-Vd2,TCR-Vg9",
  "41"="RNA: C1QA/B,CDKN1C / ADT: CD86,CD172a,CD204",
  "42"="RNA: MKI67,TOP2A,TYMS,BIRC5 (proliferating) / ADT: CD335,CD56",
  "43"="RNA: CD8A/B,ITM2C,CCL5 / ADT: CD103,CD8,Integrin-7",
  "44"="RNA: SPON2,KLRF1,FGFBP2 / ADT: CD319,VEGFR-3 (atypical)",
  "45"="RNA: AVP,PROM1,CD34 / ADT: CD117,CD133,CD34 (HSPC)",
  "46"="RNA: KIT,IL1R1,IL2RA,SPINK2 / ADT: CD117,CD25,CD161",
  "47"="RNA: GATA3,TRAC,TRBC1,LTB / ADT: CD4,CD30,CD27",
  "48"="RNA: HBB,HBA1/2,ALAS2 / ADT: CD235a/ab (erythrocyte)"
)

clusters_present <- as.character(sort(unique(as.integer(as.character(obj$wnn_clusters)))))
stopifnot(all(clusters_present %in% names(broad_lineage_map)))

obj$wnn_broad_lineage <- unname(broad_lineage_map[as.character(obj$wnn_clusters)])
obj$wnn_annotation_confidence <- unname(confidence_map[as.character(obj$wnn_clusters)])

log_msg("WNN broad lineage breakdown:")
print(sort(table(obj$wnn_broad_lineage), decreasing = TRUE))

annotation_table <- data.frame(
  cluster = names(broad_lineage_map),
  n_cells = as.integer(table(obj$wnn_clusters)[names(broad_lineage_map)]),
  broad_lineage = broad_lineage_map,
  confidence = confidence_map,
  key_markers = key_markers_map,
  row.names = NULL
)
write.csv(annotation_table, file.path(results_dir, "05_wnn_broad_annotation.csv"), row.names = FALSE)
log_msg("Saved results/05_wnn_broad_annotation.csv")

## =============================================================================
## Step 6/8 (script 30): broad-level labeled UMAP plot
## =============================================================================

log_msg("Step 6/8: Building broad-level labeled UMAP plot...")

lineage_colors <- c(
  "CD4 T" = "#4C78A8", "CD8 T" = "#72B7B2", "gdT" = "#54A24B",
  "MAIT" = "#B279A2", "NK" = "#E45756", "B cell" = "#F58518",
  "Monocyte" = "#EECA3B", "DC" = "#BAB0AC", "HSPC" = "#000000",
  "Platelet" = "#FF9DA6", "Erythrocyte" = "#8C564B", "Proliferating" = "#17BECF",
  "ILC" = "#9C755F"
)

embed <- Embeddings(obj, "wnn.umap")
df_broad <- data.frame(UMAP_1 = embed[, 1], UMAP_2 = embed[, 2], broad_lineage = obj$wnn_broad_lineage)

centroids_broad <- df_broad %>%
  group_by(broad_lineage) %>%
  summarize(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), .groups = "drop")

missing_colors <- setdiff(unique(df_broad$broad_lineage), names(lineage_colors))
if (length(missing_colors) > 0) {
  stop(sprintf("No color defined for broad_lineage categories: %s",
               paste(missing_colors, collapse = ", ")))
}

p_broad <- ggplot(df_broad, aes(UMAP_1, UMAP_2, color = broad_lineage)) +
  geom_point(size = 0.15, alpha = 0.35) +
  geom_text_repel(data = centroids_broad, aes(x = UMAP_1, y = UMAP_2, label = broad_lineage),
                   inherit.aes = FALSE, color = "black", size = 5, fontface = "bold",
                   max.overlaps = 100, bg.color = "white", bg.r = 0.15, seed = 1) +
  scale_color_manual(values = lineage_colors, name = "Broad Lineage") +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  labs(title = "WNN: Broad Level Cell Annotation", x = "UMAP 1", y = "UMAP 2") +
  theme_minimal(base_size = 13) +
  theme(panel.border = element_rect(color = "grey40", fill = NA),
        panel.grid = element_blank())

ggsave(file.path(figures_dir, "05_wnn_umap_broad_labels.png"), p_broad, width = 11, height = 8.5, dpi = 150)
log_msg("Saved results/05_wnn_umap_broad_labels.png")

## =============================================================================
## Step 7/8 (script 31): manual detailed-level (per-cluster) annotation
## =============================================================================

log_msg("Step 7/8: Manual detailed-level annotation...")

fine_label_map <- c(
  "0"  = "Monocyte - Classical",
  "1"  = "CD4 T - CXCR5+/ICOS+ (Tfh-like)",
  "2"  = "CD4 T - Naive (CCR7+)",
  "3"  = "CD4 T - CCR6+ (Th17-like)",
  "4"  = "CD4 T - Naive subset 2",
  "5"  = "NK - CD16+ CD56dim",
  "6"  = "Monocyte - Non-classical/IFN-stimulated",
  "7"  = "CD8 T - Naive",
  "8"  = "B cell - Naive IgD+",
  "9"  = "CD8 T - Naive subset 2",
  "10" = "Monocyte - Inflammatory (IL1/IL6/CCL20)",
  "11" = "NK - CD16+ (KLRF1+)",
  "12" = "Monocyte - Inflammatory (TNF+)",
  "13" = "B cell - Memory/class-switched",
  "14" = "Monocyte - HLA-DQA2+",
  "15" = "CD8 T - Cytotoxic/Effector (GZMH+EOMES+)",
  "16" = "Monocyte - Non-classical CD16+",
  "17" = "CD8 T - Effector Memory (GZMK+)",
  "18" = "CD8 T - Adaptive/tissue-resident (KIR+)",
  "19" = "Monocyte - Inflammatory subset 2",
  "20" = "CD8 T - Tissue-resident (CD103+)",
  "21" = "NK - Adaptive (KIR+)",
  "22" = "Monocyte - IFN-stimulated",
  "23" = "B cell - Naive subset 2",
  "24" = "Monocyte - Non-classical subset 2",
  "25" = "CD8 T - Terminal effector (GNLY+)",
  "26" = "CD4 T - Naive subset 3",
  "27" = "Monocyte - CYP1B1+",
  "28" = "MAIT - subset 1",
  "29" = "Monocyte - Classical subset 2",
  "30" = "B cell - Plasmablast-like (JCHAIN+)",
  "31" = "CD8 T - Naive subset 3",
  "32" = "DC - cDC2",
  "33" = "NK - CD56bright",
  "34" = "MAIT - subset 2",
  "35" = "NK - Adaptive (KIR+ CD57+)",
  "36" = "CD4 T - CTL-like *",
  "37" = "gdT - Vgamma9/Vdelta2",
  "38" = "Platelet",
  "39" = "CD4 T - Cytotoxic (GZMH+)",
  "40" = "gdT - subset 2",
  "41" = "Monocyte - Macrophage-like (C1Q+)",
  "42" = "Proliferating - NK",
  "43" = "CD8 T - Tissue-resident subset 2",
  "44" = "NK - atypical *",
  "45" = "HSPC - CD34+CD117+",
  "46" = "ILC *",
  "47" = "CD4 T - GATA3+ (Th2-like)",
  "48" = "Erythrocyte (ambient RBC contamination)"
)

stopifnot(all(clusters_present %in% names(fine_label_map)))
obj$wnn_fine_label <- unname(fine_label_map[as.character(obj$wnn_clusters)])

fine_table <- data.frame(
  cluster = names(fine_label_map),
  n_cells = as.integer(table(obj$wnn_clusters)[names(fine_label_map)]),
  broad_lineage = unname(setNames(obj$wnn_broad_lineage, obj$wnn_clusters)[names(fine_label_map)]),
  fine_label = fine_label_map,
  row.names = NULL
)
write.csv(fine_table, file.path(results_dir, "05_wnn_detailed_annotation.csv"), row.names = FALSE)
log_msg("Saved results/05_wnn_detailed_annotation.csv")

## =============================================================================
## Step 8/8 (script 32): detailed-level labeled UMAP plot
## =============================================================================

log_msg("Step 8/8: Building detailed-level labeled UMAP plot...")

df_fine <- data.frame(UMAP_1 = embed[, 1], UMAP_2 = embed[, 2],
                       cluster = as.character(obj$wnn_clusters),
                       broad_lineage = obj$wnn_broad_lineage)

centroids_fine <- df_fine %>%
  group_by(cluster) %>%
  summarize(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2),
            broad_lineage = dplyr::first(broad_lineage), .groups = "drop")

p_main <- ggplot(df_fine, aes(UMAP_1, UMAP_2, color = broad_lineage)) +
  geom_point(size = 0.15, alpha = 0.35) +
  geom_text_repel(data = centroids_fine, aes(x = UMAP_1, y = UMAP_2, label = cluster), inherit.aes = FALSE,
                   color = "black", size = 3, fontface = "bold", max.overlaps = 100,
                   bg.color = "white", bg.r = 0.15, seed = 1) +
  scale_color_manual(values = lineage_colors, name = "Broad Lineage") +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  labs(title = "WNN: Detailed Level Cell Annotation", x = "UMAP 1", y = "UMAP 2") +
  theme_minimal(base_size = 12) +
  theme(panel.border = element_rect(color = "grey40", fill = NA), panel.grid = element_blank(),
        plot.title = element_text(size = 20),
        legend.title = element_text(size = 16))

legend_df <- fine_table %>%
  mutate(cluster = as.character(cluster)) %>%
  mutate(entry = sprintf("%s: %s", cluster, fine_label)) %>%
  arrange(broad_lineage, as.integer(cluster))

legend_panel <- function(group_name) {
  d <- legend_df %>% filter(broad_lineage == group_name)
  d$y <- rev(seq_len(nrow(d)))
  ggplot(d, aes(x = 0, y = y)) +
    geom_text(aes(label = entry), hjust = 0, size = 3) +
    scale_x_continuous(limits = c(-0.1, 6)) +
    scale_y_continuous(limits = c(0.3, nrow(d) + 0.7)) +
    labs(title = group_name) +
    theme_void(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 12, hjust = 0, margin = margin(b = 4)))
}

group_order <- c("CD4 T", "CD8 T", "NK", "Monocyte", "B cell", "DC",
                 "gdT", "MAIT", "dnT", "HSPC", "Platelet", "Erythrocyte",
                 "Proliferating", "ILC")
group_order <- group_order[group_order %in% unique(legend_df$broad_lineage)]
legend_panels <- lapply(group_order, legend_panel)
legend_grid <- wrap_plots(legend_panels, ncol = 4)

p_combined <- p_main / legend_grid + plot_layout(heights = c(2, 1.6))

ggsave(file.path(figures_dir, "05_wnn_umap_detailed_labels.png"), p_combined, width = 17, height = 16, dpi = 150)
log_msg("Saved results/05_wnn_umap_detailed_labels.png")

## =============================================================================
## Final checkpoint -- replaces what was previously 05_wnn_seurat_object_detailed.rds
## =============================================================================

saveRDS(obj, file.path(results_dir, "05_wnn_seurat_object_detailed.rds"))
log_msg("Saved results/05_wnn_seurat_object_detailed.rds")
log_msg("Done: %d cells, %d WNN clusters, %d broad lineage categories.",
        ncol(obj), n_clusters, length(unique(obj$wnn_broad_lineage)))
