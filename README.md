# Integrated Analysis of Single-Cell Multiomics Data Using CITE-seq

A from-scratch, step-by-step rebuild of the CITE-seq analysis pipeline
already validated in `citeseq_learn_v2`, built incrementally in a fresh
conda environment (`citeseq-pipeline`) rather than cloning the working
environment wholesale. Each step is run and checked individually before
moving to the next. The pipeline is now complete end-to-end: raw data →
QC/filtering → independent RNA and ADT processing → WNN integration →
final combined annotation.

📊 **[CITE-seq Analysis Report](https://claude.ai/code/artifact/7966cdf7-365d-4a2b-8d4f-81e619d4259e)**

## Contents

- [1. Data](#1-data)
- [2. Computational Methods](#2-computational-methods)
- [3. Environment Setup](#3-environment-setup)
- [4. Analysis Workflow](#4-analysis-workflow)
- [5. Running the Analysis Workflow](#5-running-the-analysis-workflow)
- [6. Directory Structure](#6-directory-structure)
- [7. Output Files](#7-output-files)
- [8. Notes](#8-notes)
- [9. References](#9-references)


## 1. Data

**Source:** GEO accession [GSE164378](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164378) ^<br>
**Files:** `GSM5008737 PBMC CITE-seq RNA_3P {barcodes,features,matrix}` (RNA: 33,538 genes × 161,764 cells, ~1.05 GB) and `GSM5008738 PBMC CITE-seq ADT_3P {barcodes,features,matrix}` (ADT: 228 antibodies × 161,764 cells, ~100 MB) <br>
**Design:** 161,764 cells across 13 sequencing lanes in 2 donor pools <br>

^ Please manually download the raw data from NCBI GEO GSE164378 samples `GSM5008737` and `GSM5008738`, and place the files in the
`data/` folder.

Barcodes carry a lane prefix identifying the 13 sequencing lanes across
the two donor pools (`L_pool`: lanes `L1`–`L5` (67,090 cells); `E2_pool`: lanes `E2L1`–`E2L8` (94,674 cells))

reference/hto_doublet_calls.csv contains per-cell hashtag oligonucleotide (HTO)-based doublet/singlet calls derived from 
Hao et al. datasets (2021). It provides barcode-level annotations used in the pipeline to identify and 
remove cells classified as HTO Doublets. <br>

reference/KEGG_RIBOSOME.txt contains a list of human RPL and RPS ribosomal gene set obtained from MSigDB. It serves as a reference for identifying ribosome-related expression.

Barcodes carry a lane prefix identifying the 13 sequencing lanes across
the two donor pools (`L_pool`: lanes `L1`–`L5` (67,090 cells); `E2_pool`: lanes `E2L1`–`E2L8` (94,674 cells))



## 2. Computational Methods

`scripts/` contains the following files:

| Script | Description |
|---|---|
| `run_pipeline.sh` | Master script |
| `01_rna_qc_filter.R` | RNA QC |
| `02_adt_qc_filter.R` | ADT QC |
| `03_rna_umap.R` | RNA UMAP |
| `04_adt_umap.R` | ADT UMAP |
| `05_wnn_integration.R` | WNN-based multimodal Integration |

### RNA QC

| Step | Before | After | Removed | Notes |
|---|---|---|---|---|
| Build Seurat object (RNA+ADT) | 161,764 cells | 161,764 cells | — | 33,538 genes, 228 antibodies |
| Cell filter (≥200 genes/cell) | 161,764 cells| 161,764 cells| 0 | All cells retained |
| Gene filter (≥100 cells/gene) | 33,538 genes | 17,808 genes | 15,730 genes | — |
| Mitochondria filter (<20% `percent.mt`) | 161,764 cells | 161,764 cells | 0 | All cells retained |
| Quantile trim (2–98%, per pool) | 161,764 cells | 153,822 cells | 7,942 cells | — |
| Doublet removal (DoubletFinder + HTO, union) | 153,822 cells | **141,852 cells** | 11,970 cells | Doublets were identified and removed independently within each sequencing lane: 11,537; HTO-based detection: 508; overlap: 75 |

DoubletFinder per-lane breakdown:

| Lane | Cells | Doublets | Singlets |
|---|---|---|---|
| E2L1 | 10,895 | 817 | 10,078 |
| E2L2 | 11,531 | 865 | 10,666 |
| E2L3 | 10,789 | 809 | 9,980 |
| E2L4 | 10,867 | 815 | 10,052 |
| E2L5 | 11,536 | 865 | 10,671 |
| E2L6 | 11,478 | 861 | 10,617 |
| E2L7 | 11,448 | 859 | 10,589 |
| E2L8 | 11,485 | 861 | 10,624 |
| L1 | 13,135 | 985 | 12,150 |
| L2 | 12,462 | 935 | 11,527 |
| L3 | 11,601 | 870 | 10,731 |
| L4 | 12,721 | 954 | 11,767 |
| L5 | 13,874 | 1,041 | 12,833 |
| **Total** | **153,822** | **11,537** | **142,285** |

### ADT QC

| Step | Before | After | Removed | Notes |
|---|---|---|---|---|
| Cell filter (≥20 Ab counts/cell) | 141,852 cells | **141,852 cells** | 0 | All cells retained |
| Antibody filter (≥100 cells/Ab) | 228 Ab | **228 Ab** | 0 | All antibodies retained |

### RNA UMAP

The 141,852-cell, 17,808-gene RNA population was normalized with
`LogNormalize` (scale.factor=10000), then restricted to the top 2000
highly variable genes via `FindVariableFeatures` (`vst` method; top
genes by variance included HBB, HBA2, HBA1, CCL4L2, and CCL20). These
HVGs were scaled and reduced to 30 principal components with `RunPCA`
(seed=42), where PC1 alone explained 52.4% of variance. `RunUMAP` was
run on this `pca` reduction (dims 1:30, seed=42) to produce the
`umap.rna` embedding, and `FindNeighbors`+`FindClusters` (algorithm=3,
resolution=0.5, seed=42) grouped cells into 24 clusters. `FindAllMarkers`,
restricted to the 2000 HVGs, identified 7,446 marker rows across those
clusters, which were then manually annotated into 12 broad lineage
categories based on marker signal. The final labeled RNA UMAP plot is
saved to `results/03_rna_umap_broad_labels.png`.

### ADT UMAP

| # | Step | Result |
|---|---|---|
| 18 | `CLR`, margin=2 | 141,852 cells × 228 antibodies |
| 19 | `ScaleData` + `RunPCA`, all 228 antibodies, 30 PCs, seed=42 | PC1 explains 30.5% of variance (flatter than RNA's 52.4%) |
| 20 | `RunUMAP` on `apca`, dims 1:30, seed=42 | `umap.adt` embedding (no plot) |
| 21 | `FindNeighbors`+`FindClusters`, algorithm=3, resolution=0.5, seed=42 | 29 clusters |
| 22 | `FindAllMarkers`, all 228 antibodies | 1,290 marker rows |
| 23 | Manual marker-based annotation | 9 broad categories |
| 24 | Labeled ADT UMAP plot | `results/04_adt_umap_broad_labels.png` |

### WNN-based multimodal integration

| # | Step | Result |
|---|---|---|
| 25 | `FindMultiModalNeighbors`, `pca`+`apca`, dims 1:30 each, seed=42 | Joint `wknn`/`wsnn` graphs. RNA weight: median 0.55, ranges near-0 to 1.0 per cell |
| 26 | `RunUMAP` on the WNN graph, seed=42 | `wnn.umap` embedding (no plot) |
| 27 | `FindClusters` on `wsnn`, algorithm=3, resolution=1.2, seed=42 | 49 clusters (`wnn_clusters`, kept separate from RNA's `seurat_clusters`) |
| 28 | `FindAllMarkers`, both RNA (2000 HVGs) and ADT (228 antibodies) | 13,194 RNA + 2,320 ADT marker rows |
| 29 | Manual annotation using RNA+ADT markers together | 13 broad categories |
| 30 | Labeled WNN UMAP plot (broad) | `results/05_wnn_umap_broad_labels.png` |
| 31 | Per-cluster detailed annotation, all 49 clusters | `results/05_wnn_detailed_annotation.csv` |
| 32 | Labeled WNN UMAP plot (detailed, numbered + grouped legend) | `results/05_wnn_umap_detailed_labels.png` |

### Methodology notes

- **RNA and ADT are normalized independently** (`LogNormalize` vs. `CLR`)
  since they have fundamentally different statistical properties, then
  combined via WNN integration (script 25) — not at normalization.
- **ADT has no HVG-selection-equivalent step.** With only 228 antibodies
  total — a panel deliberately curated to be informative, unlike RNA's
  whole transcriptome — all of them are used directly as PCA input.
- **Doublet detection runs per lane, not per pool or globally.** A doublet
  physically forms within one droplet-generation run (one lane); pooling
  by sample-pool (~95k–127k cells) OOM-crashed twice in the reference
  project — per-lane (~11k–15k cells) is the scale that works.
- **Both doublet-detection methods are run and cross-checked, not just
  one.** DoubletFinder and HTO-based calls are structurally blind to
  different failure modes; they overlap on only 75 of 11,970 flagged
  cells. Cells are removed if flagged by **either** method — a true
  union (logical OR), not an intersection:
  ```r
  obj$is_doublet <- (obj$doublet_finder_call == "Doublet") |
                     (!is.na(obj$hto_doublet_call) & obj$hto_doublet_call == "Doublet")
  ```
  This choice matters: an intersection (only removing cells both
  methods agreed on) would have kept 11,895 of the 11,970 flagged cells
  in the data — union is what actually lets each method catch what the
  other one structurally can't. HTO-based calls come from
  `reference/hto_doublet_calls.csv` — a
  small, portable (~5 MB), git-friendly barcode → Doublet/Singlet
  lookup, extracted once from a 2 GB intermediate object (this repo
  doesn't implement its own HTO demultiplexing). See "Data" and
  "References" above/below for exactly where these calls originate.
- **Quantile trim is computed per pool**, matching the real, confirmed
  batch effect between `L_pool`/`E2_pool` (Wilcoxon `p < 2.2e-16`), and
  must run *after* doublet removal — percentiles are recomputed from
  whichever cells are currently present.
- **WNN clustering uses a distinct metadata column** (`wnn_clusters`),
  not the default `seurat_clusters` — since the same object already
  carried RNA-only clusters from script 14, letting `FindClusters` use
  its default column name would have silently overwritten them (same
  precaution applied to ADT's `adt_clusters` in script 21).
- **Combined RNA+ADT marker evidence at the WNN stage resolved several
  calls that were ambiguous with either modality alone** — see the WNN
  section above. This is the concrete benefit of doing single-modality
  analysis first and integration second, rather than jumping straight to
  WNN: it makes the value of integration visible and checkable, not just
  assumed.

- Single project-wide seed, **42**, used everywhere a seed is needed:
`set.seed(42)` at the top of scripts with stochastic steps, and
`seed.use = 42` explicitly passed to any Seurat function that accepts it
(`RunPCA`, `RunUMAP`, `FindClusters`'s `random.seed`).



## 3. Environment Setup

Dependencies are pinned in `renv.lock` (147 packages, R version: 4.5.3).  From the project root, open R:

```
install.package("renv")  # if not already installed
renv::restore()
```

Or from the terminal:

```
Rscript -e 'renv::restore()'
```

renv::restore() recreates the project's package environment using the versions 
recorded in `renv.lock` without modifying global R package library.


Key package versions:

| # | Package | Version | Source |
|---|---|---|---|
| 1 | `r-base` (R) | 4.5.3 | conda-forge |
| 2 | `Seurat` | 5.5.1 | conda-forge |
| 3 | `Matrix` | 1.7.5 | conda-forge |
| 4 | `DoubletFinder` | 2.0.6 | GitHub (`chris-mcginnis-ucsf/DoubletFinder`) |


## 4. Analysis Workflow


```
       RNA Matrix                      ADT Matrix
           │                               │
           └───────────────┬───────────────┘
                           │
                           ▼
                Combined Seurat Object
                           │
                           ▼
                        RNA QC
   (Cell/gene/mitochondria/quantile filter, doublet detection)
                           │
           ┌───────────────┴───────────────┐
           │                               │
           ▼                               ▼
     RNA Processing                     ADT QC
 (Log Normalization, HVG,       (Cell/antibody filter)
     scale, PCA, UMAP)                     │
           │                               │
           │                               ▼
           │                        ADT Processing
           │                  (CLR normalization, scale,
           │                          PCA, UMAP)
           │                               │
           └───────────────┬───────────────┘
                           │
                           ▼
                    WNN Integration
   (FindMultiModalNeighbors on RNA PCA + ADT PCA, UMAP)
```


## 5. Running the Analysis Workflow

After completing setup per "Environment Setup" above (raw data already
in `data/`), the pipeline can be run in two ways.

### A. Manual Execution

Run `scripts/run_pipeline.sh`.

```bash
cd multiomics/   # run from the project root

./scripts/run_pipeline.sh
```

Or run each of the scripts in order:

```bash
cd multiomics/   # run from the project root

Rscript scripts/01_rna_qc_filter.R
Rscript scripts/02_adt_qc_filter.R
Rscript scripts/03_rna_umap.R
Rscript scripts/04_adt_umap.R
Rscript scripts/05_wnn_integration.R
```

### B. Snakemake

Snakemake requires a separate environment from the R project environment
managed by `renv`. To create the environment:

```bash
mamba create -n snakemake_env -c bioconda -c conda-forge \
  --no-channel-priority snakemake-minimal=9.23.1
```

A `Snakefile` at the project root automates execution of the pipeline.  Run from the project root with the `snakemake_env` environment active:

```bash
conda activate snakemake_env

snakemake -n --cores 1                                   # dry run
snakemake --cores 2 --config rscript=/full/path/to/Rscript
```


## 6. Directory Structure

```
multiomics/
├── README.md
├── report.html                     # interactive report
├── Snakefile                       # optional Snakemake workflow
├── .Rprofile                       # activates renv for this project
├── renv.lock                       # 147 packages, exact versions, R 4.5.3
├── renv/                           # renv infrastructure
│   ├── activate.R
│   ├── settings.json                
│   └── .gitignore                   
├── data/                            
│   ├── GSM5008737 RNA_3P {barcodes,features,matrix}^  # manually download
│   └── GSM5008738 ADT_3P {barcodes,features,matrix}^  # manually download
├── reference/
│   ├── KEGG_RIBOSOME.txt          # ribosomal gene reference list
│   └── hto_doublet_calls.csv      # HTO-based doublet/singlet calls
├── results/                       
│   ├── 01_rna_qc_report_summary.png
│   ├── 01_rna_qc_filter_cell.csv
│   ├── 01_rna_qc_filter_gene.csv
│   ├── 01_rna_qc_filter_mito.csv
│   ├── 01_rna_qc_filter_quantile.csv       
│   ├── 01_rna_doublet_summary.csv          
│   ├── 02_adt_qc_filter.csv                
│   ├── 02_adt_qc_report_summary.png
│   ├── 03_rna_cluster_markers.csv
│   ├── 03_rna_broad_annotation.csv         
│   ├── 03_rna_umap_broad_labels.png
│   ├── 04_adt_cluster_markers.csv
│   ├── 04_adt_broad_annotation.csv         
│   ├── 04_adt_umap_broad_labels.png
│   ├── 05_wnn_cluster_markers_RNA.csv
│   ├── 05_wnn_cluster_markers_ADT.csv
│   ├── 05_wnn_broad_annotation.csv
│   ├── 05_wnn_umap_broad_labels.png
│   ├── 05_wnn_detailed_annotation.csv      
│   └── 05_wnn_umap_detailed_labels.png
└── scripts/
    ├── run_pipeline.sh            (master script)
    ├── 01_rna_qc_filter.R      
    ├── 02_adt_qc_filter.R      
    ├── 03_rna_umap.R           
    ├── 04_adt_umap.R           
    └── 05_wnn_integration.R    
```


## 7. Output Files

`results/` contains the following files:

- **`01_rna_*`** — RNA cell/gene/mitochondria/quantile filter (`.csv`), doublet detection filter (`.csv`), and RNA QC figure (`.png`)
- **`02_adt_*`** — ADT cell/antibody filter (`.csv`) and ADT QC figure (`.png`)
- **`03_rna_*`** — RNA cluster markers (`.csv`), broad cluster annotation (`.csv`), and RNA UMAP figure (`.png`)
- **`04_adt_*`** — ADT cluster markers (`.csv`), broad cluster annotation (`.csv`), and ADT UMAP figure (`.png`)
- **`05_wnn_*`** — RNA and ADT cluster markers (`.csv`), broad and detailed cluster annotation (`.csv`), and WNN UMAP figures (`.png`)


## 8. Notes

Code was developed with assistance from Claude, an AI coding assistant, 
based on author-defined step-by-step specifications, analytical objectives, and methodological 
decisions. The generated code was iteratively refined, reviewed by the author, and 
validated by reproducing and comparing with major results reported in the original publication.


## 9. References

Hao, Y., Hao, S., Andersen-Nissen, E., Mauck, W. M., Zheng, S., Butler,
A., Lee, M. J., Wilk, A. J., Darby, C., Zager, M., Hoffman, P.,
Stoeckius, M., Papalexi, E., Mimitou, E. P., Jain, J., Srivastava, A.,
Stuart, T., Fleming, L. B., Yeung, B., Rogers, A. J., McElrath, M. J.,
Blish, C. A., Gottardo, R., Smibert, P., & Satija, R. (2021). Integrated
analysis of multimodal single-cell data. *Cell*, *184*(13), 3573–3587.e29.
https://doi.org/10.1016/j.cell.2021.04.048
