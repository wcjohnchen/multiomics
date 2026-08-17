# Integrated Analysis of Single-Cell Multiomics Data Using CITE-seq

This project performs multimodal single-cell analysis by integrating RNA expression and antibody-derived tag (ADT) measurements from a publicly available CITE-seq dataset.

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

reference/celltype_annotations.csv contains cell-type labels for the WNN UMAP clusters based on the Hao et al. (2021) paper.


## 2. Computational Methods

`scripts/` contains the following scripts:

| Script | Description |
|---|---|
| `run_pipeline.sh` | Master script |
| `01_rna_qc_filter.R` | RNA QC |
| `02_adt_qc_filter.R` | ADT QC |
| `03_rna_umap.R` | RNA UMAP |
| `04_adt_umap.R` | ADT UMAP |
| `05_wnn_integration.R` | WNN-based multimodal Integration |
| `06_generate_report.R` | Generate report (html) |

### RNA QC

RNA data were filtered based on cell-level, gene-level, mitochondrial-content, and outlier criteria, followed by doublet detection and removal.

| Step | Before | After | Removed | Notes |
|---|---|---|---|---|
| Build Seurat object (RNA+ADT) | 161,764 cells | 161,764 cells | — | 33,538 genes, 228 antibodies |
| Cell filter (≥200 genes/cell) | 161,764 cells| 161,764 cells| 0 | All cells retained |
| Gene filter (≥100 cells/gene) | 33,538 genes | 17,808 genes | 15,730 genes | — |
| Mitochondria filter (<20% `percent.mt`) | 161,764 cells | 161,764 cells | 0 | All cells retained |
| Quantile trim (2–98%, per pool) | 161,764 cells | 153,822 cells | 7,942 cells | — |
| Doublet removal (DoubletFinder + HTO, union) | 153,822 cells | **141,852 cells** | 11,970 cells | Doublets were identified and removed independently within each sequencing lane: 11,537; HTO-based detection: 508; overlap: 75; seed=42 |

Doublet detection is a method for identifying and removing droplets containing two or more cells, which can introduce artifacts into downstream analysis. DoubletFinder per-lane breakdown:

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

ADT data measure cell-surface protein abundance and were filtered based on antibody-level and cell-level criteria.

| Step | Before | After | Removed | Notes |
|---|---|---|---|---|
| Cell filter (≥20 antibodies detected/cell) | 141,852 cells | **141,852 cells** | 0 | All cells retained |
| Antibody filter (≥100 cells/Ab) | 228 Ab | **228 Ab** | 0 | All antibodies retained |

### RNA UMAP

RNA counts were normalized using the log-normalization method (Seurat), and the top 2,000 highly variable genes (HVGs) were selected. PCA was performed for dimensionality reduction, 
and the resulting principal components used for clustering and UMAP visualization.  Cell clusters were annotated based on canonical gene markers.  The same random seed was used across analysis to ensure reproducibility.

| Parameter | Notes |
|---|---|
| `LogNormalize` | scale.factor=10000 |
| `FindVariableFeatures` | method=vst, top 2000 HVGs |
| `RunPCA` | 30 PCs, seed=42 |
| `RunUMAP` | reduction=`pca`, dims=1:30, seed=42 |
| `FindNeighbors` + `FindClusters` | algorithm=3 (SLM), resolution=0.5, seed=42 |
| `FindAllMarkers` | restricted to 2000 HVGs |

### ADT UMAP

ADT counts were normalized using the centered log-ratio (CLR) normalization method (Seurat), followed by PCA, clustering, UMAP, and annotation based on canonical protein markers.

| Parameter | Notes |
|---|---|
| `NormalizeData` | method=CLR, margin=2 |
| `ScaleData` + `RunPCA` | 30 PCs, seed=42 |
| `RunUMAP` | reduction=`apca`, dims=1:30, seed=42 |
| `FindNeighbors` + `FindClusters` | algorithm=3 (SLM), resolution=0.5, seed=42 |
| `FindAllMarkers` | 228 antibodies |

### WNN-based multimodal integration

RNA and ADT principal components were used to construct a weighted nearest neighbors (WNN)-based graph for multimodal integration, followed by UMAP visualization.  Cell clusters were annotated based on the referenece provided by Hao et al. (2021).

| Parameter | Notes |
|---|---|
| `FindMultiModalNeighbors` | `pca`+`apca`, dims 1:30 each, seed=42 |
| `RunUMAP` | reduction=weighted.nn, seed=42 |
| `FindClusters` | graph=`wsnn`, algorithm=3 (SLM), resolution=1.2, seed=42 |
| `FindAllMarkers` | RNA (2000 HVGs) and ADT (228 antibodies) |


## 3. Environment Setup

Dependencies are pinned in `renv.lock` (147 packages, R version: 4.5.3).  From the project root, open R:

```
install.packages("renv")  # if not already installed
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
 (Log normalization, HVG,       (Cell/antibody filter)
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
                           │
                           ▼
                   Generate report.html
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
Rscript scripts/06_generate_report.R
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

snakemake -n --cores 1     # dry run
snakemake --cores 2        # real run
snakemake --cores 2 --config rscript=/full/path/to/Rscript   # optional: override the auto-detected Rscript
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
│   ├── KEGG_RIBOSOME.txt              # ribosomal gene reference list
│   ├── hto_doublet_calls.csv          # HTO-based doublet/singlet calls
│   └── celltype_annotations.csv       # cell type annotation from Hao et al. 2021 paper
├── results/                       
│   ├── 01_rna_qc_report_summary.png
│   ├── 01_rna_qc_filter_cell.csv
│   ├── 01_rna_qc_filter_gene.csv
│   ├── 01_rna_qc_filter_mito.csv
│   ├── 01_rna_qc_filter_quantile.csv       
│   ├── 01_rna_doublet_summary.csv          
│   ├── 01_rna_qc_stats.csv                 
│   ├── 02_adt_qc_filter.csv                
│   ├── 02_adt_qc_report_summary.png
│   ├── 02_adt_qc_stats.csv                 
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
    ├── 05_wnn_integration.R    
    ├── 06_generate_report.R      
    └── report_template.html      
```


## 7. Output Files

`results/` contains the following files:

- **`01_rna_*`** — RNA cell/gene/mitochondria/quantile filter and statistics (`.csv`), doublet detection filter (`.csv`), and RNA QC figure (`.png`)
- **`02_adt_*`** — ADT cell/antibody filter and statistics (`.csv`) and ADT QC figure (`.png`)
- **`03_rna_*`** — RNA cluster markers (`.csv`), broad cluster annotation (`.csv`), and RNA UMAP figure (`.png`)
- **`04_adt_*`** — ADT cluster markers (`.csv`), broad cluster annotation (`.csv`), and ADT UMAP figure (`.png`)
- **`05_wnn_*`** — RNA and ADT cluster markers (`.csv`), broad and detailed cluster annotation (`.csv`), and WNN UMAP figures (`.png`)


## 8. Notes

Code was developed with assistance from Claude, an AI coding assistant,
based on author-defined step-by-step specifications, analytical
objectives, and methodological decisions. The generated code was
iteratively refined, reviewed by the author, and validated by
comparing major results reported in the original publication.


## 9. References

Hao, Y., Hao, S., Andersen-Nissen, E., Mauck, W. M., Zheng, S., Butler,
A., Lee, M. J., Wilk, A. J., Darby, C., Zager, M., Hoffman, P.,
Stoeckius, M., Papalexi, E., Mimitou, E. P., Jain, J., Srivastava, A.,
Stuart, T., Fleming, L. B., Yeung, B., Rogers, A. J., McElrath, M. J.,
Blish, C. A., Gottardo, R., Smibert, P., & Satija, R. (2021). Integrated
analysis of multimodal single-cell data. *Cell*, *184*(13), 3573–3587.e29.
https://doi.org/10.1016/j.cell.2021.04.048
