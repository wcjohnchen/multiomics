# Snakemake workflow for the CITE-seq pipeline.
# Alternative to scripts/run_pipeline.sh -- wraps the same 5 scripts,
# but gets Snakemake's DAG-based dependency tracking, partial re-runs,
# and parallelism (rule 01 and rule 02 don't depend on each other and
# can run with --cores 2+).
#
# Usage:
#   snakemake -n                                    # dry run
#   snakemake --cores 1                              # real run, auto-detects package source
#   snakemake --cores 2 --config rscript=/path/to/some/other/Rscript
#
# Auto-detects the same two package-source modes as run_pipeline.sh, and
# for the same reason: Snakemake's shell commands run with cwd = the
# project root by default, which has .Rprofile, so renv would otherwise
# activate on every rule and (since renv/settings.json's
# external.libraries is intentionally empty) hide a conda env's packages
# entirely -- confirmed by actually hitting that exact failure once.
#
#   - conda env 'citeseq-pipeline' found under ~/miniforge3, ~/miniconda3,
#     or ~/anaconda3: each rule cd's into scripts/ first (no .Rprofile
#     there), so renv never activates and packages come straight from
#     that conda env's own Rscript -- no need to activate it yourself
#     first, this points RSCRIPT at it directly.
#   - no matching conda env found: rules run from the project root
#     instead (Snakemake's default), so renv's .Rprofile activates and
#     packages come from renv/library/, i.e. whatever
#     `Rscript -e 'renv::restore()'` installed there beforehand.
#
# Override the auto-detected Rscript with --config rscript=... if needed
# (e.g. a conda env at a nonstandard path) -- this only overrides which
# Rscript binary is used, not which of the two cwd modes runs; that's
# still decided by whether the 'citeseq-pipeline' env directory exists.
#
# Data prerequisite (not automated here, same as run_pipeline.sh): the 6
# raw GEO files must already be in data/ -- see README's "Source data".

import os


def _find_conda_env(name):
    home = os.path.expanduser("~")
    for base in ("miniforge3", "miniconda3", "anaconda3"):
        candidate = os.path.join(home, base, "envs", name)
        if os.path.isdir(candidate):
            return candidate
    return None


_conda_env_path = _find_conda_env("citeseq-pipeline")

if _conda_env_path is not None:
    RSCRIPT = config.get("rscript", os.path.join(_conda_env_path, "bin", "Rscript"))
    RUN_PREFIX = "cd scripts && "
    SCRIPT_DIR = ""
else:
    RSCRIPT = config.get("rscript", "Rscript")
    RUN_PREFIX = ""
    SCRIPT_DIR = "scripts/"

RAW_RNA = expand("data/GSM5008737_RNA_3P-{part}.{ext}",
                  part=["barcodes", "features"], ext=["tsv.gz"]) + \
          ["data/GSM5008737_RNA_3P-matrix.mtx.gz"]
RAW_ADT = expand("data/GSM5008738_ADT_3P-{part}.{ext}",
                  part=["barcodes", "features"], ext=["tsv.gz"]) + \
          ["data/GSM5008738_ADT_3P-matrix.mtx.gz"]

rule all:
    input:
        "results/01_rna_qc_report_summary.png",
        "results/02_adt_qc_report_summary.png",
        "results/03_rna_umap_broad_labels.png",
        "results/04_adt_umap_broad_labels.png",
        "results/05_wnn_umap_broad_labels.png",
        "results/05_wnn_umap_detailed_labels.png"

rule rna_qc_filter:
    input:
        raw_rna = RAW_RNA,
        raw_adt = RAW_ADT,
        ribo_ref = "reference/KEGG_RIBOSOME.txt",
        hto_ref = "reference/hto_doublet_calls.csv",
        script = "scripts/01_rna_qc_filter.R"
    output:
        raw_obj = "results/01_rna_seurat_object_raw.rds",
        qc_fig = "results/01_rna_qc_report_summary.png",
        quantile_trim_obj = "results/01_rna_seurat_object_quantile_trim.rds",
        doublet_obj = "results/01_rna_seurat_object_doublet_filtered.rds",
        cell_funnel = "results/01_rna_qc_filter_cell.csv",
        gene_funnel = "results/01_rna_qc_filter_gene.csv",
        mito_funnel = "results/01_rna_qc_filter_mito.csv",
        quantile_funnel = "results/01_rna_qc_filter_quantile.csv",
        doublet_funnel = "results/01_rna_doublet_summary.csv"
    shell:
        RUN_PREFIX + RSCRIPT + " " + SCRIPT_DIR + "01_rna_qc_filter.R"

rule adt_qc_filter:
    input:
        raw_obj = "results/01_rna_seurat_object_raw.rds",
        quantile_trim_obj = "results/01_rna_seurat_object_quantile_trim.rds",
        script = "scripts/02_adt_qc_filter.R"
    output:
        adt_obj = "results/02_adt_seurat_object_filtered.rds",
        qc_fig = "results/02_adt_qc_report_summary.png",
        funnel = "results/02_adt_qc_filter.csv"
    shell:
        RUN_PREFIX + RSCRIPT + " " + SCRIPT_DIR + "02_adt_qc_filter.R"

rule rna_umap:
    input:
        doublet_obj = "results/01_rna_seurat_object_doublet_filtered.rds",
        script = "scripts/03_rna_umap.R"
    output:
        annotated_obj = "results/03_rna_seurat_object_annotated.rds",
        markers = "results/03_rna_cluster_markers.csv",
        annotation = "results/03_rna_broad_annotation.csv",
        umap_fig = "results/03_rna_umap_broad_labels.png"
    shell:
        RUN_PREFIX + RSCRIPT + " " + SCRIPT_DIR + "03_rna_umap.R"

rule adt_umap:
    input:
        rna_annotated_obj = "results/03_rna_seurat_object_annotated.rds",
        script = "scripts/04_adt_umap.R"
    output:
        adt_annotated_obj = "results/04_adt_seurat_object_annotated.rds",
        markers = "results/04_adt_cluster_markers.csv",
        annotation = "results/04_adt_broad_annotation.csv",
        umap_fig = "results/04_adt_umap_broad_labels.png"
    shell:
        RUN_PREFIX + RSCRIPT + " " + SCRIPT_DIR + "04_adt_umap.R"

rule wnn_integration:
    input:
        adt_annotated_obj = "results/04_adt_seurat_object_annotated.rds",
        script = "scripts/05_wnn_integration.R"
    output:
        wnn_obj = "results/05_wnn_seurat_object_detailed.rds",
        rna_markers = "results/05_wnn_cluster_markers_RNA.csv",
        adt_markers = "results/05_wnn_cluster_markers_ADT.csv",
        broad_annotation = "results/05_wnn_broad_annotation.csv",
        finest_annotation = "results/05_wnn_detailed_annotation.csv",
        broad_fig = "results/05_wnn_umap_broad_labels.png",
        finest_fig = "results/05_wnn_umap_detailed_labels.png"
    shell:
        RUN_PREFIX + RSCRIPT + " " + SCRIPT_DIR + "05_wnn_integration.R"
