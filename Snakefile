# Snakemake workflow for the CITE-seq pipeline.
# Alternative to scripts/run_pipeline.sh -- wraps the same 5 scripts,
# but gets Snakemake's DAG-based dependency tracking, partial re-runs,
# and parallelism (rule rna_umap and rule adt_qc_filter/adt_umap don't
# depend on each other -- both branch independently from rna_qc_filter's
# doublet-filtered checkpoint, only merging back at wnn_integration --
# and can run with --cores 2+).
#
# Usage:
#   snakemake -n                                    # dry run
#   snakemake --cores 1                              # real run, auto-detects package source
#   snakemake --cores 2 --config conda_env=my_env_name
#   snakemake --cores 2 --config rscript=/path/to/some/other/Rscript
#
# Auto-detects the same two package-source modes as run_pipeline.sh, and
# for the same reason: Snakemake's shell commands run with cwd = the
# project root by default, which has .Rprofile, so renv would otherwise
# activate on every rule and (since renv/settings.json's
# external.libraries is intentionally empty) hide a conda env's packages
# entirely -- confirmed by actually hitting that exact failure once.
#
#   - a conda env is found (default name 'citeseq-pipeline', overridable
#     with --config conda_env=..., since a given user's env doesn't have
#     to be named that) via whatever `conda` is on PATH -- asks it for its
#     own base dir with `conda info --base` rather than guessing
#     install-dir names, so this works regardless of which conda
#     distribution or install location/prefix is on this machine: each
#     rule cd's into scripts/ first (no .Rprofile there), so renv never
#     activates and packages come straight from that conda env's own
#     Rscript -- no need to activate it yourself first, this points
#     RSCRIPT at it directly.
#   - no `conda` on PATH, or no matching env under its base: rules run
#     from the project root instead (Snakemake's default), so renv's
#     .Rprofile activates and packages come from renv/library/, i.e.
#     whatever `Rscript -e 'renv::restore()'` installed there beforehand.
#
# Override the auto-detected Rscript directly with --config rscript=...
# if needed (e.g. a conda env at a nonstandard path not findable via its
# name at all) -- this only overrides which Rscript binary is used, not
# which of the two cwd modes runs; that's still decided by whether the
# conda env directory (whatever its name) was found.
#
# Data prerequisite (not automated here, same as run_pipeline.sh): the 6
# raw GEO files must already be in data/ -- see README's "Source data".

import os
import shutil
import subprocess


def _find_conda_env(name):
    # Prefer asking whatever `conda` is already on PATH for its own base
    # dir -- works for any conda distribution or install location/prefix.
    # But Snakemake may be invoked without conda's shell hooks active (e.g.
    # a non-interactive shell, or a wrapper that doesn't run `conda init`'s
    # PATH additions), so fall back to the common install locations if
    # `conda` isn't found on PATH -- same two-tier approach as
    # run_pipeline.sh, kept consistent for the same reason.
    if shutil.which("conda") is not None:
        try:
            base = subprocess.run(
                ["conda", "info", "--base"],
                capture_output=True, text=True, check=True, timeout=10,
            ).stdout.strip()
        except (subprocess.SubprocessError, OSError):
            base = None
        if base:
            candidate = os.path.join(base, "envs", name)
            if os.path.isdir(candidate):
                return candidate

    home = os.path.expanduser("~")
    for base_name in ("miniforge3", "mambaforge", "miniconda3", "anaconda3"):
        candidate = os.path.join(home, base_name, "envs", name)
        if os.path.isdir(candidate):
            return candidate
    return None


_conda_env_path = _find_conda_env(config.get("conda_env", "citeseq-pipeline"))

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
        "results/05_wnn_umap_detailed_labels.png",
        "results/05_wnn_umap_broad_labels.png"

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
        cell_filter = "results/01_rna_qc_filter_cell.csv",
        gene_filter = "results/01_rna_qc_filter_gene.csv",
        mito_filter = "results/01_rna_qc_filter_mito.csv",
        quantile_filter = "results/01_rna_qc_filter_quantile.csv",
        doublet_summary = "results/01_rna_doublet_summary.csv"
    shell:
        RUN_PREFIX + RSCRIPT + " " + SCRIPT_DIR + "01_rna_qc_filter.R"

rule adt_qc_filter:
    input:
        raw_obj = "results/01_rna_seurat_object_raw.rds",
        doublet_obj = "results/01_rna_seurat_object_doublet_filtered.rds",
        script = "scripts/02_adt_qc_filter.R"
    output:
        adt_obj = "results/02_adt_seurat_object_filtered.rds",
        qc_fig = "results/02_adt_qc_report_summary.png",
        qc_filter = "results/02_adt_qc_filter.csv"
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
        adt_filtered_obj = "results/02_adt_seurat_object_filtered.rds",
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
        rna_annotated_obj = "results/03_rna_seurat_object_annotated.rds",
        adt_annotated_obj = "results/04_adt_seurat_object_annotated.rds",
        script = "scripts/05_wnn_integration.R"
    output:
        wnn_obj = "results/05_wnn_seurat_object_detailed.rds",
        rna_markers = "results/05_wnn_cluster_markers_RNA.csv",
        adt_markers = "results/05_wnn_cluster_markers_ADT.csv",
        detailed_annotation = "results/05_wnn_detailed_annotation.csv",
        detailed_fig = "results/05_wnn_umap_detailed_labels.png",
        broad_annotation = "results/05_wnn_broad_annotation.csv",
        broad_fig = "results/05_wnn_umap_broad_labels.png"
    shell:
        RUN_PREFIX + RSCRIPT + " " + SCRIPT_DIR + "05_wnn_integration.R"
