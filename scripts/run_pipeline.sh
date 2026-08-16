#!/usr/bin/env bash
#
# Run this master script to run the full CITE-seq pipeline, in order:
#   01_rna_qc_filter.R
#   02_adt_qc_filter.R
#   03_rna_umap.R
#   04_adt_umap.R
#   05_wnn_integration.R
#   06_generate_report.R
#
# Usage:
#   ./run_pipeline.sh
#   CITESEQ_CONDA_ENV=my_env_name ./run_pipeline.sh   # option: use a different conda env name (default: 'citeseq-pipeline')

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
LOG_DIR="$PROJECT_DIR/logs"
mkdir -p "$LOG_DIR"

RUN_STAMP="$(date +%Y%m%d_%H%M%S)"

STEPS=(
  "01_rna_qc_filter.R"
  "02_adt_qc_filter.R"
  "03_rna_umap.R"
  "04_adt_umap.R"
  "05_wnn_integration.R"
  "06_generate_report.R"
)

log()  { printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*"; }
hr()   { printf -- '--------------------------------------------------------------\n'; }

# --- check raw data is present before doing anything else -------------------

DATA_DIR="$PROJECT_DIR/data"
REQUIRED_DATA_FILES=(
  "GSM5008737_RNA_3P-barcodes.tsv.gz"
  "GSM5008737_RNA_3P-features.tsv.gz"
  "GSM5008737_RNA_3P-matrix.mtx.gz"
  "GSM5008738_ADT_3P-barcodes.tsv.gz"
  "GSM5008738_ADT_3P-features.tsv.gz"
  "GSM5008738_ADT_3P-matrix.mtx.gz"
)
MISSING_DATA_FILES=()
for f in "${REQUIRED_DATA_FILES[@]}"; do
  [ -f "$DATA_DIR/$f" ] || MISSING_DATA_FILES+=("$f")
done
if [ "${#MISSING_DATA_FILES[@]}" -gt 0 ]; then
  echo "ERROR: missing raw data file(s) in $DATA_DIR:" >&2
  for f in "${MISSING_DATA_FILES[@]}"; do
    echo "  - $f" >&2
  done
  echo "Download GSM5008737 and GSM5008738 from NCBI GEO GSE164378 and place" >&2
  echo "them in data/ -- see README's \"1. Data\" section." >&2
  exit 1
fi

# --- pick a package-source mode ---------------------------------------------

# Conda env name defaults to 'citeseq-pipeline' but is overridable
#   CITESEQ_CONDA_ENV=my_env_name ./scripts/run_pipeline.sh

CONDA_ENV_NAME="${CITESEQ_CONDA_ENV:-citeseq-pipeline}"

if command -v conda >/dev/null 2>&1; then
  CONDA_BASE="$(conda info --base 2>/dev/null)"
  if [ -n "$CONDA_BASE" ] && [ -f "$CONDA_BASE/etc/profile.d/conda.sh" ]; then
    # shellcheck disable=SC1091
    source "$CONDA_BASE/etc/profile.d/conda.sh"
  fi
else
  for _base in "$HOME/miniforge3" "$HOME/mambaforge" "$HOME/miniconda3" "$HOME/anaconda3"; do
    if [ -f "$_base/etc/profile.d/conda.sh" ]; then
      # shellcheck disable=SC1091
      source "$_base/etc/profile.d/conda.sh"
      break
    fi
  done
  unset _base
fi

if command -v conda >/dev/null 2>&1 && conda env list 2>/dev/null | grep -q "^${CONDA_ENV_NAME} "; then
  conda activate "$CONDA_ENV_NAME"
  RUN_DIR="$SCRIPT_DIR"
  SCRIPT_PREFIX=""
  log "Activated conda env '$CONDA_ENV_NAME' -- running from scripts/ (renv stays inactive)"
elif command -v Rscript >/dev/null 2>&1; then
  RUN_DIR="$PROJECT_DIR"
  SCRIPT_PREFIX="scripts/"
  log "No '$CONDA_ENV_NAME' conda env found -- using renv to supply packages ($(command -v Rscript))"

  cd "$PROJECT_DIR"
  if Rscript -e 'quit(status = if (requireNamespace("Seurat", quietly = TRUE) && requireNamespace("DoubletFinder", quietly = TRUE)) 0L else 1L)' > /dev/null 2>&1; then
    log "Packages already installed (renv) -- skipping restore"
  else
    log "Packages not yet installed -- running renv::restore() automatically"
    log "(compiles from source, can take 20-40+ minutes on first run)"
    restore_log="$LOG_DIR/${RUN_STAMP}_renv_restore.log"

    if Rscript -e 'renv::restore(prompt = FALSE)' > "$restore_log" 2>&1; then
      log "renv::restore() succeeded -- see $(basename "$restore_log")"
    else

      log "First renv::restore() attempt failed -- this can happen due to a"
      log "known parallel-install race condition (see README). Retrying once..."
      if Rscript -e 'renv::restore(prompt = FALSE)' >> "$restore_log" 2>&1; then
        log "renv::restore() succeeded on retry -- see $(basename "$restore_log")"
      else
        log "FAILED: renv::restore() failed twice in a row."
        hr
        log "Last 40 lines of $restore_log:"
        tail -40 "$restore_log"
        hr
        exit 1
      fi
    fi
  fi
else
  echo "ERROR: no conda env '$CONDA_ENV_NAME' and no Rscript on PATH" >&2
  exit 1
fi

# --- run each step in order -------------------------------------------------

cd "$RUN_DIR"

pipeline_start=$(date +%s)
hr
log "Starting full pipeline run ($RUN_STAMP), ${#STEPS[@]} steps"
log "Logs: $LOG_DIR/${RUN_STAMP}_<step>.log"
hr

for script in "${STEPS[@]}"; do
  step_log="$LOG_DIR/${RUN_STAMP}_${script%.R}.log"
  log "Running $script -> $(basename "$step_log")"
  step_start=$(date +%s)

  if Rscript "${SCRIPT_PREFIX}${script}" > "$step_log" 2>&1; then
    step_end=$(date +%s)
    log "OK: $script finished in $(( (step_end - step_start) / 60 ))m $(( (step_end - step_start) % 60 ))s"
  else
    step_end=$(date +%s)
    log "FAILED: $script after $(( (step_end - step_start) / 60 ))m $(( (step_end - step_start) % 60 ))s"
    hr
    log "Last 40 lines of $step_log:"
    tail -40 "$step_log"
    hr
    log "Pipeline stopped. Fix the error above and rerun (each script re-reads its own"
    log "input checkpoint from results/, so completed steps do not need to be redone)."
    exit 1
  fi
  hr
done

pipeline_end=$(date +%s)
total=$(( pipeline_end - pipeline_start ))
log "Pipeline complete: all ${#STEPS[@]} steps finished in $(( total / 3600 ))h $(( (total % 3600) / 60 ))m $(( total % 60 ))s"
log "Final checkpoint: results/05_wnn_seurat_object_detailed.rds"
log "Report: report.html"
