#!/usr/bin/env bash
#
# Master script: runs the full CITE-seq pipeline end-to-end, in order:
#   01_rna_qc_filter.R      (steps 1-6, 9:  build object -> RNA QC -> filter -> doublet removal)
#   02_adt_qc_filter.R      (steps 7-8:     ADT QC -> filter)
#   03_rna_umap.R           (steps 10-17:   RNA normalize -> UMAP -> cluster -> annotate)
#   04_adt_umap.R           (steps 18-24:   ADT normalize -> UMAP -> cluster -> annotate)
#   05_wnn_integration.R    (steps 25-32:   WNN integrate -> UMAP -> cluster -> annotate)
#
# 02 depends only on 01's raw checkpoint; 03/04/05 form a strict chain.
# Stops immediately on the first failing step (nothing downstream would be
# trustworthy anyway). Each step's full R output goes to its own log file
# under logs/, in addition to the console.
#
# Two package-source modes, auto-detected -- no machine-specific config
# needed either way, so this file and the whole repo stay identical
# whether run here or after a fresh git clone:
#
#   - conda env 'citeseq-pipeline' found: run from scripts/ (no .Rprofile
#     there), so renv never activates and packages come straight from
#     that conda env. This is the mode used on the machine this was
#     built on.
#   - no matching conda env: run from the PROJECT ROOT instead, so renv's
#     .Rprofile activates and packages come from renv/library/, i.e.
#     whatever `Rscript -e 'renv::restore()'` installed there. This is
#     the mode a fresh GitHub clone uses.
#
# (Deliberately NOT done via renv/settings.json's external.libraries --
# a hardcoded path there would crash renv activation entirely on any
# machine where that exact path doesn't exist and can't be created,
# e.g. anyone else's clone. Keeping that setting empty is what makes
# renv::restore() safe there in the first place.)
#
# Either way: each script's own project_dir/data_dir detection is based
# on its own file location, not the working directory, so which mode
# runs doesn't affect where inputs/outputs are read or written.
#
# Expect ~1.5-2.5 hours total (01 alone is ~35-40 min; 05 is ~60 min).
#
# Usage:
#   ./run_pipeline.sh

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
)

log()  { printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*"; }
hr()   { printf -- '--------------------------------------------------------------\n'; }

# --- pick a package-source mode ---------------------------------------------

if [ -f "$HOME/miniforge3/etc/profile.d/conda.sh" ]; then
  # shellcheck disable=SC1091
  source "$HOME/miniforge3/etc/profile.d/conda.sh"
elif [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
  # shellcheck disable=SC1091
  source "$HOME/miniconda3/etc/profile.d/conda.sh"
fi

if command -v conda >/dev/null 2>&1 && conda env list 2>/dev/null | grep -q '^citeseq-pipeline '; then
  conda activate citeseq-pipeline
  RUN_DIR="$SCRIPT_DIR"
  SCRIPT_PREFIX=""
  log "Activated conda env 'citeseq-pipeline' -- running from scripts/ (renv stays inactive)"
elif command -v Rscript >/dev/null 2>&1; then
  RUN_DIR="$PROJECT_DIR"
  SCRIPT_PREFIX="scripts/"
  log "No 'citeseq-pipeline' conda env found -- using renv to supply packages ($(command -v Rscript))"

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
      # A parallel-install race (some packages start building before a
      # dependency like rlang finishes elsewhere) can fail the first
      # attempt -- confirmed reproducible, and confirmed a plain retry
      # fixes it since already-built packages are reused from cache.
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
  echo "ERROR: no conda env 'citeseq-pipeline' and no Rscript on PATH" >&2
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
