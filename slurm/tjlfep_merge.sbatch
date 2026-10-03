#!/bin/bash -l
# =============================================================================
# Re-merge a GPU scan's per-radius outputs (OUT_DIR/task_<i>.jls) into sfmin_scan.txt +
# alpha_*_crit.input. The GPU templates do this automatically at the end of the job; use this
# only if that step failed (e.g. the job hit its wall during the merge) or after re-running
# a single missing radius. Takes under a minute; a login node works too:
#
#   module load julia/1.11.7
#   OUT_DIR=runs/gpu_backfill_<jobid>_tasks GACODE_FILE=... TGLFEP_FILE=... \
#     julia --project=. slurm/common/merge_gacode_scan20_array.jl
#
# Usage (from the repo root):  OUT_DIR=runs/<run>_tasks sbatch slurm/tjlfep_merge.sbatch
# =============================================================================
#SBATCH -A m3739                   # EDIT: your NERSC CPU allocation
#SBATCH -q debug
#SBATCH -C cpu
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH -t 00:15:00
#SBATCH -J tjlfep_merge
#SBATCH -o tjlfep_merge_%j.out
#SBATCH -e tjlfep_merge_%j.err

set -uo pipefail

# ---- CONFIG ----
TJLFEP_ROOT="${TJLFEP_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
CFS_DIR="${CFS_DIR:-/global/cfs/cdirs/m3739/TJLFEP}"
TJLFEP_DEPOT="${TJLFEP_DEPOT:-}"
OUT_DIR="${OUT_DIR:?set OUT_DIR to the <run>_tasks directory holding task_*.jls}"
CASE_DIR="${CASE_DIR:-${TJLFEP_ROOT}/examples/DIIID_202017C42_500ms_v3.1}"
GACODE_FILE="${GACODE_FILE:-${CASE_DIR}/input.gacode}"           # MUST be the same inputs the scan used
TGLFEP_FILE="${TGLFEP_FILE:-${CASE_DIR}/input_scan20_nb32.TGLFEP}"
SYSIMAGE="${SYSIMAGE-${TJLFEP_SYSIMAGE:-}}"                                  # "" = JIT (fine for a merge); or $CFS_DIR/TJLFEP_cpu_sysimage.so
PROJECT="${PROJECT:-}"
# ---- end CONFIG ----

DEVICE=cpu
EXPECT_NTASKS_EQ_SCAN_N=0
source "${TJLFEP_ROOT}/slurm/common/tjlfep_env.inc.sh" \
    || { echo "ERROR: cannot source ${TJLFEP_ROOT}/slurm/common/tjlfep_env.inc.sh -- submit from the repo root or export TJLFEP_ROOT" >&2; exit 1; }
tjlfep_setup
cd "${TJLFEP_ROOT}"
tjlfep_merge
echo "=== merged -> ${OUT_DIR}/sfmin_scan.txt ==="
