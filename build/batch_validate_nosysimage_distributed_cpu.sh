#!/bin/bash -l
# File-based validation: SlurmClusterManager + pmap (one radius per worker).
# Option 1: 10 nodes, 20 workers (2 radii/node), 64 Julia threads per worker for kw scan.
# Outer parallelism = processes only (no Threads.@threads over radii); inner = @threads in kwscale_scan.
#SBATCH -A m3739
#SBATCH -q premium
#SBATCH -N 10
#SBATCH -t 01:00:00
#SBATCH -C cpu
#SBATCH -J TJLFEP_val_dist
#SBATCH -o validate_nosysimage_dist_%j.out
#SBATCH -e validate_nosysimage_dist_%j.err
#SBATCH --ntasks=20
#SBATCH --cpus-per-task=64

set -euo pipefail

module load julia/1.11.7

export JULIA_DEPOT_PATH="${PSCRATCH}/.julia"
export TJLFEP_FILE_ONLY=1
export SCAN_N=20
# Match Slurm CPUs per worker (set after allocation starts; default 64 for this script).
export JULIA_WORKER_THREADS="${JULIA_WORKER_THREADS:-${SLURM_CPUS_PER_TASK:-64}}"

TJLFEP_ROOT="${TJLFEP_ROOT:-/pscratch/sd/t/tneiser/.julia/dev/TJLFEP}"
export CASE_DIR="${CASE_DIR:-${TJLFEP_ROOT}/src/DIIIDfiles/202017C42_500ms_v3.1}"
export FILE_DIR="${FILE_DIR:-${TJLFEP_ROOT}/build/fileInput_${SLURM_JOB_ID:-local}}"

cd "${TJLFEP_ROOT}/build"

# Slurm redirects .out/.err to files → Julia uses block buffering; line-buffer + flush for live logs.
run_julia() {
  stdbuf -oL -eL julia "$@"
}

echo "=== TJLFEP validate (file-based, distributed) ==="
echo "nodes=${SLURM_NNODES:-?} SLURM_NTASKS=${SLURM_NTASKS:-?} SCAN_N=${SCAN_N} JULIA_WORKER_THREADS=${JULIA_WORKER_THREADS:-?}"
echo "JULIA_DEPOT_PATH=${JULIA_DEPOT_PATH}"

# Single driver: precompile TJLFEP+TJLF on manager, then spawn workers (load pkgimages only).
run_julia --project="${TJLFEP_ROOT}" validate_nosysimage_distributed.jl

echo "=== validate distributed finished ==="
