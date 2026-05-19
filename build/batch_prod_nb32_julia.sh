#!/bin/bash -l
# Production Julia: N_BASIS=32, SCAN_N=20, 10 nodes, 20 workers (match Fortran prod).
#SBATCH -A m3739
#SBATCH -q premium
#SBATCH -N 10
#SBATCH -t 06:00:00
#SBATCH -C cpu
#SBATCH -J TJLFEP_prod32
#SBATCH -o prod_nb32_julia_%j.out
#SBATCH -e prod_nb32_julia_%j.err
#SBATCH --ntasks=20
#SBATCH --cpus-per-task=64

set -euo pipefail

module load julia/1.11.7
export JULIA_DEPOT_PATH="${PSCRATCH}/.julia"
export SCAN_N=20
export JULIA_WORKER_THREADS="${JULIA_WORKER_THREADS:-${SLURM_CPUS_PER_TASK:-64}}"
export TJLFEP_DEBUG=0
export TJLFEP_FILE_ONLY=1

TJLFEP_ROOT="${TJLFEP_ROOT:-/pscratch/sd/t/tneiser/.julia/dev/TJLFEP}"
export CASE_DIR="${CASE_DIR:-${TJLFEP_ROOT}/src/DIIIDfiles/202017C42_500ms_v3.1}"
export TGLFEP_FILE="${TJLFEP_ROOT}/build/debug_prod/input.TGLFEP"
export FILE_DIR="${FILE_DIR:-${TJLFEP_ROOT}/build/fileInput_prod_${SLURM_JOB_ID:-local}}"

cd "${TJLFEP_ROOT}/build"
run_julia() { stdbuf -oL -eL julia "$@"; }

echo "=== TJLFEP prod: N_BASIS=32 SCAN_N=${SCAN_N} ==="
echo "nodes=${SLURM_NNODES:-?} SLURM_NTASKS=${SLURM_NTASKS:-?} JULIA_WORKER_THREADS=${JULIA_WORKER_THREADS:-?}"
echo "TGLFEP_FILE=${TGLFEP_FILE}"
echo "FILE_DIR=${FILE_DIR}"

run_julia --project="${TJLFEP_ROOT}" validate_nosysimage_distributed.jl

echo "=== done ==="
