#!/bin/bash -l
# Julia: N_BASIS=32, SCAN_N=20, 10 nodes, 20 workers (SlurmClusterManager).
#SBATCH -A m3739
#SBATCH -q premium
#SBATCH -N 10
#SBATCH -t 06:00:00
#SBATCH -C cpu
#SBATCH -J TJLFEP_nb32_s20
#SBATCH -o debug_nb32_julia20_10n_%j.out
#SBATCH -e debug_nb32_julia20_10n_%j.err
#SBATCH --ntasks=20
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=64

set -euo pipefail

module load julia/1.11.7
export JULIA_DEPOT_PATH="${PSCRATCH}/.julia${JULIA_DEPOT_PATH:+:${JULIA_DEPOT_PATH}}"
export TJLFEP_FILE_ONLY=1
export SCAN_N=20
export JULIA_WORKER_THREADS="${JULIA_WORKER_THREADS:-${SLURM_CPUS_PER_TASK:-64}}"
export TJLFEP_DEBUG=0

TJLFEP_ROOT="${TJLFEP_ROOT:-/pscratch/sd/t/tneiser/.julia/dev/TJLFEP}"
export CASE_DIR="${CASE_DIR:-${TJLFEP_ROOT}/src/DIIIDfiles/202017C42_500ms_v3.1}"
export GACODE_DUMP="${GACODE_DUMP:-${CASE_DIR}/input.gacode}"
export TGLFEP_FILE="${TJLFEP_ROOT}/build/debug_nb32/input_scan20.TGLFEP"
export FILE_DIR="${TJLFEP_ROOT}/build/debug_nb32/fileInput_scan20_10n_${SLURM_JOB_ID:-local}"

cd "${TJLFEP_ROOT}/build"
# shellcheck source=julia_sysimage.inc.sh
source "${TJLFEP_ROOT}/build/julia_sysimage.inc.sh"

echo "=== Julia nb32 SCAN_N=20 on ${SLURM_NNODES:-10} nodes ==="
stdbuf -oL -eL julia --startup-file=no --project="${TJLFEP_ROOT}" \
    "${JULIA_SYSIMAGE_ARGS[@]}" \
    debug_compare_nb32_scan20_distributed.jl
echo "=== Julia nb32 scan20 distributed finished ==="
