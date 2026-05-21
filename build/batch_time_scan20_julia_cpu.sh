#!/bin/bash -l
# Timing: Julia CPU SCAN_N=20, 10 nodes, 20 workers (2/node), SlurmClusterManager + pmap.
# Same layout as batch_debug_nb6_julia_scan20_10n.sh (faithful production path).
#SBATCH -A m3739
#SBATCH -q premium
#SBATCH -N 10
#SBATCH -t 02:00:00
#SBATCH -C cpu
#SBATCH -J time_s20_jcpu
#SBATCH -o time_scan20_julia_cpu_%j.out
#SBATCH -e time_scan20_julia_cpu_%j.err
#SBATCH --ntasks=20
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=64

set -euo pipefail

module load julia/1.11.7
export JULIA_DEPOT_PATH="${PSCRATCH}/.julia${JULIA_DEPOT_PATH:+:${JULIA_DEPOT_PATH}}"
export TJLFEP_FILE_ONLY=1
export TJLFEP_DEBUG=0
export SCAN_N=20
export JULIA_WORKER_THREADS="${JULIA_WORKER_THREADS:-${SLURM_CPUS_PER_TASK:-64}}"

TJLFEP_ROOT="${TJLFEP_ROOT:-/pscratch/sd/t/tneiser/.julia/dev/TJLFEP}"
export CASE_DIR="${CASE_DIR:-${TJLFEP_ROOT}/src/DIIIDfiles/202017C42_500ms_v3.1}"
export GACODE_DUMP="${GACODE_DUMP:-${CASE_DIR}/input.gacode}"
export TGLFEP_FILE="${TJLFEP_ROOT}/build/debug_nb6/input_scan20.TGLFEP"
export FILE_DIR="${TJLFEP_ROOT}/build/debug_nb6/fileInput_scan20_10n_${SLURM_JOB_ID}"

cd "${TJLFEP_ROOT}/build"
# shellcheck source=julia_sysimage.inc.sh
source "${TJLFEP_ROOT}/build/julia_sysimage.inc.sh"

echo "=== Julia CPU timing: SCAN_N=${SCAN_N} nodes=${SLURM_NNODES:-?} ntasks=${SLURM_NTASKS:-?} (SlurmClusterManager) ==="
echo "FILE_DIR=${FILE_DIR}"

stdbuf -oL -eL julia --startup-file=no --project="${TJLFEP_ROOT}" \
    "${JULIA_SYSIMAGE_ARGS[@]}" \
    -t 8 time_scan20_julia_cpu.jl

echo "=== Julia CPU distributed timing finished ==="
