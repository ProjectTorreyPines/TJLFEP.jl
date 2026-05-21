#!/bin/bash -l
# nb6 SCAN_N=1 from input.gacode + input.TGLFEP only (no MTGLF/EXPRO/dump).
#SBATCH -A m3739
#SBATCH -q premium
#SBATCH -N 1
#SBATCH -t 01:00:00
#SBATCH -C cpu
#SBATCH -J TJLFEP_gacode_nb6
#SBATCH -o run_gacode_nb6_%j.out
#SBATCH -e run_gacode_nb6_%j.err
#SBATCH --cpus-per-task=8

set -euo pipefail

module load julia/1.11.7
export JULIA_DEPOT_PATH="${PSCRATCH}/.julia${JULIA_DEPOT_PATH:+:${JULIA_DEPOT_PATH}}"
export TJLFEP_FILE_ONLY=1
export TJLFEP_DEBUG=0

TJLFEP_ROOT="${TJLFEP_ROOT:-/pscratch/sd/t/tneiser/.julia/dev/TJLFEP}"
cd "${TJLFEP_ROOT}/build"

echo "=== runTHD_from_gacode nb6 (input.gacode + input.TGLFEP only) ==="
stdbuf -oL -eL julia --startup-file=no --project="${TJLFEP_ROOT}" -t "${SLURM_CPUS_PER_TASK:-8}" run_gacode_nb6.jl
echo "=== done ==="
