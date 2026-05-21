#!/bin/bash -l
# Submit SCAN_N=20 GPU scan + in-job merge (m3739_g, premium).
set -euo pipefail

TJLFEP_ROOT="${TJLFEP_ROOT:-/pscratch/sd/t/tneiser/.julia/dev/TJLFEP}"
cd "${TJLFEP_ROOT}/build"

SCAN_JOB=$(sbatch --parsable batch_run_gacode_scan20_gpu_5nodes.sh)
echo "Submitted 5-node GPU scan+merge job ${SCAN_JOB}"
echo "OUT_DIR=${TJLFEP_ROOT}/build/gacode_scan20_${SCAN_JOB}_tasks"
echo "Monitor: squeue -j ${SCAN_JOB}"
