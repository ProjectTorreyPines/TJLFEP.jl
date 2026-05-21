#!/bin/bash -l
# Submit SCAN_N=20 GPU scan (5 nodes, 20 tasks) + merge job (m3739_g, premium).
set -euo pipefail

TJLFEP_ROOT="${TJLFEP_ROOT:-/pscratch/sd/t/tneiser/.julia/dev/TJLFEP}"
cd "${TJLFEP_ROOT}/build"

SCAN_JOB=$(sbatch --parsable batch_run_gacode_scan20_gpu_5nodes.sh)
echo "Submitted 5-node GPU scan job ${SCAN_JOB}"

export OUT_DIR="${TJLFEP_ROOT}/build/gacode_scan20_${SCAN_JOB}_tasks"
MERGE_JOB=$(sbatch --parsable --dependency=afterok:"${SCAN_JOB}" \
    --export=ALL,OUT_DIR="${OUT_DIR}" \
    batch_merge_gacode_scan20.sh)
echo "Submitted merge job ${MERGE_JOB} (afterok:${SCAN_JOB})"
echo "OUT_DIR=${OUT_DIR}"
echo "Monitor: squeue -j ${SCAN_JOB},${MERGE_JOB}"
