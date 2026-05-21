#!/bin/bash
# Submit sysimage build, then GPU smoke test when build succeeds.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

BUILD_OUT=$(sbatch --parsable batch_build_sysimage.sh)
echo "Submitted sysimage build: job ${BUILD_OUT}"

# Perlmutter: use --gpus-per-node=1 (not --gpus-per-task without -n)
TEST_OUT=$(sbatch --parsable --dependency=afterok:"${BUILD_OUT}" batch_smoke_test.sh)
echo "Submitted smoke test (after build): job ${TEST_OUT}"

echo ""
echo "Monitor:"
echo "  squeue -u \$USER"
echo "  tail -f ${SCRIPT_DIR}/build_sysimage_${BUILD_OUT}.out"
echo "  tail -f ${SCRIPT_DIR}/smoke_test_${TEST_OUT}.out"
