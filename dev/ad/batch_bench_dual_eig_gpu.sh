#!/bin/bash -l
#SBATCH -A m3739_g
#SBATCH -q debug
#SBATCH -N 1
#SBATCH -t 00:20:00
#SBATCH -C gpu
#SBATCH -J dual_eig_bench
#SBATCH -o ad/bench_dual_eig_gpu_%j.out
#SBATCH -e ad/bench_dual_eig_gpu_%j.err
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --gpus-per-node=1

set -uo pipefail
module load cudatoolkit/12.9
module load julia/1.11.7
export JULIA_DEPOT_PATH="${PSCRATCH}/.julia${JULIA_DEPOT_PATH:+:${JULIA_DEPOT_PATH}}"
export TJLFEP_FILE_ONLY=1 USE_GPU=1 JULIA_CUDA_USE_COMPAT=false

TJLFEP_ROOT="${TJLFEP_ROOT:-/pscratch/sd/t/tneiser/.julia/dev/TJLFEP}"
GPU_SYSIMG="${TJLFEP_GPU_SYSIMAGE:-${TJLFEP_ROOT}/dev/TJLFEP_gpu_sysimage.so}"
if [[ -f "${GPU_SYSIMG}" ]]; then SYS=(--sysimage="${GPU_SYSIMG}"); echo "sysimage: ${GPU_SYSIMG}"; else SYS=(); echo "sysimage: JIT"; fi

cd "${TJLFEP_ROOT}/dev"
stdbuf -oL -eL julia --startup-file=no --project="${TJLFEP_ROOT}" "${SYS[@]}" -t 8 ad/bench_dual_eig_gpu.jl
echo "=== batch done ==="
