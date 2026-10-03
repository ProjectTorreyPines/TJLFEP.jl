#!/bin/bash
# Generalized per-node MPS launch + GPU-pinning wrapper for the production scan, parameterized
# by GPUS_PER_RADIUS (1, 2, or 4) so the same code drives the 20N / 10N / 5N layouts:
#   GPUS_PER_RADIUS=4  -> 1 task/node  (20 nodes): 1 radius/node on all 4 A100s
#   GPUS_PER_RADIUS=2  -> 2 tasks/node (10 nodes): 2 radii/node, 2 A100s each
#   GPUS_PER_RADIUS=1  -> 4 tasks/node ( 5 nodes): 4 radii/node, 1 A100 each
# All are a single wave: total tasks = total radii (20), SCAN_INDEX = global procid + 1.
#
# LAUNCH-ORDER FIX (critical): the MPS control daemon MUST be up before ANY client calls cuInit
# or clients fall back to private contexts and time-slice. Node-local rank 0 starts the daemon;
# ALL ranks BLOCK on the control pipe before exec'ing the (CUDA-initializing) command.
#
#   GPUS_PER_RADIUS=N srun -n 20 --ntasks-per-node=K [--gpus-per-node=4] mps-scan-wrapper.sh <cmd>...

set -uo pipefail

export CUDA_MPS_PIPE_DIRECTORY="${CUDA_MPS_PIPE_DIRECTORY:-/tmp/nvidia-mps}"
export CUDA_MPS_LOG_DIRECTORY="${CUDA_MPS_LOG_DIRECTORY:-/tmp/nvidia-log}"
CONTROL="${CUDA_MPS_PIPE_DIRECTORY}/control"

LOCALID="${SLURM_LOCALID:-0}"
G="${GPUS_PER_RADIUS:-2}"

# Node-local rank 0 starts the daemon with ALL node GPUs visible (one server/GPU on first connect).
if [ "${LOCALID}" -eq 0 ]; then
    if [ ! -e "${CONTROL}" ]; then
        mkdir -p "${CUDA_MPS_PIPE_DIRECTORY}" "${CUDA_MPS_LOG_DIRECTORY}"
        CUDA_VISIBLE_DEVICES="${SLURM_JOB_GPUS:-${SLURM_STEP_GPUS:-0,1,2,3}}" \
            nvidia-cuda-mps-control -d
    fi
fi

# BARRIER: every rank waits until the control pipe exists (daemon up) before any CUDA init.
for _ in $(seq 1 120); do
    [ -e "${CONTROL}" ] && break
    sleep 0.5
done
if [ ! -e "${CONTROL}" ]; then
    echo "[mps-scan-wrapper localid=${LOCALID}] ERROR: MPS control pipe ${CONTROL} never appeared" >&2
    exit 1
fi
sleep 3   # small extra margin so the daemon is fully ready to accept clients

# This task's GPUs: contiguous block of G devices starting at localid*G (node-local indices).
START=$(( LOCALID * G ))
DEVS=""
for ((i=0; i<G; i++)); do
    DEVS="${DEVS:+${DEVS},}$(( START + i ))"
done
export CUDA_VISIBLE_DEVICES="${DEVS}"
export TEAM_GPUS="${DEVS}"
# One radius per global task.
export SCAN_INDEX="$(( ${SLURM_PROCID:-0} + 1 ))"

echo "[mps-scan-wrapper] node=$(hostname) localid=${LOCALID} procid=${SLURM_PROCID:-0} G=${G} -> SCAN_INDEX=${SCAN_INDEX} GPUs=${DEVS}"
exec "$@"
