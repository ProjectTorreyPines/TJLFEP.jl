#!/bin/bash -l
set -euo pipefail
module load julia/1.11.7 2>/dev/null || true
export JULIA_DEPOT_PATH="${PSCRATCH}/.julia${JULIA_DEPOT_PATH:+:${JULIA_DEPOT_PATH}}"
cd "$(dirname "$0")"
julia --project=.. plot_scan20_timing.jl
