#!/bin/bash
# Orchestrate the -gpu tag: pull the CI-built image, bake the GPU-traced
# sysimage on an A100, layer it in, GPU-accept, and (optionally) publish.
# Run on a Perlmutter LOGIN node:
#
#   TJLFEP_ENVIRONMENT=v2.0.14 ./deploy/container/gpu-sysimage/bake_and_publish.sh
#
# Steps (each idempotent; rerun after a failure and completed steps are fast):
#   1. podman-hpc pull ghcr.io/projecttorreypines/tjlfep:<v> + migrate
#      (skipped when localhost/tjlfep:<v> already exists, e.g. a local build)
#   2. sbatch --wait batch_bake_gpu_sysimage.sh   (A100, m5377_g premium, ~1-4 h)
#   3. podman-hpc build Containerfile.gpu -> localhost/tjlfep:<v>-gpu + migrate
#   4. CPU acceptance on the -gpu image (asserts the sysimage is loaded)
#   5. submit test_slurm_gpu.sbatch for GPU acceptance
#   6. PUBLISH=1 only: push :<v>-gpu to ghcr (publish_ghcr.sh TAG_SUFFIX=-gpu)
#      — run this only after step 5's job passed.

set -euo pipefail

version="${TJLFEP_ENVIRONMENT:?set TJLFEP_ENVIRONMENT (e.g. v2.0.14)}"
scriptdir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
container_dir="$(dirname "$scriptdir")"
out="${SYSIMAGE_OUT:-$PSCRATCH/tjlfep_sysimage_out}"
ghcr="ghcr.io/projecttorreypines/tjlfep"

command -v podman-hpc >/dev/null || { echo "ERROR: podman-hpc not found — run on Perlmutter" >&2; exit 1; }

echo "### [1/6] Ensure localhost/tjlfep:$version is present + migrated"
if ! podman-hpc image exists "localhost/tjlfep:$version"; then
    podman-hpc pull "$ghcr:$version"
    podman-hpc tag "$ghcr:$version" "localhost/tjlfep:$version"
fi
podman-hpc migrate "localhost/tjlfep:$version"

if [[ ! -f "$out/sys_tjlfep_gpu.so" || "${REBAKE:-0}" == "1" ]]; then
    echo "### [2/6] Baking GPU sysimage on an A100 (sbatch --wait, ~1-4 h)"
    TJLFEP_ENVIRONMENT="$version" SYSIMAGE_OUT="$out" GPU_SYSIMAGE_DIR="$scriptdir" \
        sbatch --wait --export=ALL "$scriptdir/batch_bake_gpu_sysimage.sh"
else
    echo "### [2/6] Reusing existing $out/sys_tjlfep_gpu.so (REBAKE=1 to force)"
fi
[[ -f "$out/sys_tjlfep_gpu.so" ]] || { echo "ERROR: bake produced no sysimage" >&2; exit 1; }

echo "### [3/6] Layering sysimage -> localhost/tjlfep:$version-gpu"
podman-hpc build -t "localhost/tjlfep:$version-gpu" \
    -f "$scriptdir/Containerfile.gpu" \
    --build-arg BASE="localhost/tjlfep:$version" \
    "$out"
podman-hpc migrate "localhost/tjlfep:$version-gpu"

echo "### [4/6] CPU acceptance on the -gpu image"
TJLFEP_ENVIRONMENT="$version-gpu" "$container_dir/acceptance.sh"

echo "### [5/6] Submitting GPU acceptance job (check its output before publishing)"
TJLFEP_ENVIRONMENT="$version-gpu" sbatch --export=ALL "$container_dir/test_slurm_gpu.sbatch"

if [[ "${PUBLISH:-0}" == "1" ]]; then
    echo "### [6/6] Publishing $ghcr:$version-gpu"
    TJLFEP_ENVIRONMENT="$version" TAG_SUFFIX=-gpu "$container_dir/publish_ghcr.sh"
else
    echo "### [6/6] Skipped publish. After the GPU acceptance job passes, run:"
    echo "    TJLFEP_ENVIRONMENT=$version TAG_SUFFIX=-gpu $container_dir/publish_ghcr.sh"
fi
