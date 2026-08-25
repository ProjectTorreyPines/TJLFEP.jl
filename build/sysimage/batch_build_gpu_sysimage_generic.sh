#!/bin/bash -l
# Build the GENERIC GPU sysimage (CUDA + TJLF + TJLFEP full + FUSE/IMAS stack) on a GPU node.
# Works for both the file-based scan and the IMAS/FUSE actor path on GPU.
# Output: build/TJLFEP_gpu_generic_sysimage.so
#
#   cd build && sbatch sysimage/batch_build_gpu_sysimage_generic.sh
#
#SBATCH -A m3739_g
#SBATCH -q regular
#SBATCH -N 1
#SBATCH -t 02:00:00
#SBATCH -C gpu
#SBATCH -J TJLFEP_gpu_gen_sysimg
#SBATCH -o build_gpu_generic_sysimage_%j.out
#SBATCH -e build_gpu_generic_sysimage_%j.err
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --gpus-per-node=1

set -euo pipefail
# New files/dirs (CFS publishes, depot writes) must stay m3739 group-readable.
umask 007

module load cudatoolkit/12.9
module load julia/1.11.7
# TJLFEP_DEPOT: base depot for the bake. For group-shareable images this MUST be the CFS
# depot (paths are baked into the image); defaults to the private scratch depot for
# personal builds.
# TJLFEP_BAKE_DEPOT: optional isolated primary depot for the bake (stacked in front of
# the shared depot, which stays read-only for packages/artifacts). Avoids cross-flavor
# precompile-cache collisions (pkgimages=no emit vs pkgimages=true shared caches, or two
# projects resolving different versions of the same package into one compiled/ dir).
DEPOT_BASE="${TJLFEP_DEPOT:-${PSCRATCH}/.julia}"
if [[ -n "${TJLFEP_BAKE_DEPOT:-}" ]]; then
    export JULIA_DEPOT_PATH="${TJLFEP_BAKE_DEPOT}:${DEPOT_BASE}"
else
    export JULIA_DEPOT_PATH="${DEPOT_BASE}"
fi
mkdir -p "${JULIA_DEPOT_PATH%%:*}/compiled"

# GENERIC image: IMAS/GACODE/TurbulentTransport are weak deps of TJLFEP and load the
# TJLFEPIMASExt extension (the dd/FUSE actor path) automatically when present. The build
# project (TJLFEP_ROOT) precompiles the ext because the weak deps are in its manifest, and
# the precompile_execution file `using`s FUSE -> the ext + FUSE are baked into the image.
# JULIA_CUDA_USE_COMPAT=false matches the runtime env.
export JULIA_CUDA_USE_COMPAT=false

TJLFEP_ROOT="${TJLFEP_ROOT:-/pscratch/sd/t/tneiser/.julia/dev/TJLFEP}"
# The generic image bakes the FUSE/IMAS stack, which can only be resolved from the FUSE
# project (FUSE is not -- and cannot be -- a TJLFEP dep). build_gpu_sysimage_generic.jl
# activates FUSE_ROOT and stacks TJLFEP_ROOT for PackageCompiler.
export FUSE_ROOT="${FUSE_ROOT:-${TJLFEP_ROOT}/../FUSE}"
cd "${TJLFEP_ROOT}/build"

echo "=== TJLFEP GENERIC GPU sysimage build ==="
echo "host: $(hostname)  date: $(date)"
julia --version
nvidia-smi -L 2>/dev/null | head -1 || true
echo "TJLF: $(git -C "${TJLFEP_ROOT}/../TJLF" rev-parse --abbrev-ref HEAD 2>/dev/null) $(git -C "${TJLFEP_ROOT}/../TJLF" log -1 --oneline 2>/dev/null)"
echo "TJLFEP: $(git -C "${TJLFEP_ROOT}" rev-parse --abbrev-ref HEAD 2>/dev/null) $(git -C "${TJLFEP_ROOT}" log -1 --oneline 2>/dev/null)"

# Instantiate + precompile dependencies (FUSE/IMAS/CUDA/TurbulentTransport/...). The weak
# deps resolve from the FUSE project (not TJLFEP's own manifest); loading them there builds
# TJLFEPIMASExt. No forced recompile needed with the extension model.
julia --project="${FUSE_ROOT}" -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
julia --project="${TJLFEP_ROOT}" -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

julia --project="${TJLFEP_ROOT}" -t "${SLURM_CPUS_PER_TASK:-32}" sysimage/build_gpu_sysimage_generic.jl

SO="${TJLFEP_ROOT}/build/TJLFEP_gpu_generic_sysimage.so"
if [[ ! -f "${SO}" ]]; then
    echo "ERROR: sysimage not found at ${SO}"
    exit 1
fi
ls -lh "${SO}"

# Self-check: load the freshly built image and assert it is GENERIC -- the TJLFEPIMASExt
# extension is loaded (dd/FUSE actor path), FUSE is baked, and runTHD has the ::IMAS.dd
# method. Fail loudly otherwise so we never ship a silently file-only "generic" image.
echo "=== verifying generic image ==="
julia --startup-file=no --sysimage="${SO}" --project="${TJLFEP_ROOT}" -e '
    fuse = Base.PkgId(Base.UUID("e64856f0-3bb8-4376-b4b7-c03396503992"), "FUSE")
    ext = Base.get_extension(TJLFEP, :TJLFEPIMASExt)
    @assert ext !== nothing  "GENERIC build FAILED: TJLFEPIMASExt not loaded (IMAS/GACODE/TurbulentTransport not baked)"
    @assert haskey(Base.loaded_modules, fuse)  "GENERIC build FAILED: FUSE not baked into image"
    @assert length(methods(TJLFEP.runTHD)) >= 2  "GENERIC build FAILED: runTHD(::IMAS.dd) method missing (actor path not compiled in)"
    println("generic image OK: TJLFEPIMASExt loaded, FUSE baked, runTHD methods=", length(methods(TJLFEP.runTHD)))'

# Publish to the shared CFS location with a .sha sidecar recording per-repo
# `git describe --always --dirty` provenance + the julia version the image requires.
CFS_DIR="${TJLFEP_CFS_DIR:-/global/cfs/cdirs/m3739/TJLFEP}"
if mkdir -p "${CFS_DIR}" 2>/dev/null; then
    cp -f "${SO}" "${CFS_DIR}/"
    chgrp m3739 "${CFS_DIR}/$(basename "${SO}")" 2>/dev/null || true
    SHA_FILE="${CFS_DIR}/$(basename "${SO}").sha"
    {
        for repo in "${TJLFEP_ROOT}" "${TJLFEP_ROOT}/../TJLF" "${FUSE_ROOT}" \
                    "${FUSE_ROOT}/../IMASdd" "${FUSE_ROOT}/../IMASggd"; do
            printf '%s=%s\n' "$(basename "${repo}")" \
                "$(git -C "${repo}" describe --always --dirty --tags 2>/dev/null || echo unknown)"
        done
        printf 'julia=%s\n' "$(julia --version | awk '{print $NF}')"
    } > "${SHA_FILE}"
    chgrp m3739 "${SHA_FILE}" 2>/dev/null || true
    echo "=== published to ${CFS_DIR}/$(basename "${SO}") ==="
    cat "${SHA_FILE}"
else
    echo "WARNING: could not write ${CFS_DIR}; sysimage left at ${SO} only"
fi

echo "=== GENERIC GPU sysimage build OK: ${SO} ==="
