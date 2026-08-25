#!/bin/bash -l
# Build the FILE-ONLY GPU sysimage (CUDA + TJLF + TJLFEP standalone, no FUSE/IMAS) on a GPU node.
# This is what a TGLF-EP user running the file-based scan path gets before going FUSE-native:
# leaner + faster-loading than the generic image because the IMAS/FUSE stack is not baked.
# Output: build/TJLFEP_gpu_sysimage.so
#
#   cd build && sbatch sysimage/batch_build_gpu_sysimage_fileonly.sh
#
#SBATCH -A m3739_g
#SBATCH -q premium
#SBATCH -N 1
#SBATCH -t 02:00:00
#SBATCH -C gpu
#SBATCH -J TJLFEP_gpu_fo_sysimg
#SBATCH -o build_gpu_fileonly_sysimage_%j.out
#SBATCH -e build_gpu_fileonly_sysimage_%j.err
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --gpus-per-node=1

set -euo pipefail
# New files/dirs (CFS publishes, depot writes) must stay m3739 group-readable.
umask 007

module load cudatoolkit/12.9
module load julia/1.11.7
# TJLFEP_DEPOT: base depot for the bake. For group-shareable images this MUST be the CFS
# depot (paths are baked into the image); defaults to the private scratch depot.
export JULIA_DEPOT_PATH="${TJLFEP_DEPOT:-${PSCRATCH}/.julia}"
mkdir -p "${JULIA_DEPOT_PATH}/compiled"

# JULIA_CUDA_USE_COMPAT=false matches the runtime worker env.
export JULIA_CUDA_USE_COMPAT=false

TJLFEP_ROOT="${TJLFEP_ROOT:-/pscratch/sd/t/tneiser/.julia/dev/TJLFEP}"
cd "${TJLFEP_ROOT}/build"

echo "=== TJLFEP FILE-ONLY GPU sysimage build ==="
echo "host: $(hostname)  date: $(date)"
julia --version
nvidia-smi -L 2>/dev/null | head -1 || true
echo "TJLF: $(git -C "${TJLFEP_ROOT}/../TJLF" rev-parse --abbrev-ref HEAD 2>/dev/null) $(git -C "${TJLFEP_ROOT}/../TJLF" log -1 --oneline 2>/dev/null)"
echo "TJLFEP: $(git -C "${TJLFEP_ROOT}" rev-parse --abbrev-ref HEAD 2>/dev/null) $(git -C "${TJLFEP_ROOT}" log -1 --oneline 2>/dev/null)"

# Instantiate + precompile the TJLFEP project deps (CUDA/TJLF/...). The extension only builds/
# loads when IMAS/GACODE/TurbulentTransport are loaded, which the file-only workload never does.
julia --project="${TJLFEP_ROOT}" -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

julia --project="${TJLFEP_ROOT}" -t "${SLURM_CPUS_PER_TASK:-32}" sysimage/build_gpu_sysimage_fileonly.jl

SO="${TJLFEP_ROOT}/build/TJLFEP_gpu_sysimage.so"
if [[ ! -f "${SO}" ]]; then
    echo "ERROR: sysimage not found at ${SO}"
    exit 1
fi
ls -lh "${SO}"

# Self-check: load the freshly built image and assert it is FILE-ONLY -- the TJLFEPIMASExt
# extension is NOT loaded and FUSE is NOT baked. Fail loudly otherwise so we never ship a
# silently generic image under the file-only name.
echo "=== verifying file-only image ==="
julia --startup-file=no --sysimage="${SO}" --project="${TJLFEP_ROOT}" -e '
    fuse = Base.PkgId(Base.UUID("e64856f0-3bb8-4376-b4b7-c03396503992"), "FUSE")
    ext = Base.get_extension(TJLFEP, :TJLFEPIMASExt)
    @assert ext === nothing  "FILE-ONLY build FAILED: TJLFEPIMASExt is loaded (IMAS/GACODE/TurbulentTransport got baked)"
    @assert !haskey(Base.loaded_modules, fuse)  "FILE-ONLY build FAILED: FUSE is baked into image"
    println("file-only image OK: TJLFEPIMASExt dormant, FUSE not baked, runTHD methods=", length(methods(TJLFEP.runTHD)))'

# Publish to the shared CFS location with a .sha provenance sidecar (TJLF/TJLFEP only:
# the file-only image bakes no FUSE/IMAS stack).
CFS_DIR="${TJLFEP_CFS_DIR:-/global/cfs/cdirs/m3739/TJLFEP}"
if mkdir -p "${CFS_DIR}" 2>/dev/null; then
    cp -f "${SO}" "${CFS_DIR}/"
    chgrp m3739 "${CFS_DIR}/$(basename "${SO}")" 2>/dev/null || true
    SHA_FILE="${CFS_DIR}/$(basename "${SO}").sha"
    {
        for repo in "${TJLFEP_ROOT}" "${TJLFEP_ROOT}/../TJLF"; do
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

echo "=== FILE-ONLY GPU sysimage build OK: ${SO} ==="
