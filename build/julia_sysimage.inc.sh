# Source from batch scripts: sets SYSIMAGE and JULIA_SYSIMAGE_ARGS.
TJLFEP_ROOT="${TJLFEP_ROOT:-/pscratch/sd/t/tneiser/.julia/dev/TJLFEP}"
SYSIMAGE="${TJLFEP_ROOT}/build/noTJLF_TJLFEP_sysimage.so"
export TJLFEP_SYSIMAGE="${SYSIMAGE}"
if [[ -f "${SYSIMAGE}" ]]; then
    JULIA_SYSIMAGE_ARGS=(--sysimage="${SYSIMAGE}")
    echo "Using sysimage: ${SYSIMAGE} ($(du -h "${SYSIMAGE}" | cut -f1))"
else
    JULIA_SYSIMAGE_ARGS=()
    echo "WARNING: sysimage not found at ${SYSIMAGE}; running without --sysimage"
fi
