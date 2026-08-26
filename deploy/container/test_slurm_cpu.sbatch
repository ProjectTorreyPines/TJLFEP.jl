#!/bin/bash
#SBATCH --job-name=tjlfep-container-cpu
#SBATCH --account=m3739
#SBATCH --constraint=cpu
#SBATCH --qos=debug
#SBATCH --time=0:25:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --output=tjlfep_container_cpu_%j.out
#
# Smoke-test the migrated TJLFEP podman-hpc image inside a Slurm job on a
# Perlmutter CPU node (AMD Milan, znver3 — one of the pkgimage clone targets).
# Pull/build + migrate first (see README), then:
#
#   TJLFEP_ENVIRONMENT=v2.0.14 sbatch --export=ALL deploy/container/test_slurm_cpu.sbatch
#
# Optional: SQUASH_DIR=... to test a shared image store instead of the
# per-user one.
#
# Pass: Sys.CPU_NAME reports znver3, the single-radius DIII-D critical-gradient
# solve completes offline on CPU, ends with SLURM CPU SMOKE OK.

set -euo pipefail

version="${TJLFEP_ENVIRONMENT:?set TJLFEP_ENVIRONMENT (e.g. v2.0.14)}"

podman=(podman-hpc)
[[ -n "${SQUASH_DIR:-}" ]] && podman=(podman-hpc --squash-dir "$SQUASH_DIR")

echo "node: $(hostname)"
lscpu | grep 'Model name'

"${podman[@]}" run --rm --network none \
    -e JULIA_NUM_THREADS="$SLURM_CPUS_PER_TASK" \
    "localhost/tjlfep:$version" tjlfep -e '
        println("cpu: ", Sys.CPU_NAME, "  threads: ", Threads.nthreads())
        @time import TJLFEP
        import CUDA
        println("CUDA.functional() = ", CUDA.functional(), " (CPU node — expect false)")
        ex = joinpath(pkgdir(TJLFEP), "examples", "DIIID_202017C42_500ms_v3.1")
        TJLFEP.run_gacode_scan_task(
            joinpath(ex, "input.gacode"),
            joinpath(ex, "input_singleradius_nb6.TGLFEP"),
            1; out_dir=mktempdir(), use_gpu=false, printout=true)
        println("SLURM CPU SMOKE OK  TJLFEP ", pkgversion(TJLFEP))'
