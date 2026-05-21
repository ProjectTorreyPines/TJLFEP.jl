# One radius of SCAN_N=20 (gacode + TGLFEP only).
#
#   sbatch batch_run_gacode_scan20_gpu_5nodes.sh   # 5 nodes, srun -n 20 (recommended)
#   sbatch batch_run_gacode_scan20_gpu_array.sh  # 20 separate array tasks
#   sbatch --dependency=afterok:<JOBID> batch_merge_gacode_scan20.sh

ENV["TJLFEP_FILE_ONLY"] = "1"

using Pkg
Pkg.activate(normpath(@__DIR__, ".."))

# Load CUDA before TJLFEP so TJLFEPCUDAExt registers GPU eigensolver hooks.
if get(ENV, "USE_GPU", "") == "1"
    using CUDA
end

using TJLFEP
using TJLF
using LinearAlgebra

BLAS.set_num_threads(1)

const ROOT = normpath(@__DIR__, "..")
const CASE = get(ENV, "CASE_DIR", joinpath(ROOT, "src", "DIIIDfiles", "202017C42_500ms_v3.1"))
const GACODE = get(ENV, "GACODE_FILE", joinpath(CASE, "input.gacode"))
const TGLFEP = get(ENV, "TGLFEP_FILE", joinpath(ROOT, "build", "debug_nb6", "input_scan20.TGLFEP"))
const OUT_DIR = get(ENV, "OUT_DIR", joinpath(@__DIR__, "gacode_scan20_$(get(ENV, "SLURM_JOB_ID", "local"))_tasks"))

@assert isfile(GACODE) "missing $GACODE"
@assert isfile(TGLFEP) "missing $TGLFEP"

opts, _, _ = preprocess_gacode_inputs(GACODE, TGLFEP)
scan_n = opts.SCAN_N

task0 = slurm_array_task_id()
scan_index = task0 + 1
@assert 1 <= scan_index <= scan_n "array task $task0 -> scan_index=$scan_index invalid for SCAN_N=$scan_n"

use_gpu = get(ENV, "USE_GPU", "") == "1" || TJLF.pick_device(:auto) === :gpu
printout = get(ENV, "TJLFEP_PRINTOUT", "0") == "1"

println("=== gacode scan task ===")
println("task_id=$task0 (array=$(get(ENV, "SLURM_ARRAY_TASK_ID", "—")) procid=$(get(ENV, "SLURM_PROCID", "—"))) scan_index=$scan_index / $scan_n")
println("OUT_DIR=$OUT_DIR")
println("device: ", use_gpu ? "GPU" : "CPU")
println("CUDA: functional=", TJLF._cuda_functional(), " solve=", TJLF._CUDA_SOLVE[] !== nothing)

t0 = time()
result = run_gacode_scan_task(
    GACODE, TGLFEP, scan_index;
    out_dir=OUT_DIR,
    use_gpu=use_gpu,
    printout=printout,
)
println("OK scan_index=$(result.scan_index) ir=$(result.ir) sfmin=$(result.sfmin) in $(round(time() - t0; digits=1)) s")
