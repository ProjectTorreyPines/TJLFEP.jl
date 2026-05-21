# Run TJLFEP from input.gacode only (no dump.gacode, MTGLF, or EXPRO).
#
#   TJLFEP_FILE_ONLY=1 julia --startup-file=no --project=.. build/run_gacode_nb6.jl

ENV["TJLFEP_FILE_ONLY"] = "1"

using Pkg
Pkg.activate(normpath(@__DIR__, ".."))

if get(ENV, "USE_GPU", "") == "1"
    using CUDA
end

using TJLFEP
using TJLF

const ROOT = normpath(@__DIR__, "..")
const CASE = joinpath(ROOT, "src", "DIIIDfiles", "202017C42_500ms_v3.1")
const GACODE = joinpath(CASE, "input.gacode")
const TGLFEP = joinpath(ROOT, "build", "debug_nb6", "input.TGLFEP")

@assert isfile(GACODE)
@assert isfile(TGLFEP)

println("=== preprocess_gacode_inputs ===")
opts, prof, expro = preprocess_gacode_inputs(GACODE, TGLFEP)
println("NR=$(prof.NR) NS=$(prof.NS) SCAN_N=$(opts.SCAN_N) IR_EXP=$(opts.IR_EXP) IS_EP=$(opts.IS_EP)")

use_gpu = get(ENV, "USE_GPU", "") == "1" || TJLF.pick_device(:auto) === :gpu
println("\n=== runTHD_from_gacode (SCAN_N from input.TGLFEP) ===")
println("device: ", use_gpu ? "GPU" : "CPU", "  pick_device(:auto)=", TJLF.pick_device(:auto))
println("CUDA functional: ", TJLF._cuda_functional(), "  _CUDA_SOLVE set: ", TJLF._CUDA_SOLVE[] !== nothing)
width, kymark, SFmin, dpdr, dndr = runTHD_from_gacode(GACODE, TGLFEP; printout=false, use_gpu=use_gpu, parallel=:threads)
println("SFmin = ", SFmin)
println("done")
