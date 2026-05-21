# Run TJLFEP from input.gacode only — N_BASIS=32, SCAN_N=1.
#
#   USE_GPU=1 julia --startup-file=no --project=.. build/run_gacode_nb32.jl

ENV["TJLFEP_FILE_ONLY"] = "1"

using Pkg
Pkg.activate(normpath(@__DIR__, ".."))

if get(ENV, "USE_GPU", "") == "1"
    using CUDA
end

using TJLFEP
using TJLF

const ROOT = normpath(@__DIR__, "..")
const CASE = get(ENV, "CASE_DIR", joinpath(ROOT, "src", "DIIIDfiles", "202017C42_500ms_v3.1"))
const GACODE = get(ENV, "GACODE_FILE", joinpath(CASE, "input.gacode"))
const TGLFEP = get(ENV, "TGLFEP_FILE", joinpath(ROOT, "build", "debug_nb32", "input.TGLFEP"))

@assert isfile(GACODE)
@assert isfile(TGLFEP)

println("=== preprocess_gacode_inputs (N_BASIS=32) ===")
opts, prof, expro = preprocess_gacode_inputs(GACODE, TGLFEP)
println("NR=$(prof.NR) NS=$(prof.NS) SCAN_N=$(opts.SCAN_N) N_BASIS=$(opts.N_BASIS) IR_EXP=$(opts.IR_EXP)")

use_gpu = get(ENV, "USE_GPU", "") == "1" || TJLF.pick_device(:auto) === :gpu
println("\n=== runTHD_from_gacode ===")
println("device: ", use_gpu ? "GPU" : "CPU")
width, kymark, SFmin, dpdr, dndr = runTHD_from_gacode(GACODE, TGLFEP; printout=false, use_gpu=use_gpu, parallel=:threads)
println("SFmin = ", SFmin)
println("done")
