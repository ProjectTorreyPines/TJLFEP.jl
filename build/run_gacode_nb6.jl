# Run TJLFEP from input.gacode only (no dump.gacode, MTGLF, or EXPRO).
#
#   TJLFEP_FILE_ONLY=1 julia --startup-file=no --project=.. build/run_gacode_nb6.jl

ENV["TJLFEP_FILE_ONLY"] = "1"

using Pkg
Pkg.activate(normpath(@__DIR__, ".."))

using TJLFEP

const ROOT = normpath(@__DIR__, "..")
const CASE = joinpath(ROOT, "src", "DIIIDfiles", "202017C42_500ms_v3.1")
const GACODE = joinpath(CASE, "input.gacode")
const TGLFEP = joinpath(ROOT, "build", "debug_nb6", "input.TGLFEP")

@assert isfile(GACODE)
@assert isfile(TGLFEP)

println("=== preprocess_gacode_inputs ===")
opts, prof, expro = preprocess_gacode_inputs(GACODE, TGLFEP)
println("NR=$(prof.NR) NS=$(prof.NS) SCAN_N=$(opts.SCAN_N) IR_EXP=$(opts.IR_EXP) IS_EP=$(opts.IS_EP)")

println("\n=== runTHD_from_gacode (SCAN_N from input.TGLFEP) ===")
width, kymark, SFmin, dpdr, dndr = runTHD_from_gacode(GACODE, TGLFEP; printout=false, parallel=:threads)
println("SFmin = ", SFmin)
println("done")
