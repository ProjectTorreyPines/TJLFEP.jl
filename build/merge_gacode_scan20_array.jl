# Merge Slurm array task outputs and write alpha profiles.
#
#   OUT_DIR=.../gacode_scan20_<jobid>_tasks julia --project=.. merge_gacode_scan20_array.jl

ENV["TJLFEP_FILE_ONLY"] = "1"

using Pkg
Pkg.activate(normpath(@__DIR__, ".."))

using TJLFEP

const ROOT = normpath(@__DIR__, "..")
const CASE = get(ENV, "CASE_DIR", joinpath(ROOT, "src", "DIIIDfiles", "202017C42_500ms_v3.1"))
const GACODE = get(ENV, "GACODE_FILE", joinpath(CASE, "input.gacode"))
const TGLFEP = get(ENV, "TGLFEP_FILE", joinpath(ROOT, "build", "debug_nb6", "input_scan20.TGLFEP"))
const OUT_DIR = get(() -> error("set OUT_DIR to gacode_scan20_<jobid>_tasks directory"), ENV, "OUT_DIR")

@assert isfile(GACODE)
@assert isfile(TGLFEP)
@assert isdir(OUT_DIR)

println("=== finalize_gacode_scan ===")
println("OUT_DIR=$OUT_DIR")

t0 = time()
width, kymark, SFmin, dpdr, dndr = finalize_gacode_scan(GACODE, TGLFEP, OUT_DIR; printout=true)
println("SFmin = ", SFmin)
println("done in $(round(time() - t0; digits=1)) s")
