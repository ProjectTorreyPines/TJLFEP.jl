# Merge the per-radius task outputs of a scan into the final SFmin profile + alpha files.
# Run with an explicit --project (the TJLFEP checkout, or the published env_lean/env_full
# that matches a prebuilt sysimage); the slurm/ templates call this automatically.
#
#   OUT_DIR=runs/<run>_tasks GACODE_FILE=... TGLFEP_FILE=... \
#     julia --project=. slurm/common/merge_gacode_scan20_array.jl

using TJLFEP

const ROOT = normpath(@__DIR__, "..", "..")
const CASE = get(ENV, "CASE_DIR", joinpath(ROOT, "examples", "DIIID_202017C42_500ms_v3.1"))
# NB: `GACODE_PATH`, not `GACODE` — the generic GPU sysimage bakes the `GACODE` package module
# into Main, so a top-level `const GACODE = ...` collides ("invalid redefinition of constant").
const GACODE_PATH = get(ENV, "GACODE_FILE", joinpath(CASE, "input.gacode"))
const TGLFEP = get(ENV, "TGLFEP_FILE", joinpath(CASE, "input_scan20_nb32.TGLFEP"))
const OUT_DIR = get(() -> error("set OUT_DIR to the <run>_tasks directory holding task_*.jls"), ENV, "OUT_DIR")

@assert isfile(GACODE_PATH)
@assert isfile(TGLFEP)
@assert isdir(OUT_DIR)

println("=== finalize_gacode_scan ===")
println("OUT_DIR=$OUT_DIR")

t0 = time()
width, kymark, SFmin, dpdr, dndr = finalize_gacode_scan(GACODE_PATH, TGLFEP, OUT_DIR; printout=true)
println("SFmin = ", SFmin)
println("done in $(round(time() - t0; digits=1)) s")
