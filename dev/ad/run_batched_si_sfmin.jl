using Pkg; Pkg.activate(normpath(@__DIR__, "..", ".."))
import CUDA
using TJLFEP, TJLF, Printf
using TJLFEP: preprocess_gacode_inputs, kwscale_scan

# Run the wired inner=:batched_si (fixed-shift hybrid) solver over the full 20-radius scan at a
# given N_BASIS and write "idx IR sfmin" lines (same format as gacode_*_tasks/sfmin_scan.txt) so it
# can be overlaid on the stored nb32 grid/Fortran reference. No dense golden is recomputed here —
# the reference values already live in the repo (build/gacode_nb32_scan20_jgpu_*/sfmin_scan.txt).
CASE = get(ENV, "CASE_DIR", normpath(@__DIR__, "..", "..", "examples", "DIIID_202017C42_500ms_v3.1"))
NB   = parse(Int, get(ENV, "NB", "32"))
GAC  = joinpath(CASE, "input.gacode")
TGL  = joinpath(CASE, "input_scan20_nb$(NB).TGLFEP")
USE_GPU = get(ENV, "USE_GPU", "1") == "1"
OUT  = get(ENV, "OUT", normpath(@__DIR__, "batched_si_sfmin_nb$(NB).txt"))
grid = (nfactor = parse(Int, get(ENV,"NFACTOR","8")), nefwid = parse(Int, get(ENV,"NEFWID","8")),
        nkyhat  = parse(Int, get(ENV,"NKYHAT","4")), k_max  = parse(Int, get(ENV,"KMAX","4")))

base_ep, prof, _ = preprocess_gacode_inputs(GAC, TGL)
scan_n = Int(base_ep.SCAN_N)
# RADII_IDX = subset of scan indices (1..scan_n) this process handles (for multi-GPU sharding).
idxs = let r = get(ENV, "RADII_IDX", "")
    isempty(r) ? collect(1:scan_n) : parse.(Int, split(r, ","))
end
@printf("batched_si sfmin run: nb=%d grid=%s use_gpu=%s idxs=%s -> %s\n",
        NB, grid, USE_GPU, join(idxs, ","), basename(OUT))

twall = @elapsed open(OUT, "w") do io
    for i in idxs
        ep = deepcopy(base_ep); ep.IR = base_ep.IR_EXP[i]
        ep.WIDTH_IN_FLAG = false; ep.MODE_IN = 2; ep.KY_MODEL = 3; ep.PROCESS_IN = 5
        ep.FACTOR_IN = Float64(base_ep.FACTOR[i])
        t = @elapsed begin
            _g, epo, = kwscale_scan(ep, prof, false; use_gpu=USE_GPU, inner=:batched_si, grid...)
        end
        sf = Float64(epo.FACTOR_IN)
        @printf(io, "%d %d %.16g\n", i, base_ep.IR_EXP[i], sf); flush(io)
        @printf("  ir=%3d  sfmin=%9.4g  ky=%.3f w=%.3f  %5.0fs\n",
                base_ep.IR_EXP[i], sf, Float64(epo.KYMARK), Float64(epo.WIDTH_IN), t); flush(stdout)
    end
end
@printf("wrote %s  TOTAL_WALL_S=%.1f  (%d radii)\n", OUT, twall, length(idxs))
