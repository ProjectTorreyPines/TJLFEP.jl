# Regenerate the production AD sfmin(radius) profile for DIII-D N_BASIS=32, SCAN_N=20 using
# the exact production code path (run_gacode_scan_task with inner=:threads, use_gpu=true).
#   SOLVER         (default robust_ad) — :ad (fast descent) or :robust_ad (robust,
#                  global-min over the (ky,w) grid; tracks the Fortran/grid sfmin).
#   REFINE_ROUNDS  (default 1) — accuracy/speed knob for :robust_ad (rounds of (ky,w)
#                  window narrowing; 0=coarse grid min, higher=better off-node resolution).
# Writes "ad_threads_sfmin_nb32_<solver>[_r<refine>].txt" (cols: scan_index ir sfmin) so the
# fast and robust profiles can be compared side by side. Run JIT on a GPU node (no sysimage).

using CUDA
using TJLF
using TJLFEP
using Printf

const CASE    = joinpath(@__DIR__, "..", "..", "examples", "DIIID_202017C42_500ms_v3.1")
const GACODE  = joinpath(CASE, "input.gacode")
const TGLFEP  = joinpath(CASE, "input_scan20_nb32.TGLFEP")
const SOLVER  = Symbol(get(ENV, "SOLVER", "robust_ad"))
const REFINE  = parse(Int, get(ENV, "REFINE_ROUNDS", "1"))
const EXTMODE = get(ENV, "AD_EXTEND_MODE", "locate")   # :ad extend strategy (locate|wide|only)
# Tag the output by solver, and for :ad also by extend mode (so :only/:wide do not
# clobber the default :ad/:locate profile). Legacy :ad/:locate keeps the bare "ad" tag.
const TAG     = SOLVER === :robust_ad ? "robust_ad_r$(REFINE)" :
                SOLVER === :ad ? (EXTMODE == "locate" ? "ad" : "ad_$(EXTMODE)") :
                String(SOLVER)
const OUT_DIR = joinpath(@__DIR__, "ad_threads_sfmin_nb32_$(TAG)_tasks")
const OUT_TXT = joinpath(@__DIR__, "ad_threads_sfmin_nb32_$(TAG).txt")

function main()
    @assert CUDA.functional() "run on a GPU node"
    opts, _, _ = preprocess_gacode_inputs(GACODE, TGLFEP)
    scan_n = opts.SCAN_N
    # SI_LIST (e.g. "8,9,10") re-runs only those scan indices and UPDATES the existing OUT_TXT in
    # place (keeping all other radii); empty = full 1:scan_n sweep (overwrite, as before).
    si_env = strip(get(ENV, "SI_LIST", ""))
    si_list = isempty(si_env) ? collect(1:scan_n) : parse.(Int, split(si_env, ','))
    @printf("DIII-D N_BASIS=%d  SCAN_N=%d  solver=%s refine_rounds=%d  GPU=%s  si_list=%s\n",
            opts.N_BASIS, scan_n, SOLVER, REFINE, CUDA.name(first(CUDA.devices())), string(si_list))
    flush(stdout)

    # seed the table from disk so a subset re-run preserves untouched radii
    rows = Dict{Int,Tuple{Int,Float64}}()
    if isfile(OUT_TXT)
        for line in eachline(OUT_TXT)
            p = split(strip(line)); length(p) >= 3 || continue
            rows[parse(Int, p[1])] = (parse(Int, p[2]), parse(Float64, p[3]))
        end
    end

    for si in si_list
        t = @elapsed r = run_gacode_scan_task(GACODE, TGLFEP, si;
                out_dir = OUT_DIR, use_gpu = true, inner = :threads,
                solver = SOLVER, refine_rounds = REFINE)
        rows[si] = (r.ir, r.sfmin)
        @printf("done si=%2d ir=%3d sfmin=%.6g  in %.1fs\n", si, r.ir, r.sfmin, t)
        flush(stdout)
        open(OUT_TXT, "w") do io
            for k in sort(collect(keys(rows)))
                ir, sf = rows[k]
                @printf(io, "%d %d %.10g\n", k, ir, sf)
            end
        end
    end
    println("=== wrote ", OUT_TXT, " ===")
end

main()
