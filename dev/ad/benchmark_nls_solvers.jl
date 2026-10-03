#!/usr/bin/env julia
# Benchmark the derivative-free nonlinear-optimizer solvers (:multistart, :nlopt) against
# the reference tiers (:grid = Fortran-equivalent w>=1 wide box, :ad :locate = narrow-width AD) on
# the 20-radius DIII-D 202017C42_500ms_v3.1 scan. For each solver x radius it records the critical
# factor sfmin, the marking (kyhat, width), the eigensolve/confirm eval count, and wall time, then
# prints a comparison table (ratios vs :ad :locate and :grid), writes a CSV, and an sfmin-vs-IR plot.
#
# The point is a HEAD-TO-HEAD of the derivative-free search (multistart+BOBYQA, DIRECT-L+BOBYQA)
# vs the AD path: does it track :ad :locate at the core, capture the
# narrow-width edge modes the :grid box misses (IR >~ 65), and at what cost / robustness.
#
# Usage (module load julia/1.11.7; JULIA_DEPOT_PATH set per the workspace rule):
#   julia --project=. dev/ad/benchmark_nls_solvers.jl
# Env knobs:
#   NB=32            N_BASIS (input_scan20_nb{NB}.TGLFEP; default 32)
#   RADII=1,8,20     1-based scan indices to run (default: all SCAN_N)
#   SOLVERS=grid,ad,multistart,nlopt   which solvers (default grid,ad,multistart,nlopt)
#   USE_GPU=1        use the GPU eigensolve path (imports+inits CUDA; default off = CPU threads)
#   NLS_LOCAL_ALGO=LN_BOBYQA local descent algo for :multistart (LN_SBPLX|LN_NELDERMEAD)
#   NLOPT_ALGO=GN_DIRECT_L   global NLopt algorithm for :nlopt
#   NLOPT_MAXEVAL=40         :nlopt global budget
using Pkg
Pkg.activate(normpath(@__DIR__, "..", ".."))
using TJLFEP
using Printf
using Plots
using Plots.PlotMeasures: mm

# GPU eigensolve path. TJLF's GPU solve lives in a package extension that only activates once
# `CUDA` is imported into the session; a bare `USE_GPU=1` without the import leaves it dormant and
# `_gpu_solve!` throws "CUDA extension not loaded". Import CUDA at TOP LEVEL here (guarded by the
# env flag) so the extension's methods are visible in the running world — doing the import inside a
# function and immediately calling the freshly-loaded methods hits a world-age error instead.
const _GPU_ENV = get(ENV, "USE_GPU", "")
const USE_GPU  = isempty(_GPU_ENV) ? false : (_GPU_ENV != "0")
if USE_GPU
    import CUDA
    CUDA.functional() || error("USE_GPU=$(_GPU_ENV) but CUDA.functional() is false on this host")
    CUDA.device!(first(CUDA.devices()))
end

const CASE = normpath(@__DIR__, "..", "..", "examples", "DIIID_202017C42_500ms_v3.1")
const NB   = parse(Int, get(ENV, "NB", "32"))
const GAC  = joinpath(CASE, "input.gacode")
const TGL  = joinpath(CASE, "input_scan20_nb$(NB).TGLFEP")

const INNER   = Symbol(get(ENV, "INNER", "threads"))   # threads | batched_si | mps_team

# Prepare one radius's Options exactly as mainsub's PROCESS_IN=5 path does, then hand it to a
# low-level solver. Returns a fresh deepcopy per call so solvers don't cross-contaminate.
function _prep(base_ep, i)
    ep = deepcopy(base_ep)
    ep.IR = base_ep.IR_EXP[i]
    ep.WIDTH_IN_FLAG = false
    ep.MODE_IN = 2
    ep.KY_MODEL = 3
    ep.PROCESS_IN = 5
    ep.FACTOR_IN = Float64(base_ep.FACTOR[i])   # per-radius scan ceiling
    ep
end

# Each runner returns (; sfmin, ky, w, evals, wall). `evals` is a combined eigensolve-equivalent
# count (grid: all IFLUX=true; the others: mostly cheap IFLUX=false + a few IFLUX=true confirms).
function run_grid(ep, prof)
    t = @elapsed begin
        _g, epo, _p, _mq, _sb, _wf = TJLFEP.kwscale_scan(ep, prof, false; use_gpu=USE_GPU, inner=INNER)
    end
    # kwscale_scan defaults: nkyhat*nefwid*nfactor*k_max full (IFLUX=true) combos.
    (; sfmin=Float64(epo.FACTOR_IN), ky=Float64(epo.KYMARK), w=Float64(epo.WIDTH_IN),
       evals=4*8*8*4, wall=t)
end

function run_ad(ep, prof)
    local res
    t = @elapsed begin
        res = critical_factor_optimize(ep, prof; faithful_confirm=true, extend_width=true,
                    extend_mode=:locate, scan_lo=Float64(ep.FACTOR_IN)/512.0,
                    inner=INNER, use_gpu=USE_GPU)
    end
    sf = (res.faithful !== nothing && res.faithful.binding != :none &&
          isfinite(res.faithful.factor_faithful)) ? res.faithful.factor_faithful : res.sfmin
    (; sfmin=sf, ky=res.kyhat, w=res.width, evals=res.evals, wall=t)
end

function run_multistart(ep, prof)
    lalgo = Symbol(get(ENV, "NLS_LOCAL_ALGO", "LN_BOBYQA"))
    local res
    t = @elapsed begin
        res = critical_factor_multistart(ep, prof; local_algo=lalgo,
                    scan_lo=Float64(ep.FACTOR_IN)/512.0, inner=INNER, use_gpu=USE_GPU)
    end
    (; sfmin=res.sfmin, ky=res.kyhat, w=res.width,
       evals=res.total_evals_full + res.total_evals_eig, wall=t)
end

function run_nlopt(ep, prof)
    algo = Symbol(get(ENV, "NLOPT_ALGO", "GN_DIRECT_L"))
    mev  = parse(Int, get(ENV, "NLOPT_MAXEVAL", "40"))
    local res
    t = @elapsed begin
        res = critical_factor_nlopt(ep, prof; algo=algo, max_evals=mev,
                    scan_lo=Float64(ep.FACTOR_IN)/512.0, inner=INNER, use_gpu=USE_GPU)
    end
    (; sfmin=res.sfmin, ky=res.kyhat, w=res.width,
       evals=res.total_evals_full + res.total_evals_eig, wall=t)
end

const RUNNERS = Dict(:grid=>run_grid, :ad=>run_ad, :multistart=>run_multistart,
                     :nlopt=>run_nlopt)
const LABELS  = Dict(:grid=>"grid (wide box)", :ad=>"ad :locate (narrow)",
                     :multistart=>"multistart", :nlopt=>"nlopt")

function main()
    @assert isfile(GAC) && isfile(TGL) "missing case files:\n  $GAC\n  $TGL"
    base_ep, prof, _ = preprocess_gacode_inputs(GAC, TGL)
    scan_n = Int(base_ep.SCAN_N)

    radii = let r = get(ENV, "RADII", "")
        isempty(r) ? collect(1:scan_n) : parse.(Int, split(r, ","))
    end
    solvers = let s = get(ENV, "SOLVERS", "")
        isempty(s) ? [:grid, :ad, :multistart, :nlopt] : Symbol.(split(s, ","))
    end

    @printf("Benchmark: DIII-D 202017C42_500ms_v3.1  nb=%d  use_gpu=%s  radii=%s\n",
            NB, USE_GPU, join(radii, ","))
    println("solvers: ", join(string.(solvers), ", "))

    # results[solver] :: Vector over `radii` of NamedTuples (with .ir added)
    results = Dict{Symbol,Vector{NamedTuple}}(s => NamedTuple[] for s in solvers)
    for i in radii
        ir = base_ep.IR_EXP[i]
        for s in solvers
            r = RUNNERS[s](_prep(base_ep, i), prof)
            push!(results[s], (; ir=ir, r...))
            @printf("  [i=%2d ir=%3d] %-8s sfmin=%10.4g  ky=%.3f  w=%.3f  evals=%6d  wall=%6.1fs\n",
                    i, ir, string(s), r.sfmin, r.ky, r.w, r.evals, r.wall)
            flush(stdout)
        end
    end

    irs = [base_ep.IR_EXP[i] for i in radii]

    # ── comparison table (ratios vs ad:locate and grid where available) ──
    println("\n================ sfmin comparison ================")
    hdr = "  IR " * join([@sprintf("%12s", LABELS[s]) for s in solvers])
    println(hdr)
    for (k, ir) in enumerate(irs)
        row = @sprintf("%4d", ir)
        for s in solvers
            row *= @sprintf("%12.4g", results[s][k].sfmin)
        end
        println(row)
    end

    if haskey(results, :ad)
        println("\n---- ratio to ad:locate (narrow-width benchmark) ----")
        for (k, ir) in enumerate(irs)
            adv = results[:ad][k].sfmin
            row = @sprintf("%4d", ir)
            for s in solvers
                row *= @sprintf("%12.3f", results[s][k].sfmin / adv)
            end
            println(row)
        end
    end

    # ── per-solver summary: median|max ratio to ad:locate, total wall, mean evals ──
    println("\n================ summary ================")
    for s in solvers
        walls = [x.wall for x in results[s]]
        evs   = [x.evals for x in results[s]]
        if haskey(results, :ad) && s != :ad
            ratios = [results[s][k].sfmin / results[:ad][k].sfmin for k in eachindex(irs)
                      if isfinite(results[:ad][k].sfmin) && results[:ad][k].sfmin > 0]
            matchN = count(x -> x <= 1.05, ratios)   # within 5% or lower (i.e. as/more unstable)
            @printf("  %-8s  Σwall=%7.1fs  mean_evals=%6.0f  ratio/ad median=%.3f max=%.3f  #(<=1.05)/N=%d/%d\n",
                    string(s), sum(walls), sum(evs)/length(evs),
                    _median(ratios), maximum(ratios), matchN, length(ratios))
        else
            @printf("  %-8s  Σwall=%7.1fs  mean_evals=%6.0f\n", string(s), sum(walls), sum(evs)/length(evs))
        end
    end

    # ── CSV ──
    csv = joinpath(@__DIR__, "benchmark_nls_solvers_nb$(NB).csv")
    open(csv, "w") do io
        println(io, "ir,solver,sfmin,ky,width,evals,wall_s")
        for s in solvers, x in results[s]
            @printf(io, "%d,%s,%.6g,%.6g,%.6g,%d,%.3f\n", x.ir, s, x.sfmin, x.ky, x.w, x.evals, x.wall)
        end
    end
    println("\nWrote ", csv)

    # ── plot sfmin vs IR ──
    colors = Dict(:grid=>:gray55, :ad=>:seagreen, :multistart=>:darkorange,
                  :nlopt=>:royalblue)
    marks  = Dict(:grid=>:diamond, :ad=>:utriangle, :multistart=>:rect,
                  :nlopt=>:circle)
    default(legendfontsize=9, guidefontsize=11, tickfontsize=9, dpi=200)
    p = plot(xlabel="radial grid index IR", ylabel="sfmin (critical factor)",
             title="Derivative-free solvers vs grid/ad @ nb=$(NB): DIII-D 202017C42_500ms",
             yscale=:log10, legend=:topleft, size=(900, 560), left_margin=4mm, bottom_margin=4mm)
    for s in solvers
        ys = [x.sfmin for x in results[s]]
        plot!(p, irs, ys; label=LABELS[s], color=get(colors, s, :black),
              marker=get(marks, s, :xcross), markersize=5, linewidth=2)
    end
    png = joinpath(@__DIR__, "benchmark_nls_solvers_nb$(NB).png")
    savefig(p, png)
    println("Wrote ", png)
end

_median(v) = isempty(v) ? NaN : (s = sort(v); n = length(s); iseven(n) ? 0.5*(s[n÷2]+s[n÷2+1]) : s[(n+1)÷2])

main()
