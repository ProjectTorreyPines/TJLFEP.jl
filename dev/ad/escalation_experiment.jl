# Validate the adf1->escalation production policy under the CANONICAL ky>=0.25 floor.
# Per radius runs each heavy method ONCE and synthesizes the escalation outcome (no redundant
# wrapper re-runs):
#   adf1      : fast default (+ trust diagnostics cheap_gap / feasible_frac)
#   direct40  : canonical-bounds DIRECT-40  <-- THE UNTESTED CONFIG (refs came from unbounded ky)
#   robust    : canonical grid-zoom (robust_ad) = the :grid escalation target
# Then prints, per radius: the trust flag + what escalate_to=:direct vs :grid would deliver, so we
# can choose the default escalation target on evidence (esp. the sparse IR=95).
#
# Env: USE_GPU(1) NB(32) RADII(22,38,48,95) DIRECT_EVALS(40) KY_LO(0.25) INNER(mps_team) MPS_TEAM(4)

using TJLF, TJLFEP, NLopt, Printf, Distributed

const USE_GPU  = get(ENV, "USE_GPU", "1") == "1"
const INNER    = Symbol(get(ENV, "INNER", "threads"))
const MPS_TEAM = parse(Int, get(ENV, "MPS_TEAM", "0"))
const THREADS_PER_WORKER = parse(Int, get(ENV, "JULIA_WORKER_THREADS", "2"))
const USE_MPS  = USE_GPU && INNER === :mps_team && MPS_TEAM > 0

if USE_GPU
    using CUDA
    @assert CUDA.functional() "USE_GPU=1 but no functional GPU"
end

if USE_MPS
    let root = normpath(@__DIR__, "..", ".."),
        team_gpus = String.(split(get(ENV, "TEAM_GPUS", get(ENV, "CUDA_VISIBLE_DEVICES", "0")), ',', keepempty=false)),
        base_env = Dict{String,String}()
        for k in ("JULIA_DEPOT_PATH", "CUDA_MPS_PIPE_DIRECTORY", "CUDA_MPS_LOG_DIRECTORY",
                  "JULIA_CUDA_USE_COMPAT", "JULIA_CUDA_MEMORY_POOL")
            haskey(ENV, k) && (base_env[k] = ENV[k])
        end
        base_env["JULIA_PKG_PRECOMPILE_AUTO"] = "0"
        for w in 1:MPS_TEAM
            env = copy(base_env)
            env["CUDA_VISIBLE_DEVICES"] = team_gpus[(w - 1) % length(team_gpus) + 1]
            addprocs(1; exeflags=`--project=$(root) -t $(THREADS_PER_WORKER) --startup-file=no`, env=env)
        end
    end
    @everywhere begin
        using CUDA, TJLFEP, TJLF, LinearAlgebra
        BLAS.set_num_threads(1)
        CUDA.functional() && CUDA.device!(first(CUDA.devices()))
    end
end

include(joinpath(@__DIR__, "direct_solver.jl"))

const CASE   = joinpath(@__DIR__, "..", "..", "examples", "DIIID_202017C42_500ms_v3.1")
const GACODE = joinpath(CASE, "input.gacode")
const NB     = parse(Int, get(ENV, "NB", "32"))
const TGLFEP = joinpath(CASE, "input_scan20_nb$(NB).TGLFEP")
const RADII  = parse.(Int, split(get(ENV, "RADII", "22,38,48,95"), ',', keepempty=false))
const D_EVALS = parse(Int, get(ENV, "DIRECT_EVALS", "40"))
const KY_LO  = parse(Float64, get(ENV, "KY_LO", "0.25"))
const GAP_THRESH  = parse(Float64, get(ENV, "GAP_THRESH", "1.5"))
const FEAS_THRESH = parse(Float64, get(ENV, "FEAS_THRESH", "0.25"))

# Canonical-range references: grid (=robust_ad, ky>=0.25) is the trusted truth on hard radii;
# dense/D40 came from UNBOUNDED ky and are shown only for context.
const GRID    = Dict(22 => 0.17577, 38 => 0.01953, 48 => 0.039062, 95 => 2.636713)
const DENSE   = Dict(22 => 0.16246, 38 => 0.019531, 48 => 0.028553, 95 => 4.2007)

pct(x, ref) = (isfinite(x) && isfinite(ref) && ref != 0) ? @sprintf("%+6.1f%%", 100*(x-ref)/ref) : "   n/a "

function main()
    opts, prof, _ = preprocess_gacode_inputs(GACODE, TGLFEP)
    opts.N_BASIS = NB
    team = USE_MPS ? workers() : nothing
    dev = USE_GPU ? CUDA.name(first(CUDA.devices())) : "CPU"
    @printf("DIII-D N_BASIS=%d  %s  inner=%s team=%s  DIRECT-%d  ky_lo=%.3g  gate: gap>%.2g | feas<%.2g\n",
            NB, dev, String(INNER), team===nothing ? "-" : string(length(team)), D_EVALS, KY_LO, GAP_THRESH, FEAS_THRESH)
    flush(stdout)

    let ep = deepcopy(opts); ep.IR = 38
        critical_factor_ad_f1seed(ep, prof; ky_lo=KY_LO, inner=INNER, team=team, use_gpu=USE_GPU)
        critical_factor_direct(ep, prof; max_evals=12, ky_lo=KY_LO, inner=INNER, team=team, use_gpu=USE_GPU)
        critical_factor_robust(ep, prof; inner=INNER, team=team, use_gpu=USE_GPU)
    end

    for ir in RADII
        ep = deepcopy(opts); ep.IR = ir
        g = get(GRID, ir, NaN); d = get(DENSE, ir, NaN)

        ta = @elapsed ra = critical_factor_ad_f1seed(ep, prof; ky_lo=KY_LO, inner=INNER, team=team, use_gpu=USE_GPU)
        td = @elapsed rd = critical_factor_direct(ep, prof; max_evals=D_EVALS, ky_lo=KY_LO, inner=INNER, team=team, use_gpu=USE_GPU)
        tr = @elapsed rr = critical_factor_robust(ep, prof; inner=INNER, team=team, use_gpu=USE_GPU)

        reasons = Symbol[]
        (ra.status === :no_onset)        && push!(reasons, :no_onset)
        (ra.status === :cap)             && push!(reasons, :cap)
        (ra.cheap_gap > GAP_THRESH)      && push!(reasons, :cheap_gap)
        (ra.feasible_frac < FEAS_THRESH) && push!(reasons, :sparse)
        flagged = !isempty(reasons)
        esc_direct = flagged ? min(ra.sfmin, rd.sfmin) : ra.sfmin
        esc_grid   = flagged ? min(ra.sfmin, rr.sfmin) : ra.sfmin

        @printf("\nIR=%-3d  grid=%-9.5g dense=%-9.5g  [refs: grid=canonical truth]\n", ir, g, d)
        @printf("  adf1     sfmin=%-10.5g vsgrid=%s  gap=%.2f feas=%.2f %6.1fs st=%-9s flagged=%s%s\n",
                ra.sfmin, pct(ra.sfmin,g), ra.cheap_gap, ra.feasible_frac, ta, String(ra.status),
                flagged, isempty(reasons) ? "" : " "*string(reasons))
        @printf("  direct40 sfmin=%-10.5g vsgrid=%s  full=%-4d eig=%-5d %6.1fs nconf=%d/%d st=%s\n",
                rd.sfmin, pct(rd.sfmin,g), rd.total_evals_full, rd.total_evals_eig, td, rd.n_confirm, rd.n_samples, String(rd.status))
        @printf("  robust   sfmin=%-10.5g vsgrid=%s  full=%-4d eig=%-5d %6.1fs st=%s\n",
                rr.sfmin, pct(rr.sfmin,g), rr.total_evals_full, rr.total_evals_eig, tr, String(rr.status))
        @printf("  ==> escalate:direct => %-10.5g (%s) | escalate:grid => %-10.5g (%s)\n",
                esc_direct, pct(esc_direct,g), esc_grid, pct(esc_grid,g))
        flush(stdout)
    end
    println("\n=== escalation experiment done ===")
end

main()
