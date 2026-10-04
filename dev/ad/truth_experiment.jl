# Validation of the extended-box + separable-nbasis "truth" protocol and its triggered production
# wrapper, over four representative radii on DIII-D n_scan=20:
#   IR=22  clean              — trigger should NOT fire; truth must match canonical (no regression)
#   IR=38  floor-pinned       — trigger may fire (width floor); truth should match canonical onset
#   IR=48  interior-hard      — expect narrow-width bowl; truth ≈ 0.0195 (≪ canonical w≥1)
#   IR=95  sparse-hard        — expect narrow-width bowl; truth ≈ 0.21 (≪ canonical ~2.64)
#
# For each radius prints: canonical adf1 (w∈[1,2]) sfmin/(ky,w); trigger decision+reasons; extended
# truth sfmin_work/(ky,w); nbasis table + geometric extrapolation; eig-eval counts; per-stage wall time.
#
# Env: USE_GPU(1) NB(32) INNER(mps_team) MPS_TEAM(4) RADII(22,38,48,95)

using TJLF, TJLFEP, Printf, Distributed

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

# critical_factor_ad_f1seed / critical_factor_truth / critical_factor_triggered are now in core
# (exported by TJLFEP); no local include needed.

const CASE   = joinpath(@__DIR__, "..", "..", "examples", "DIIID_202017C42_500ms_v3.1")
const GACODE = joinpath(CASE, "input.gacode")
const NB     = parse(Int, get(ENV, "NB", "32"))
const TGLFEP = joinpath(CASE, "input_scan20_nb$(NB).TGLFEP")
const RADII  = parse.(Int, split(get(ENV, "RADII", "22,38,48,95"), ','))

function run_radius(opts, prof, ir; team, use_gpu)
    ep = deepcopy(opts); ep.IR = ir; ep.N_BASIS = NB
    gth = TJLFEP._gamma_thresh_for(ep, prof); shi = Float64(ep.FACTOR_IN); slo = shi/512.0
    kw = (; gamma_thresh=gth, scan_lo=slo, scan_hi=shi, inner=INNER, team=team, use_gpu=use_gpu)

    @printf("\n========================= IR=%d =========================\n", ir); flush(stdout)

    # (1) canonical adf1 on w∈[WIDTH_MIN,WIDTH_MAX] (current production answer + trigger diagnostics)
    tc = @elapsed base = critical_factor_ad_f1seed(ep, prof; ky_lo=0.25, nseed_ky=4, nseed_w=8,
            n_eig_seed=12, kw...)
    @printf("  canonical adf1 : sfmin=%.5g  (ky=%.3g, w=%.3g)  bind=%s  feas=%.2f cheap_gap=%.2g  [%.1fs]\n",
            base.sfmin, base.kyhat, base.width, String(base.binding), base.feasible_frac, base.cheap_gap, tc)
    flush(stdout)

    # trigger decision (same logic as critical_factor_triggered)
    wmin = Float64(ep.WIDTH_MIN)
    reasons = Symbol[]
    (isfinite(base.width) && base.width <= wmin*1.05) && push!(reasons, :width_floor)
    (base.status === :no_onset)        && push!(reasons, :no_onset)
    (base.status === :cap)             && push!(reasons, :cap)
    (base.cheap_gap > 1.5)             && push!(reasons, :cheap_gap)
    (base.feasible_frac < 0.25)        && push!(reasons, :sparse)
    @printf("  trigger        : %s  reasons=%s\n", isempty(reasons) ? "NO" : "YES", string(reasons))
    flush(stdout)

    # (2) extended truth (run unconditionally here so we can verify clean radii don't change)
    tt = @elapsed tr = critical_factor_truth(ep, prof; ky_lo=0.05, w_lo=0.05,
            nseed_ky=6, nseed_w=10, n_eig_seed=12, k_descend=6, nb_steps=[32,40,48], kw...)
    @printf("  extended truth : sfmin_work=%.5g  (ky=%.3g, w=%.3g)  bind=%s  [%.1fs]\n",
            tr.sfmin_work, tr.kyhat, tr.width, String(tr.binding), tt)
    nbs, vals = tr.nb_table
    @printf("    nbasis: "); for i in eachindex(nbs); @printf("nb%d=%.5g ", nbs[i], vals[i]); end; println()
    @printf("    nb→limit=%.5g  ratio=%.3g  converged=%s\n", tr.nb_limit, tr.nb_ratio, tr.nb_converged)

    # production answer = critical_factor_triggered's policy: min(canonical, truth) when flagged.
    prod = (!isempty(reasons) && isfinite(tr.sfmin)) ? min(base.sfmin, tr.sfmin) : base.sfmin
    prod_src = (prod == base.sfmin) ? :canonical : :truth
    ratio = (isfinite(prod) && isfinite(base.sfmin) && prod>0) ? base.sfmin/prod : NaN
    @printf("  => PRODUCTION sfmin=%.5g (%s)   canonical/production = %.3g×\n", prod, String(prod_src), ratio)
    @printf("  evals: canonical eig=%d full=%d | truth eig=%d full=%d\n",
            base.total_evals_eig, base.total_evals_full, tr.total_evals_eig, tr.total_evals_full)
    flush(stdout)
    return (; ir=ir, base=base.sfmin, base_w=base.width, truth=tr.sfmin, truth_w=tr.width,
            prod=prod, prod_src=prod_src, nb_limit=tr.nb_limit, nb_conv=tr.nb_converged,
            trig=!isempty(reasons), reasons=reasons, t_canon=tc, t_truth=tt)
end

function main()
    opts, prof, _ = preprocess_gacode_inputs(GACODE, TGLFEP)
    opts.N_BASIS = NB
    team = USE_MPS ? workers() : nothing
    dev = USE_GPU ? CUDA.name(first(CUDA.devices())) : "CPU"
    @printf("DIII-D %s  inner=%s team=%s  NB=%d  TRUTH-PROTOCOL validation  radii=%s\n",
            dev, String(INNER), team===nothing ? "-" : string(length(team)), NB, string(RADII))
    flush(stdout)

    # warmup (compile both paths)
    let ep = deepcopy(opts); ep.IR = RADII[1]
        critical_factor_ad_f1seed(ep, prof; ky_lo=0.25, nseed_ky=2, nseed_w=2, n_eig_seed=8,
            inner=INNER, team=team, use_gpu=USE_GPU)
    end

    rows = NamedTuple[]
    for ir in RADII; push!(rows, run_radius(opts, prof, ir; team=team, use_gpu=USE_GPU)); end

    println("\n===================== SUMMARY =====================")
    @printf("  %-5s %-10s %-8s %-10s %-8s %-10s %-9s %-6s %s\n",
            "IR","canon","w*","truth","w*","PROD","src","conv","reasons")
    for r in rows
        @printf("  %-5d %-10.5g %-8.3g %-10.5g %-8.3g %-10.5g %-9s %-6s %s\n",
                r.ir, r.base, r.base_w, r.truth, r.truth_w, r.prod, String(r.prod_src), r.nb_conv, string(r.reasons))
    end
    @printf("\n  wallclock: canonical Σ=%.1fs  truth Σ=%.1fs\n",
            sum(r.t_canon for r in rows), sum(r.t_truth for r in rows))
    println("\n=== truth experiment done ===")
end

main()
