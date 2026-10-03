# Ext-box round 4: higher-nbasis convergence at the narrow-width optima to distinguish
#   slow-convergence (finite limit => grid genuinely overestimates) vs runaway (->0 => artifact,
#   WIDTH_MIN=1 justified). nbasis in {32,48,64,96}. Includes a STANDARD-width control (w=1.5) to
#   confirm nb=32 is already converged there (the regime the production scan uses).
#
# Env: USE_GPU(1) NB(32) INNER(mps_team) MPS_TEAM(4)

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

const CASE   = joinpath(@__DIR__, "..", "..", "examples", "DIIID_202017C42_500ms_v3.1")
const GACODE = joinpath(CASE, "input.gacode")
const NB     = parse(Int, get(ENV, "NB", "32"))
const TGLFEP = joinpath(CASE, "input_scan20_nb$(NB).TGLFEP")
const NBS    = [32, 48, 64, 96]

function nbasis_sweep(opts, prof, label, ir, ky, w; team, use_gpu)
    @printf("  %-26s @ IR=%d (ky=%.4g, w=%.3g):\n", label, ir, ky, w)
    @printf("    %-6s %-12s %-12s %-14s %-8s\n", "nb", "sfmin", "g_AE(f=1)", "binding", "secs")
    prev = NaN
    for nb in NBS
        ep = deepcopy(opts); ep.IR = ir; ep.N_BASIS = nb
        gth = TJLFEP._gamma_thresh_for(ep, prof); shi = Float64(ep.FACTOR_IN); slo = shi/512.0
        local f, g1, bind, t
        try
            t = @elapsed r = TJLFEP.marginal_factor_faithful(ep, prof; kyhat=ky, width=w, gamma_thresh=gth,
                    scan_lo=slo, scan_hi=shi, threaded=true, inner=INNER, team=team, use_gpu=use_gpu)
            f = (r.binding !== :none && isfinite(r.factor_faithful)) ? r.factor_faithful : Inf
            bind = String(r.binding)
        catch e
            # SingularException / ill-conditioned basis at this (nb,width) — itself a convergence signal.
            f = NaN; bind = "SINGULAR/"*string(typeof(e)); t = 0.0
        end
        g1 = NaN
        try
            epg = deepcopy(ep); epg.KYHAT_IN = ky; epg.WIDTH_IN = w
            g1 = Float64(TJLFEP._gamma_lead_dfactor(epg, prof, 1.0; use_gpu=use_gpu, ae_band=true)[1])
        catch; end
        d = (isnan(prev) || isnan(f)) ? "" : @sprintf("Δ=%+.4g", f - prev)
        isfinite(f) && (prev = f)
        @printf("    %-6d %-12.5g %-12.4g %-22s %-8.1f %s\n", nb, f, g1, bind, t, d)
        flush(stdout)
    end
end

function main()
    opts, prof, _ = preprocess_gacode_inputs(GACODE, TGLFEP)
    opts.N_BASIS = NB
    team = USE_MPS ? workers() : nothing
    dev = USE_GPU ? CUDA.name(first(CUDA.devices())) : "CPU"
    @printf("DIII-D %s  inner=%s team=%s  EXT-BOX-4 (high-nbasis convergence) nbs=%s\n",
            dev, String(INNER), team===nothing ? "-" : string(length(team)), string(NBS))
    flush(stdout)

    let ep = deepcopy(opts); ep.IR = 48
        TJLFEP.marginal_factor_faithful(ep, prof; kyhat=0.25, width=0.6,
            gamma_thresh=TJLFEP._gamma_thresh_for(ep, prof),
            scan_lo=Float64(ep.FACTOR_IN)/512, scan_hi=Float64(ep.FACTOR_IN),
            threaded=false, inner=:serial, team=nothing, use_gpu=USE_GPU)
    end

    println()
    nbasis_sweep(opts, prof, "IR48 NARROW (w=0.6)", 48, 0.25, 0.6;  team=team, use_gpu=USE_GPU)
    nbasis_sweep(opts, prof, "IR48 STANDARD (w=1.5)", 48, 0.25, 1.5; team=team, use_gpu=USE_GPU)
    nbasis_sweep(opts, prof, "IR95 NARROW (w=0.1)", 95, 0.8, 0.1;   team=team, use_gpu=USE_GPU)
    println("\n=== extbox4 experiment done ===")
end

main()
