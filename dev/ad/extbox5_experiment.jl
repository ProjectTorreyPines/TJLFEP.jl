# Ext-box round 5 — IR=95 ground-truth question: is there a finite minimum as width->0, or a runaway
# into the basis-breakdown corner?
#   (A) sfmin(width) fine sweep at fixed nb=32, ky in {0.5,0.8}, width 0.05..1.0 — does it bottom out
#       or keep dropping until the Hermite overlap matrix (p0/bp) goes singular?
#   (B) finer nbasis sweep {8,16,24,32,40,48,56} at the narrow optimum (0.8,0.1) — plateau within the
#       usable rank before the singular ceiling?
# All evals wrapped in try/catch so a SingularException (rank-deficient basis) is logged, not fatal.
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
const IR     = parse(Int, get(ENV, "IR", "95"))
const WIDTHS = [0.05, 0.075, 0.1, 0.125, 0.15, 0.2, 0.3, 0.5, 0.75, 1.0]
const KYS    = [0.5, 0.8]
const NBS    = [8, 16, 24, 32, 40, 48, 56]

# faithful sfmin at one (ky,w,nb); returns (sfmin, binding-string)
function faithful_pt(opts, prof, ir, ky, w, nb; team, use_gpu)
    ep = deepcopy(opts); ep.IR = ir; ep.N_BASIS = nb
    gth = TJLFEP._gamma_thresh_for(ep, prof); shi = Float64(ep.FACTOR_IN); slo = shi/512.0
    try
        r = TJLFEP.marginal_factor_faithful(ep, prof; kyhat=ky, width=w, gamma_thresh=gth,
                scan_lo=slo, scan_hi=shi, threaded=true, inner=INNER, team=team, use_gpu=use_gpu)
        f = (r.binding !== :none && isfinite(r.factor_faithful)) ? r.factor_faithful : Inf
        return (f, String(r.binding))
    catch e
        return (NaN, "SINGULAR")
    end
end

function main()
    opts, prof, _ = preprocess_gacode_inputs(GACODE, TGLFEP)
    opts.N_BASIS = NB
    team = USE_MPS ? workers() : nothing
    dev = USE_GPU ? CUDA.name(first(CUDA.devices())) : "CPU"
    @printf("DIII-D %s  inner=%s team=%s  EXT-BOX-5 IR=%d (width-trend + finer nbasis)\n",
            dev, String(INNER), team===nothing ? "-" : string(length(team)), IR)
    flush(stdout)

    let ; faithful_pt(opts, prof, IR, 0.5, 0.5, NB; team=team, use_gpu=USE_GPU); end  # warmup

    # (A) sfmin(width) at fixed nb=NB
    @printf("\n========== (A) sfmin(width) @ IR=%d, nb=%d ==========\n", IR, NB)
    @printf("  %-8s", "ky\\w"); for w in WIDTHS; @printf("%9.4g", w); end; println()
    for ky in KYS
        @printf("  %-8.4g", ky)
        prev = NaN
        for w in WIDTHS
            f, b = faithful_pt(opts, prof, IR, ky, w, NB; team=team, use_gpu=USE_GPU)
            if b == "SINGULAR"; @printf("%9s", "SING")
            elseif !isfinite(f); @printf("%9s", "·")
            else; @printf("%9.4g", f); end
        end
        println()
    end
    @printf("  (decreasing monotonically toward small w => width->0 runaway; flattening => finite floor)\n")
    flush(stdout)

    # (B) finer nbasis at the narrow optimum
    @printf("\n========== (B) nbasis sweep @ IR=%d (ky=0.8, w=0.1) ==========\n", IR)
    @printf("    %-6s %-12s %-14s\n", "nb", "sfmin", "binding")
    prev = NaN
    for nb in NBS
        f, b = faithful_pt(opts, prof, IR, 0.8, 0.1, nb; team=team, use_gpu=USE_GPU)
        d = (isnan(prev) || isnan(f)) ? "" : @sprintf("Δ=%+.4g", f - prev)
        isfinite(f) && (prev = f)
        @printf("    %-6d %-12.5g %-14s %s\n", nb, f, b, d)
        flush(stdout)
    end

    println("\n=== extbox5 experiment done ===")
end

main()
