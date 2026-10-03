# Ext-box round 3:
#   (A) IR=95 corner extension — push ky up to 1.0 and width down to 0.1 (last min pinned at ky=0.5,
#       w=0.2 corner) to find where sfmin actually bottoms out.
#   (B) nbasis convergence sweep at the narrow-width optima (IR=48 @ ky=0.25,w=0.6; IR=95 @ Part-A
#       best): sfmin(nbasis) and gamma_AE(nbasis, factor=1) for nbasis in {8,16,24,32}. If flat by
#       ~16, narrow modes are over-resolved at 32 (cheaper is possible) and the values are trustworthy.
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

const KYHATS = [0.25, 0.4, 0.5, 0.65, 0.8, 1.0]    # IR=95 corner: push ky UP
const WIDTHS = [0.1, 0.15, 0.2, 0.3, 0.45, 0.6]     #              push width DOWN
const NBS    = [8, 16, 24, 32]

function eval_box(ep0, prof; team, use_gpu)
    gth = TJLFEP._gamma_thresh_for(ep0, prof)
    shi = Float64(ep0.FACTOR_IN); slo = shi/512.0
    pts = [(ky, w) for ky in KYHATS for w in WIDTHS]
    res = TJLFEP._ad_pmap(idx -> begin
            ky, w = pts[idx]
            r = TJLFEP.marginal_factor_faithful(ep0, prof; kyhat=ky, width=w, gamma_thresh=gth,
                    scan_lo=slo, scan_hi=shi, threaded=false, inner=:serial, team=nothing, use_gpu=use_gpu)
            f = (r.binding !== :none && isfinite(r.factor_faithful)) ? r.factor_faithful : Inf
            (; ky=ky, w=w, f=f, binding=r.binding)
        end, length(pts); inner=(team===nothing ? :threads : :mps_team), team=team)
    return reshape(res, length(WIDTHS), length(KYHATS))
end

function nbasis_sweep(opts, prof, ir, ky, w; team, use_gpu)
    @printf("  nbasis sweep @ IR=%d (ky=%.4g, w=%.3g):\n", ir, ky, w)
    @printf("    %-6s %-12s %-12s %-14s\n", "nb", "sfmin", "g_AE(f=1)", "binding")
    for nb in NBS
        ep = deepcopy(opts); ep.IR = ir; ep.N_BASIS = nb
        gth = TJLFEP._gamma_thresh_for(ep, prof); shi = Float64(ep.FACTOR_IN); slo = shi/512.0
        r = TJLFEP.marginal_factor_faithful(ep, prof; kyhat=ky, width=w, gamma_thresh=gth,
                scan_lo=slo, scan_hi=shi, threaded=true, inner=INNER, team=team, use_gpu=use_gpu)
        epg = deepcopy(ep); epg.KYHAT_IN = ky; epg.WIDTH_IN = w
        g1 = Float64(TJLFEP._gamma_lead_dfactor(epg, prof, 1.0; use_gpu=use_gpu, ae_band=true)[1])
        f = (r.binding !== :none && isfinite(r.factor_faithful)) ? r.factor_faithful : Inf
        @printf("    %-6d %-12.5g %-12.4g %-14s\n", nb, f, g1, String(r.binding))
        flush(stdout)
    end
end

function main()
    opts, prof, _ = preprocess_gacode_inputs(GACODE, TGLFEP)
    opts.N_BASIS = NB
    team = USE_MPS ? workers() : nothing
    dev = USE_GPU ? CUDA.name(first(CUDA.devices())) : "CPU"
    @printf("DIII-D %s  inner=%s team=%s  EXT-BOX-3 (IR95 corner + nbasis convergence)\n",
            dev, String(INNER), team===nothing ? "-" : string(length(team)))
    flush(stdout)

    let ep = deepcopy(opts); ep.IR = 95
        TJLFEP.marginal_factor_faithful(ep, prof; kyhat=0.5, width=0.2,
            gamma_thresh=TJLFEP._gamma_thresh_for(ep, prof),
            scan_lo=Float64(ep.FACTOR_IN)/512, scan_hi=Float64(ep.FACTOR_IN),
            threaded=false, inner=:serial, team=nothing, use_gpu=USE_GPU)
    end

    # (A) IR=95 corner extension
    ep95 = deepcopy(opts); ep95.IR = 95
    @printf("\n========== (A) IR=95 corner extension (nb=%d) ==========\n", NB)
    @printf("  kyhat=%s width=%s\n", string(KYHATS), string(WIDTHS))
    t = @elapsed M = eval_box(ep95, prof; team=team, use_gpu=USE_GPU)
    @printf("  %-6s", "w\\ky"); for ky in KYHATS; @printf("%9.4g", ky); end; println()
    for (iw, w) in enumerate(WIDTHS)
        @printf("  %-6.3g", w)
        for ik in 1:length(KYHATS)
            r = M[iw, ik]
            if !isfinite(r.f); @printf("%9s", "·")
            else; mark = r.binding === :ae_band_growth ? " " : "*"; @printf("%8.4g%s", r.f, mark); end
        end
        println()
    end
    fin = [r for r in vec(M) if isfinite(r.f)]
    b95 = fin[argmin([r.f for r in fin])]
    @printf("  MIN sfmin=%.5g at (kyhat=%.4g, width=%.3g) binding=%s  ky_edge=%s w_edge=%s  (%.1fs)\n",
            b95.f, b95.ky, b95.w, String(b95.binding),
            b95.ky == KYHATS[1] || b95.ky == KYHATS[end], b95.w == WIDTHS[1] || b95.w == WIDTHS[end], t)

    # (B) nbasis convergence at the narrow-width optima
    println("\n========== (B) nbasis convergence ==========")
    nbasis_sweep(opts, prof, 48, 0.25, 0.6; team=team, use_gpu=USE_GPU)
    nbasis_sweep(opts, prof, 95, b95.ky, b95.w; team=team, use_gpu=USE_GPU)

    println("\n=== extbox3 experiment done ===")
end

main()
