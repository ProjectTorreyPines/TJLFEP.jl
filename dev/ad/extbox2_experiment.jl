# Extended-box round 2: (1) push WIDTH down to 0.2 to find where the narrow-width optimum bottoms out
# (IR=95 pinned at w=0.5 last time), and (2) an EP-DRIVE check at each radius's low-width optimum:
# evaluate the AE-band leading growth gamma_AE(factor) as the EP drive is scaled from ~0 to nominal.
#   EP-DRIVEN  => gamma_AE < gth at factor->0 and crosses above as factor grows (genuine EP onset).
#   BACKGROUND => gamma_AE > gth even at factor->0 (mode unstable without EP; lowering WIDTH_MIN would
#                 then WRONGLY drop the critical factor — a false positive, not an EP threshold).
#
# Env: USE_GPU(1) NB(32) RADII(48,95) INNER(mps_team) MPS_TEAM(4)

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
const RADII  = parse.(Int, split(get(ENV, "RADII", "48,95"), ',', keepempty=false))

const KYHATS = [0.005, 0.01, 0.02, 0.05, 0.1, 0.15, 0.25, 0.5]
const WIDTHS = [0.2, 0.3, 0.4, 0.5, 0.6, 0.75, 0.9, 1.0]
const FACTORS = [1.0e-4, 1.0e-3, 0.01953, 0.05, 0.1, 0.3, 1.0, 3.0, 10.0]   # EP-drive scan

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
    return reshape(res, length(WIDTHS), length(KYHATS)), slo, shi
end

# AE-band leading growth gamma_AE(factor) at fixed (ky,w) — same band the marginal scan tracks.
function ep_drive_curve(ep0, prof, ky, w; team, use_gpu)
    ep = deepcopy(ep0); ep.KYHAT_IN = ky; ep.WIDTH_IN = w
    res = TJLFEP._ad_pmap(i -> begin
            g = TJLFEP._gamma_lead_dfactor(ep, prof, FACTORS[i]; use_gpu=use_gpu, ae_band=true)[1]
            (; f=FACTORS[i], g=Float64(g))
        end, length(FACTORS); inner=(team===nothing ? :threads : :mps_team), team=team)
    return res
end

function main()
    opts, prof, _ = preprocess_gacode_inputs(GACODE, TGLFEP)
    opts.N_BASIS = NB
    team = USE_MPS ? workers() : nothing
    dev = USE_GPU ? CUDA.name(first(CUDA.devices())) : "CPU"
    @printf("DIII-D N_BASIS=%d  %s  inner=%s team=%s  EXT-BOX-2 (narrow width + EP-drive check)\n",
            NB, dev, String(INNER), team===nothing ? "-" : string(length(team)))
    @printf("  kyhat = %s\n  width = %s\n", string(KYHATS), string(WIDTHS))
    flush(stdout)

    let ep = deepcopy(opts); ep.IR = RADII[1]
        TJLFEP.marginal_factor_faithful(ep, prof; kyhat=0.25, width=0.6,
            gamma_thresh=TJLFEP._gamma_thresh_for(ep, prof),
            scan_lo=Float64(ep.FACTOR_IN)/512, scan_hi=Float64(ep.FACTOR_IN),
            threaded=false, inner=:serial, team=nothing, use_gpu=USE_GPU)
    end

    for ir in RADII
        ep = deepcopy(opts); ep.IR = ir
        gth = TJLFEP._gamma_thresh_for(ep, prof)
        t = @elapsed (M, slo, shi) = eval_box(ep, prof; team=team, use_gpu=USE_GPU)
        @printf("\n================= IR=%d  (gth=%.2g scan_lo=%.4g, %.1fs) =================\n", ir, gth, slo, t)
        @printf("  %-6s", "w\\ky")
        for ky in KYHATS; @printf("%9.4g", ky); end
        println()
        for (iw, w) in enumerate(WIDTHS)
            @printf("  %-6.3g", w)
            for ik in 1:length(KYHATS)
                r = M[iw, ik]
                if !isfinite(r.f); @printf("%9s", "·")
                else; mark = r.binding === :ae_band_growth ? " " : "*"; @printf("%8.4g%s", r.f, mark); end
            end
            println()
        end
        flat = vec(M); fin = [r for r in flat if isfinite(r.f)]
        if isempty(fin); println("  (no faithful onset in box)"); continue; end
        b = fin[argmin([r.f for r in fin])]
        on_w_edge = b.w == WIDTHS[1] || b.w == WIDTHS[end]
        on_ky_edge = b.ky == KYHATS[1] || b.ky == KYHATS[end]
        @printf("  MIN sfmin=%.5g at (kyhat=%.4g, width=%.3g) binding=%s  ky_edge=%s w_edge=%s\n",
                b.f, b.ky, b.w, String(b.binding), on_ky_edge, on_w_edge)

        # EP-drive check at the low-width optimum
        cur = ep_drive_curve(ep, prof, b.ky, b.w; team=team, use_gpu=USE_GPU)
        @printf("  EP-drive gamma_AE(factor) @ (ky=%.4g, w=%.3g), gth=%.2g:\n", b.ky, b.w, gth)
        @printf("    factor: %s\n", join([@sprintf("%9.4g", c.f) for c in cur], ""))
        @printf("    g_AE  : %s\n", join([@sprintf("%9.3g", c.g) for c in cur], ""))
        g_lo = cur[1].g; g_hi = cur[end].g
        verdict = (g_lo <= gth && g_hi > gth) ? "EP-DRIVEN (stable at factor->0, unstable at nominal)" :
                  (g_lo > gth) ? "BACKGROUND (unstable even at factor->0 — NOT an EP threshold!)" :
                  "AE-STABLE at nominal (g_hi<=gth)"
        @printf("    => %s\n", verdict)
        flush(stdout)
    end
    println("\n=== extbox2 experiment done ===")
end

main()
