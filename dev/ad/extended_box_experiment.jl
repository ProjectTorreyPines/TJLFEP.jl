# Extended-box PHYSICS diagnostic: evaluate the FAITHFUL marginal factor on a (kyhat,width) mesh that
# deliberately steps OUTSIDE the input/grid bounds — kyhat far below the grid's first node (0.25) and
# down past 0.001, width below WIDTH_MIN=1.0 and above WIDTH_MAX=2.0. Answers:
#   - does sfmin keep dropping as kyhat->0, or stabilize (is there a real physical low-ky floor)?
#   - does the width optimum sit below 1.0 (is the WIDTH_MIN=1 floor cutting the true minimum)?
#   - do the modes stay genuine AE-band instabilities out there (binding==ae_band_growth) or break down?
# marginal_factor_faithful sets KYHAT_IN/WIDTH_IN directly (KY_MODEL forced 3 => KY ∝ kyhat); the
# single-point eval is NOT clamped to [WIDTH_MIN,WIDTH_MAX], so we can probe outside the box.
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

# Extended mesh — steps OUTSIDE the [0.25..1]x[1..2] grid box on purpose.
const KYHATS = [0.0005, 0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.15, 0.25, 0.5, 0.75, 1.0]
const WIDTHS = [0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5]

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

function main()
    opts, prof, _ = preprocess_gacode_inputs(GACODE, TGLFEP)
    opts.N_BASIS = NB
    team = USE_MPS ? workers() : nothing
    dev = USE_GPU ? CUDA.name(first(CUDA.devices())) : "CPU"
    @printf("DIII-D N_BASIS=%d  %s  inner=%s team=%s  EXTENDED BOX (grid box = ky[0.25..1] x w[1..2])\n",
            NB, dev, String(INNER), team===nothing ? "-" : string(length(team)))
    @printf("  kyhat = %s\n  width = %s\n", string(KYHATS), string(WIDTHS))
    flush(stdout)

    # warmup
    let ep = deepcopy(opts); ep.IR = RADII[1]
        TJLFEP.marginal_factor_faithful(ep, prof; kyhat=0.25, width=1.5,
            gamma_thresh=TJLFEP._gamma_thresh_for(ep, prof),
            scan_lo=Float64(ep.FACTOR_IN)/512, scan_hi=Float64(ep.FACTOR_IN),
            threaded=false, inner=:serial, team=nothing, use_gpu=USE_GPU)
    end

    for ir in RADII
        ep = deepcopy(opts); ep.IR = ir
        t = @elapsed (M, slo, shi) = eval_box(ep, prof; team=team, use_gpu=USE_GPU)
        @printf("\n================= IR=%d  (scan_lo=%.4g scan_hi=%.4g, %.1fs) =================\n", ir, slo, shi, t)
        # sfmin matrix: rows=width, cols=kyhat. '*' marks a non-AE-band binding (mode changed character);
        # 'pin' (Inf) marks no faithful onset.
        @printf("  %-6s", "w\\ky")
        for ky in KYHATS; @printf("%9.4g", ky); end
        println()
        for (iw, w) in enumerate(WIDTHS)
            @printf("  %-6.3g", w)
            for (ik, _) in enumerate(KYHATS)
                r = M[iw, ik]
                if !isfinite(r.f)
                    @printf("%9s", "·")
                else
                    mark = r.binding === :ae_band_growth ? " " : "*"
                    @printf("%8.4g%s", r.f, mark)
                end
            end
            println()
        end
        # global min + boundary diagnostics
        flat = vec(M); fin = [r for r in flat if isfinite(r.f)]
        if isempty(fin)
            println("  (no faithful onset anywhere in extended box)")
        else
            b = fin[argmin([r.f for r in fin])]
            on_ky_edge = b.ky == KYHATS[1] || b.ky == KYHATS[end]
            on_w_edge  = b.w  == WIDTHS[1]  || b.w  == WIDTHS[end]
            @printf("  MIN sfmin=%.5g at (kyhat=%.4g, width=%.3g) binding=%s  ky_edge=%s w_edge=%s\n",
                    b.f, b.ky, b.w, String(b.binding), on_ky_edge, on_w_edge)
            # low-ky trend at the min's width row: is sfmin still dropping toward ky->0?
            iw = findfirst(==(b.w), WIDTHS)
            row = [M[iw, ik].f for ik in 1:length(KYHATS)]
            @printf("  low-ky trend @w=%.3g: %s\n", b.w, join([@sprintf("%.3g", x) for x in row], " "))
            n_nonae = count(r -> isfinite(r.f) && r.binding !== :ae_band_growth, flat)
            @printf("  non-AE-band feasible nodes: %d / %d\n", n_nonae, length(fin))
        end
        flush(stdout)
    end
    println("\n=== extended-box experiment done ===")
end

main()
