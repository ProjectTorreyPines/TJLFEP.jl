# Head-to-head for the fast-AND-accurate (ky,w) critical-factor solver, scored against the DIRECT-40
# accuracy ceiling (job 54490874):
#   direct20 : critical_factor_direct(max_evals=20)        — lean DIRECT (option 2: cut DIRECT's cost)
#   adf1     : critical_factor_ad_f1seed                   — f1 seed grid + :ad descent (option 1)
# Both faithful. Question: which reaches DIRECT-40 accuracy (esp. the off-node spikes IR=48/95) at
# the lowest cost? eig COUNT is the hardware-independent cost; wallclock is reported under the run's
# INNER (default mps_team=4, matching how DIRECT is timed best).
#
# Env: USE_GPU (0/1, default 1), NB (default 32), RADII (csv), DIRECT20_EVALS (default 20),
#      ADF1_NSEED_KY/ADF1_NSEED_W/ADF1_NEIG (seed grid), INNER (threads|mps_team), MPS_TEAM.

using TJLF
using TJLFEP
using NLopt
using Printf
using Distributed

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
const D20_EVALS = parse(Int, get(ENV, "DIRECT20_EVALS", "20"))
const A_NKY  = parse(Int, get(ENV, "ADF1_NSEED_KY", "4"))
const A_NW   = parse(Int, get(ENV, "ADF1_NSEED_W", "8"))
const A_NEIG = parse(Int, get(ENV, "ADF1_NEIG", "12"))
const KY_LO  = parse(Float64, get(ENV, "KY_LO", "0.25"))   # canonical kwscale_scan kyhat floor

# References (nb=32) from the DIRECT experiment job 54490874.
const GRID    = Dict(22 => 0.17577, 38 => 0.01953, 48 => 0.039062, 95 => 2.636713)
const DENSE   = Dict(22 => 0.16246, 38 => 0.019531, 48 => 0.028553, 95 => 4.2007)
const DIRECT40= Dict(22 => 0.16476, 38 => 0.019531, 48 => 0.026728, 95 => 1.3812)

pct(x, ref) = (isfinite(x) && isfinite(ref) && ref != 0) ? @sprintf("%+6.1f%%", 100*(x-ref)/ref) : "   n/a "

function main()
    opts, prof, _ = preprocess_gacode_inputs(GACODE, TGLFEP)
    opts.N_BASIS = NB
    team = USE_MPS ? workers() : nothing
    dev = USE_GPU ? CUDA.name(first(CUDA.devices())) : "CPU"
    @printf("DIII-D N_BASIS=%d  %s  head-to-head  inner=%s team=%s  direct20=%d  adf1=%dx%d neig=%d  ky_lo=%.3g\n",
            NB, dev, String(INNER), team===nothing ? "-" : string(length(team)), D20_EVALS, A_NKY, A_NW, A_NEIG, KY_LO)
    flush(stdout)

    let ep = deepcopy(opts); ep.IR = 38
        critical_factor_direct(ep, prof; max_evals=12, ky_lo=KY_LO, inner=INNER, team=team, use_gpu=USE_GPU)
        critical_factor_ad_f1seed(ep, prof; nseed_ky=A_NKY, nseed_w=A_NW, n_eig_seed=A_NEIG,
                                  ky_lo=KY_LO, inner=INNER, team=team, use_gpu=USE_GPU)
    end

    for ir in RADII
        ep = deepcopy(opts); ep.IR = ir
        g = get(GRID, ir, NaN); d = get(DENSE, ir, NaN); d40 = get(DIRECT40, ir, NaN)

        t2 = @elapsed r2 = critical_factor_direct(ep, prof; max_evals=D20_EVALS, ky_lo=KY_LO,
                inner=INNER, team=team, use_gpu=USE_GPU)
        t1 = @elapsed r1 = critical_factor_ad_f1seed(ep, prof; nseed_ky=A_NKY, nseed_w=A_NW,
                n_eig_seed=A_NEIG, ky_lo=KY_LO, inner=INNER, team=team, use_gpu=USE_GPU)

        @printf("\nIR=%-3d  grid=%-9.5g dense=%-9.5g direct40=%-9.5g\n", ir, g, d, d40)
        @printf("  %-8s sfmin=%-10.5g  vsgrid=%s vsdense=%s vsD40=%s  full=%-4d eig=%-5d %6.1fs  %-14s nconf=%d/%d st=%s\n",
                "direct20", r2.sfmin, pct(r2.sfmin,g), pct(r2.sfmin,d), pct(r2.sfmin,d40),
                r2.total_evals_full, r2.total_evals_eig, t2, String(r2.binding), r2.n_confirm, r2.n_samples, String(r2.status))
        @printf("  %-8s sfmin=%-10.5g  vsgrid=%s vsdense=%s vsD40=%s  full=%-4d eig=%-5d %6.1fs  %-14s ndesc=%d nconf=%d st=%s\n",
                "adf1", r1.sfmin, pct(r1.sfmin,g), pct(r1.sfmin,d), pct(r1.sfmin,d40),
                r1.total_evals_full, r1.total_evals_eig, t1, String(r1.binding), r1.n_descend, r1.n_confirm, String(r1.status))
        flush(stdout)
    end
    println("\n=== headtohead experiment done ===")
end

main()
