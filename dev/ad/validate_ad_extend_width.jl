# Sanity check for the width-extended `:ad` path (critical_factor_optimize extend_width=true).
#
# The canonical `:ad` descent only searches w∈[WIDTH_MIN,WIDTH_MAX] (w≥1) and so biases sfmin high
# at near-marginal radii whose true critical factor sits at narrow width (w≪1). extend_width=true
# folds in the SAME `_locate_extended` log-w locate `:robust_ad` uses (seeded on the cheap :ad
# descent, no faithful grid). For a narrow-mode radius (IR=95) and a cleaner radius we compare:
#   pure :ad (extend_width=false)  vs  extended :ad (extend_width=true)  vs  :robust_ad (reference).
# Expectations (NOT a strict-parity test — extended :ad is meant to be cheaper & less accurate than
# robust_ad): at the narrow radius extended-:ad sfmin ≪ pure-:ad sfmin and is in the ballpark of
# robust_ad (≤ ~2× of it), at far fewer eigensolves than robust_ad; clean radius ≈ agrees.
#
#   module load julia/1.11.7 && export JULIA_DEPOT_PATH=$PSCRATCH/.julia
#   julia -t 16 --project=. dev/ad/validate_ad_extend_width.jl     # env: RADII=22,95 NB=8

using TJLF, TJLFEP, Printf

const USE_GPU = get(ENV, "USE_GPU", "0") == "1"
if USE_GPU
    using CUDA
    @assert CUDA.functional() "USE_GPU=1 but no functional GPU"
    CUDA.device!(first(CUDA.devices()))
end

const AX_CASE   = joinpath(@__DIR__, "..", "..", "examples", "DIIID_202017C42_500ms_v3.1")
const AX_GACODE = joinpath(AX_CASE, "input.gacode")
const AX_RADII  = parse.(Int, split(get(ENV, "RADII", "22,95"), ','))
const AX_NB     = parse(Int, get(ENV, "NB", "8"))

function run_case(opts, prof, ir)
    ep = deepcopy(opts); ep.IR = ir
    gth = TJLFEP._gamma_thresh_for(ep, prof); shi = Float64(ep.FACTOR_IN); slo = shi / 512.0
    kw = (; gamma_thresh = gth, scan_lo = slo, scan_hi = shi, inner = :threads, team = nothing, use_gpu = USE_GPU)

    @printf("\n===== IR=%d  NB=%d  (FACTOR_IN=%.4g, WIDTH_MIN=%.3g) =====\n",
            ir, ep.N_BASIS, shi, Float64(ep.WIDTH_MIN)); flush(stdout)

    pure = critical_factor_optimize(ep, prof; faithful_confirm = true, extend_width = false, kw...)
    ext  = critical_factor_optimize(ep, prof; faithful_confirm = true, extend_width = true,  kw...)
    rob  = critical_factor_robust(ep, prof; extend_width = true, kw...)

    pf = (pure.faithful !== nothing && pure.faithful.binding != :none) ? pure.faithful.factor_faithful : pure.sfmin
    ef = (ext.faithful  !== nothing && ext.faithful.binding  != :none) ? ext.faithful.factor_faithful  : ext.sfmin
    @printf("  pure :ad (w≥1)    sfmin=%.5g  at (ky=%.4f, w=%.4f)  evals=%d\n", pf, pure.kyhat, pure.width, pure.evals)
    @printf("  ext  :ad (w-ext)  sfmin=%.5g  at (ky=%.4f, w=%.4f)  evals=%d  n_ext_confirm=%d\n",
            ef, ext.kyhat, ext.width, ext.evals, ext.n_ext_confirm)
    @printf("  robust_ad (ref)   sfmin=%.5g  at (ky=%.4f, w=%.4f)  evals_full=%d evals_eig=%d\n",
            rob.sfmin, rob.kyhat, rob.width, rob.total_evals_full, rob.total_evals_eig)
    @printf("    ext/pure = %.3f   ext/robust = %.3f   robust_evals/ext_evals = %.2fx\n",
            ef / max(pf, eps()), ef / max(rob.sfmin, eps()),
            (rob.total_evals_full + rob.total_evals_eig) / max(ext.evals, 1))
    flush(stdout)
    return (; ir, pure = pf, ext = ef, robust = rob.sfmin, ext_w = ext.width,
            ext_evals = ext.evals, rob_evals = rob.total_evals_full + rob.total_evals_eig)
end

function main()
    @printf("width-extended :ad sanity (%s)  radii=%s  NB=%d  threads=%d\n",
            USE_GPU ? "GPU" : "CPU", string(AX_RADII), AX_NB, Threads.nthreads()); flush(stdout)
    tglfep = joinpath(AX_CASE, "input_scan20_nb$(AX_NB).TGLFEP")
    opts, prof, _ = preprocess_gacode_inputs(AX_GACODE, tglfep)
    opts.N_BASIS = AX_NB

    rows = NamedTuple[]
    for ir in AX_RADII
        push!(rows, run_case(opts, prof, ir))
    end

    println("\n===================== SUMMARY =====================")
    @printf("  %-5s %-12s %-12s %-12s %-8s %-8s %-8s\n",
            "IR", "pure(w≥1)", "ext(w-ext)", "robust", "ext_w", "ext_ev", "rob_ev")
    for r in rows
        @printf("  %-5d %-12.5g %-12.5g %-12.5g %-8.3g %-8d %-8d\n",
                r.ir, r.pure, r.ext, r.robust, r.ext_w, r.ext_evals, r.rob_evals)
    end
end

main()
