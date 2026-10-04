# Pure-AD vs confirm-AD width-extended `:ad` comparison.
#
# The production width-extended `:ad` (faithful_confirm=true) faithful-confirms (IFLUX=true) the
# located narrow candidates, so its sfmin is the faithful keep onset and matches :robust_ad bitwise
# at the located mode. The "pure AD" variant (faithful_confirm=false) stays on the cheap AE-band
# onset surface end-to-end (NO IFLUX=true confirm anywhere) — much cheaper, but reports the AE-band
# onset (≤ the faithful keep onset) so it does NOT match :robust_ad bitwise and can dip below it.
#
# Over the 20-radius DIII-D scan we compare, on identical hardware:
#   confirm-ext :ad  (faithful_confirm=true, extend_width=true)   ← production
#   pure-ext    :ad  (faithful_confirm=false, extend_width=true)  ← this experiment
#   robust_ad        (extend_width=true)                          ← faithful reference
# recording sfmin, wallclock, and eigensolve counts for each.
#
#   module load julia/1.11.7 && export JULIA_DEPOT_PATH=$PSCRATCH/.julia
#   USE_GPU=1 NB=8 julia --project=. dev/ad/profile_ad_pure_vs_confirm.jl
# env: RADII (default all 20), NB (default 8), MODES=pure,confirm,robust (subset to skip robust)

using TJLF, TJLFEP, Printf

const USE_GPU = get(ENV, "USE_GPU", "0") == "1"
if USE_GPU
    using CUDA
    @assert CUDA.functional() "USE_GPU=1 but no functional GPU"
    CUDA.device!(first(CUDA.devices()))
end

const AX_CASE   = joinpath(@__DIR__, "..", "..", "examples", "DIIID_202017C42_500ms_v3.1")
const AX_GACODE = joinpath(AX_CASE, "input.gacode")
const AX_NB     = parse(Int, get(ENV, "NB", "8"))
const AX_MODES  = Set(split(get(ENV, "MODES", "pure,confirm,robust"), ','))

faithful_of(r) = (r.faithful !== nothing && r.faithful.binding != :none) ? r.faithful.factor_faithful : r.sfmin

function run_case(opts, prof, ir)
    ep = deepcopy(opts); ep.IR = ir
    gth = TJLFEP._gamma_thresh_for(ep, prof); shi = Float64(ep.FACTOR_IN); slo = shi / 512.0
    kw = (; gamma_thresh = gth, scan_lo = slo, scan_hi = shi, inner = :threads, team = nothing, use_gpu = USE_GPU)

    @printf("\n===== IR=%d  NB=%d  (FACTOR_IN=%.4g, WIDTH_MIN=%.3g) =====\n",
            ir, ep.N_BASIS, shi, Float64(ep.WIDTH_MIN)); flush(stdout)

    row = Dict{Symbol,Any}(:ir => ir)

    if "orig" in AX_MODES
        # The ORIGINAL fast :ad: single w≥1 descent + 1 faithful confirm, NO _locate_extended grid.
        t = @elapsed o = critical_factor_optimize(ep, prof; faithful_confirm = true, extend_width = false, kw...)
        of = faithful_of(o)
        row[:orig] = of; row[:orig_w] = o.width; row[:orig_ev] = o.evals; row[:orig_t] = t
        @printf("  orig :ad (w≥1)    sfmin=%.5g  at (ky=%.4f, w=%.4f)  evals=%d  (no locate grid)  %.1fs\n",
                of, o.kyhat, o.width, o.evals, t); flush(stdout)
    end
    if "wide" in AX_MODES
        # The "why not just widen the box?" test: the SAME cheap solver (single seed-scan + ONE
        # projected-gradient descent + 1 faithful confirm), but with the descent/seed box widened to
        # w∈[0.05, WIDTH_MAX] and denser seeds — extend_width=false, so NO multi-descent locate grid.
        wmax = Float64(ep.WIDTH_MAX)
        t = @elapsed wd = critical_factor_optimize(ep, prof; faithful_confirm = true, extend_width = false,
                w_bounds = (0.05, wmax), nseed_ky = 6, nseed_w = 10, kw...)
        wf = faithful_of(wd)
        row[:wide] = wf; row[:wide_w] = wd.width; row[:wide_ev] = wd.evals; row[:wide_t] = t
        @printf("  wide :ad (1 desc) sfmin=%.5g  at (ky=%.4f, w=%.4f)  evals=%d  (extended box, single descent)  %.1fs\n",
                wf, wd.kyhat, wd.width, wd.evals, t); flush(stdout)
    end
    if "confirm" in AX_MODES
        t = @elapsed c = critical_factor_optimize(ep, prof; faithful_confirm = true, extend_width = true, kw...)
        cf = faithful_of(c)
        row[:confirm] = cf; row[:confirm_w] = c.width; row[:confirm_ev] = c.evals; row[:confirm_t] = t
        @printf("  confirm-ext :ad   sfmin=%.5g  at (ky=%.4f, w=%.4f)  evals=%d  n_ext=%d  %.1fs\n",
                cf, c.kyhat, c.width, c.evals, c.n_ext_confirm, t); flush(stdout)
    end
    if "pure" in AX_MODES
        t = @elapsed p = critical_factor_optimize(ep, prof; faithful_confirm = false, extend_width = true, kw...)
        row[:pure] = p.sfmin; row[:pure_w] = p.width; row[:pure_ev] = p.evals; row[:pure_t] = t
        @printf("  pure-ext :ad      sfmin=%.5g  at (ky=%.4f, w=%.4f)  evals=%d  (AE-band, no confirm)  %.1fs\n",
                p.sfmin, p.kyhat, p.width, p.evals, t); flush(stdout)
    end
    if "robust" in AX_MODES
        t = @elapsed rob = critical_factor_robust(ep, prof; extend_width = true, kw...)
        row[:robust] = rob.sfmin; row[:robust_w] = rob.width
        row[:robust_ev] = rob.total_evals_full + rob.total_evals_eig; row[:robust_t] = t
        @printf("  robust_ad (ref)   sfmin=%.5g  at (ky=%.4f, w=%.4f)  evals=%d  %.1fs\n",
                rob.sfmin, rob.kyhat, rob.width, row[:robust_ev], t); flush(stdout)
    end
    return row
end

function main()
    # IR list: default = the 20-radius scan grid for this case.
    radii_env = get(ENV, "RADII", "")
    @printf("pure vs confirm width-extended :ad (%s)  NB=%d  modes=%s  threads=%d\n",
            USE_GPU ? "GPU" : "CPU", AX_NB, join(sort(collect(AX_MODES)), ","), Threads.nthreads()); flush(stdout)
    tglfep = joinpath(AX_CASE, "input_scan20_nb$(AX_NB).TGLFEP")
    opts, prof, _ = preprocess_gacode_inputs(AX_GACODE, tglfep)
    opts.N_BASIS = AX_NB

    # Mirror the production scan grid: IR_EXP[1:SCAN_N] (falling back to ir_exp_from_scan when unset).
    ir_exp = (!isempty(opts.IR_EXP) && !all(iszero, opts.IR_EXP)) ? Int.(opts.IR_EXP) :
             TJLFEP.ir_exp_from_scan(prof.NR, prof.IRS, opts.SCAN_N)
    radii = isempty(radii_env) ? Int.(ir_exp[1:opts.SCAN_N]) : parse.(Int, split(radii_env, ','))
    @printf("radii (%d): %s\n", length(radii), string(radii)); flush(stdout)

    # Warm up compilation on the first radius so per-radius timings exclude JIT.
    run_case(opts, prof, first(radii))

    rows = Dict{Symbol,Any}[]
    for ir in radii
        push!(rows, run_case(opts, prof, ir))
    end

    println("\n===================== SUMMARY (sfmin) =====================")
    @printf("  %-5s %-11s %-11s %-11s %-11s %-11s\n", "IR", "orig(w≥1)", "wide(1desc)", "confirm", "pure", "robust")
    for r in rows
        @printf("  %-5d %-11.5g %-11.5g %-11.5g %-11.5g %-11.5g\n",
                r[:ir], get(r, :orig, NaN), get(r, :wide, NaN), get(r, :confirm, NaN),
                get(r, :pure, NaN), get(r, :robust, NaN))
    end

    println("\n===================== SUMMARY (cost: per-radius wallclock s / evals) =====================")
    @printf("  %-5s %-9s %-9s %-9s %-9s %-9s | %-8s %-8s %-8s\n",
            "IR", "orig_t", "wide_t", "conf_t", "pure_t", "rob_t", "orig_ev", "wide_ev", "conf_ev")
    to = tw = tc = tp = tr = 0.0
    for r in rows
        ot = get(r, :orig_t, 0.0); wt = get(r, :wide_t, 0.0); ct = get(r, :confirm_t, 0.0)
        pt = get(r, :pure_t, 0.0); rt = get(r, :robust_t, 0.0)
        to += ot; tw += wt; tc += ct; tp += pt; tr += rt
        @printf("  %-5d %-9.2f %-9.2f %-9.2f %-9.2f %-9.2f | %-8d %-8d %-8d\n",
                r[:ir], ot, wt, ct, pt, rt, get(r, :orig_ev, 0), get(r, :wide_ev, 0), get(r, :confirm_ev, 0))
    end
    @printf("  %-5s %-9.2f %-9.2f %-9.2f %-9.2f %-9.2f\n", "SUM", to, tw, tc, tp, tr)
    to > 0 && @printf("  wide/orig = %.2fx   pure/orig = %.2fx   confirm/orig = %.2fx   (vs the fast single descent)\n",
                      tw / max(to, eps()), tp / max(to, eps()), tc / max(to, eps()))
end

main()
