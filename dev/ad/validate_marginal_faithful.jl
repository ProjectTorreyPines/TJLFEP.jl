# Validate the AD-accelerated FAITHFUL marginal-factor finder (Phase 2):
# `marginal_factor_faithful` must reproduce the production keep-onset (smallest
# factor where any mode is KEPT under the full IFLUX=true filters) on BOTH:
#   - DIII-D : band-entry / growth-threshold binding (binding=:ae_band_growth)
#   - ITER   : a secondary filter (th-pinch / QL-ratio) binds, pushing the onset
#              well above the frequency band-entry (binding ∈ pinch/QL flags)
#
# Reference: a brute fine faithful factor sweep (full TJLFEP_ky per sample) whose
# lower kept-window edge is refined by bisection. We compare the AD value to that
# ground truth and report the eigensolve/full-eval economy.
#
#   module load julia/1.11.7 && export JULIA_DEPOT_PATH=$PSCRATCH/.julia
#   julia -t 8 --project=. dev/ad/validate_marginal_faithful.jl

using TJLFEP
using Printf

const DIIID = joinpath(@__DIR__, "..", "..", "examples", "DIIID_202017C42_500ms_v3.1")
const ITER  = joinpath(@__DIR__, "..", "..", "examples", "ITER")

# Fair, well-resolved reference: brute LINEAR sweep over a given [lo,hi] window
# (faithful IFLUX=true) returning the lowest kept factor, bisected to xtol.
function brute_onset(ep0, prof; kyhat, width, lo, hi, npts, xtol = 1e-4)
    nfull = Ref(0)
    keep(f) = (nfull[] += 1; TJLFEP.keep_at(ep0, prof, f; kyhat = kyhat, width = width).any_keep)
    fs = range(lo, hi; length = npts)
    a = NaN; b = NaN; preva = lo
    for f in fs
        if keep(f)
            b = f; a = preva
            break
        end
        preva = f
    end
    isnan(b) && return (NaN, nfull[])
    while (b - a) > xtol
        m = 0.5 * (a + b)
        keep(m) ? (b = m) : (a = m)
    end
    return (b, nfull[])
end

function run_case(label, ep0, prof; kyhat, width, lo, hi)
    @printf("\n========== %s  (IR=%d kyhat=%.3f width=%.3f) ==========\n", label, ep0.IR, kyhat, width)

    r = marginal_factor_faithful(ep0, prof; kyhat = kyhat, width = width,
                                 scan_lo = lo, scan_hi = hi, verbose = true)
    @printf("  AD faithful onset    = %.5e   binding=%s   pinned_lo=%s\n", r.factor_faithful, r.binding, r.pinned_lo)
    @printf("    AE unstable hull   = (%.5e, %.5e)   fast AD onset f1 = %.5e\n", r.window[1], r.window[2], r.factor_fast)
    @printf("    eigensolves(IFLUX=false) = %d ;  full evals(IFLUX=true) = %d ;  total = %d\n",
            r.evals_eig, r.evals_full, r.evals_eig + r.evals_full)
    @printf("    kept_modes=%s\n", r.kept_modes)

    # Independent reference: fine LINEAR brute sweep over the SAME physical hull.
    if r.converged
        f1, f2 = r.window
        gt, gt_full = brute_onset(ep0, prof; kyhat = kyhat, width = width, lo = f1, hi = f2, npts = 120)
        rel = abs(r.factor_faithful - gt) / max(abs(gt), 1e-12)
        @printf("  brute ref (120-pt linear in hull) onset = %.5e  (%d full evals)\n", gt, gt_full)
        @printf("    |Δ AD vs brute| = %.3e   rel = %.3e\n", abs(r.factor_faithful - gt), rel)
    end
    flush(stdout)
    return r
end

function diiid_inputs()
    opts, prof, _ = preprocess_gacode_inputs(joinpath(DIIID, "input.gacode"), joinpath(DIIID, "input.TGLFEP"))
    opts.IR = 38
    return opts, prof
end

function iter_inputs()
    prof = TJLFEP.readMTGLF(joinpath(ITER, "input.MTGLF"))
    profile = prof[1]; ir_exp = prof[2]
    opts = TJLFEP.readTGLFEP(joinpath(ITER, "input.TGLFEP"), ir_exp)
    ni, Ti, dlnnidr, dlntidr, cs, rmin_ex, gammaE, gammap, omegaGAM =
        TJLFEP.read_expro_for_alpha(joinpath(ITER, "input.EXPRO"), profile, opts.IS_EP; gacode_file = nothing)
    expro = (ni=ni, Ti=Ti, dlnnidr=dlnnidr, dlntidr=dlntidr, cs=cs, rmin_ex=rmin_ex,
             gammaE=gammaE, gammap=gammap, omegaGAM=omegaGAM)
    TJLFEP._apply_runthd_expro_setup!(opts, profile, expro)
    opts.IR = opts.IR_EXP[3]
    return opts, profile
end

function main()
    @printf("threads = %d\n", Threads.nthreads())

    dep, dpr = diiid_inputs()
    @printf("DIII-D N_BASIS=%d\n", dep.N_BASIS)
    run_case("DIII-D band-entry", dep, dpr; kyhat = 0.25, width = 1.571, lo = 1e-4, hi = 0.3)

    iep, ipr = iter_inputs()
    @printf("ITER N_BASIS=%d\n", iep.N_BASIS)
    run_case("ITER secondary-filter", iep, ipr; kyhat = 0.151, width = 1.494, lo = 0.5, hi = 35.0)
end

main()
