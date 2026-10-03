# Validate the IFT-gradient (kyhat,width) optimizer (Phase 3):
# `critical_factor_optimize` must find the same min of the AE-band marginal factor
# f★(ky,w) as a brute (ky,w) grid scan of `marginal_factor`, at far fewer
# eigensolves. The IFT gradient ∂f★/∂(ky,w) = -(∂γ/∂(ky,w))/(∂γ/∂f) comes free
# from one gamma_grad pass at each marginal point.
#
#   module load julia/1.11.7 && export JULIA_DEPOT_PATH=$PSCRATCH/.julia
#   julia -t 8 --project=. dev/ad/validate_critical_factor_optimize.jl

using TJLFEP
using Printf

const DIIID = joinpath(@__DIR__, "..", "..", "examples", "DIIID_202017C42_500ms_v3.1")
const ITER  = joinpath(@__DIR__, "..", "..", "examples", "ITER")

# Brute (ky,w) grid min of the AE-band marginal factor (same objective the
# optimizer descends). Returns (fmin, ky*, w*, total_eigensolves).
function brute_grid(ep0, prof; gth, nky, nw, scan_lo, scan_hi)
    wlo = Float64(ep0.WIDTH_MIN); whi = Float64(ep0.WIDTH_MAX)
    kys = [(1.0 / nky) * i for i in 1:nky]
    ws  = [wlo + (whi - wlo) * (i - 1) / (nw - 1) for i in 1:nw]
    fmin = Inf; kymin = NaN; wmin = NaN; tot = 0
    for ky in kys, w in ws
        ep = deepcopy(ep0); ep.KYHAT_IN = ky; ep.WIDTH_IN = w
        mf = marginal_factor(ep, prof; gamma_thresh = gth, ae_band = true,
                             scan_lo = scan_lo, scan_hi = scan_hi)
        tot += mf.evals
        if mf.converged && mf.factor < fmin
            fmin = mf.factor; kymin = ky; wmin = w
        end
    end
    return (fmin, kymin, wmin, tot)
end

function run_case(label, ep0, prof; scan_lo, scan_hi, nky = 8, nw = 8)
    @printf("\n========== %s  (IR=%d) ==========\n", label, ep0.IR)
    gth = TJLFEP._gamma_thresh_for(ep0, prof)

    fb, kyb, wb, ngrid = brute_grid(ep0, prof; gth = gth, nky = nky, nw = nw, scan_lo = scan_lo, scan_hi = scan_hi)
    @printf("  brute %dx%d grid:   f★min = %.5e at (ky=%.4f, w=%.4f)   [%d eigensolves]\n",
            nky, nw, fb, kyb, wb, ngrid)

    r = critical_factor_optimize(ep0, prof; scan_lo = scan_lo, scan_hi = scan_hi, verbose = false)
    @printf("  IFT optimizer:     sfmin = %.5e at (ky=%.4f, w=%.4f)   [%d eigensolves, %d iters, conv=%s]\n",
            r.sfmin, r.kyhat, r.width, r.evals, r.iters, r.converged)
    @printf("    seed-grid min    = %.5e at (ky=%.4f, w=%.4f)\n", r.f_seedmin, r.ky_seed, r.w_seed)
    rel = abs(r.sfmin - fb) / max(abs(fb), 1e-12)
    @printf("    |Δ opt vs grid|  = %.3e   rel = %.3e   speedup(eigensolves) = %.1fx\n",
            abs(r.sfmin - fb), rel, ngrid / max(r.evals, 1))
    flush(stdout)
    return r
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
    iep, ipr = iter_inputs()
    @printf("ITER N_BASIS=%d\n", iep.N_BASIS)
    run_case("ITER AE-band onset min", iep, ipr; scan_lo = 1e-2, scan_hi = 35.0, nky = 8, nw = 8)
end

main()
