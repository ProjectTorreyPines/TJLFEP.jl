# Validate the AD-exact marginal_factor on the ITER case (ROTATIONAL_SUPPRESSION_FLAG=1,
# finite GAMMA_THRESH), where the marginal is a smooth grid-convergent root.
#
# For the two genuinely-unstable radii (i=3 IR=83, i=4 IR=123), at the operating
# point (kyhat,width) the converged brute scan settled on, we compare:
#   - AD-Newton marginal_factor(ae_band=true, auto rotational γ*)   [few eigensolves]
#   - a dense factor scan of AE-band γ_keep, root-found to γ*       [reference]
# They use the same γ(factor), so agreement validates the Newton root-finder; the
# AD value reaching the converged sfmin (~20, ~112) validates fidelity to truth.
#
#   module load julia/1.11.7 && export JULIA_DEPOT_PATH=$PSCRATCH/.julia
#   julia -t 8 --project=. dev/ad/validate_iter_marginal.jl

using TJLFEP
using Printf

const DIR = joinpath(@__DIR__, "..", "..", "examples", "ITER")

function build_inputs()
    prof = TJLFEP.readMTGLF(joinpath(DIR, "input.MTGLF"))
    profile = prof[1]; ir_exp = prof[2]
    Options = TJLFEP.readTGLFEP(joinpath(DIR, "input.TGLFEP"), ir_exp)
    ni, Ti, dlnnidr, dlntidr, cs, rmin_ex, gammaE, gammap, omegaGAM =
        TJLFEP.read_expro_for_alpha(joinpath(DIR, "input.EXPRO"), profile, Options.IS_EP; gacode_file = nothing)
    expro = (ni=ni, Ti=Ti, dlnnidr=dlnnidr, dlntidr=dlntidr, cs=cs, rmin_ex=rmin_ex,
             gammaE=gammaE, gammap=gammap, omegaGAM=omegaGAM)
    TJLFEP._apply_runthd_expro_setup!(Options, profile, expro)
    return Options, profile
end

# dense reference: scan AE-band kept γ over factor, linear-interp first up-crossing of γ*
function ref_marginal(base, prof, ky, w, gstar; fmax = 200.0, n = 60)
    fu = -abs(prof.omegaGAM[base.IR])
    fs = exp.(range(log(1.0e-2), log(fmax); length = n))
    prevf = NaN; prevg = NaN; fcross = NaN
    for f in fs
        ep = deepcopy(base); ep.KYHAT_IN = ky; ep.WIDTH_IN = w; ep.FACTOR_IN = f
        r = gamma_dgamma_dfactor(ep, prof)
        cand = findall(<(fu), r.freq)
        gk = isempty(cand) ? 0.0 : maximum(@view r.gamma[cand])
        if isnan(fcross) && !isnan(prevg) && prevg < gstar && gk >= gstar
            fcross = prevf + (gstar - prevg) * (f - prevf) / (gk - prevg)
        end
        prevf = f; prevg = gk
    end
    return fcross
end

function main()
    @printf("threads = %d\n", Threads.nthreads())
    opts, prof = build_inputs()

    # (scan index i, IR, converged kymark, converged width, converged-ref sfmin)
    cases = [(3, opts.IR_EXP[3], 0.151, 1.494, 20.1),
             (4, opts.IR_EXP[4], 0.156, 1.250, 112.6)]

    for (i, ir, ky, w, sref) in cases
        ep = deepcopy(opts); ep.IR = ir; ep.KYHAT_IN = ky; ep.WIDTH_IN = w
        gstar = TJLFEP._gamma_thresh_for(ep, prof)
        @printf("\n===== ITER i=%d IR=%d  kyhat=%.3f width=%.3f =====\n", i, ir, ky, w)
        @printf("  rotational γ* = %.6e   (converged brute sfmin ≈ %.1f)\n", gstar, sref)

        rm = marginal_factor(ep, prof; scan_lo = 1.0e-2, scan_hi = 200.0,
                             nscan = 16, ae_band = true)
        @printf("  AD-Newton   f★ = %.5e   (evals=%d, iters=%d, converged=%s)\n",
                rm.factor, rm.evals, rm.iters, rm.converged)
        flush(stdout)
        fref = ref_marginal(ep, prof, ky, w, gstar)
        @printf("  dense-scan  f★ = %.5e\n", fref)
        if !isnan(fref) && fref != 0
            @printf("  AD/dense ratio = %.4f\n", rm.factor / fref)
        end
    end
end

main()
