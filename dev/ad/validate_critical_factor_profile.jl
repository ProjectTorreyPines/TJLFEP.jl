# Validate the AD radial driver (Phase 4): `critical_factor_profile` computes the
# AE-band critical factor sfmin(IR) over the SCAN_N radii with cross-radius
# CONTINUATION (each optimum warm-starts the next). We check that
#   (a) the continuation profile matches a COLD per-radius optimization
#       (seed grid at every radius), and
#   (b) continuation costs materially fewer eigensolves (only radius 1 pays the
#       seed grid; the rest converge from the neighbor's optimum).
#
#   module load julia/1.11.7 && export JULIA_DEPOT_PATH=$PSCRATCH/.julia
#   julia -t 8 --project=. dev/ad/validate_critical_factor_profile.jl

using TJLFEP
using Printf

const ITER = joinpath(@__DIR__, "..", "..", "examples", "ITER")

function iter_inputs()
    prof = TJLFEP.readMTGLF(joinpath(ITER, "input.MTGLF"))
    profile = prof[1]; ir_exp = prof[2]
    opts = TJLFEP.readTGLFEP(joinpath(ITER, "input.TGLFEP"), ir_exp)
    ni, Ti, dlnnidr, dlntidr, cs, rmin_ex, gammaE, gammap, omegaGAM =
        TJLFEP.read_expro_for_alpha(joinpath(ITER, "input.EXPRO"), profile, opts.IS_EP; gacode_file = nothing)
    expro = (ni=ni, Ti=Ti, dlnnidr=dlnnidr, dlntidr=dlntidr, cs=cs, rmin_ex=rmin_ex,
             gammaE=gammaE, gammap=gammap, omegaGAM=omegaGAM)
    TJLFEP._apply_runthd_expro_setup!(opts, profile, expro)
    return opts, profile
end

const NSEED = 3   # 3x3 seed grid (cheaper validation; production would use more)

function cold_profile(ep0, prof, radii; scan_lo, scan_hi)
    n = length(radii)
    sf = fill(Inf, n); kys = fill(NaN, n); ws = fill(NaN, n); tot = 0
    for (i, ir) in enumerate(radii)
        ep = deepcopy(ep0); ep.IR = ir
        r = critical_factor_optimize(ep, prof; seed = nothing, nseed_ky = NSEED, nseed_w = NSEED,
                                     scan_lo = scan_lo, scan_hi = scan_hi)
        sf[i] = r.sfmin; kys[i] = r.kyhat; ws[i] = r.width; tot += r.evals
    end
    return (; sfmin = sf, kyhat = kys, width = ws, total_evals = tot)
end

function main()
    @printf("threads = %d\n", Threads.nthreads())
    ep0, prof = iter_inputs()
    radii = [43, 83, 123]   # interior subset (continuation across neighbors)
    @printf("ITER N_BASIS=%d  SCAN_N=%d  validate radii=%s\n", ep0.N_BASIS, ep0.SCAN_N, radii)
    scan_lo = 1e-2; scan_hi = 35.0

    @printf("\n--- WARM continuation (critical_factor_profile) ---\n")
    warm = critical_factor_profile(ep0, prof; radii = radii, nseed_ky = NSEED, nseed_w = NSEED,
                                   scan_lo = scan_lo, scan_hi = scan_hi, verbose = false)
    for i in eachindex(radii)
        @printf("  IR=%3d  sfmin=%.5e  ky=%.4f  w=%.4f  evals=%4d  conv=%s\n",
                radii[i], warm.sfmin[i], warm.kyhat[i], warm.width[i], warm.evals[i], warm.converged[i])
    end
    @printf("  total eigensolves (warm) = %d\n", warm.total_evals)
    flush(stdout)

    @printf("\n--- COLD per-radius (seed grid each) ---\n")
    cold = cold_profile(ep0, prof, radii; scan_lo = scan_lo, scan_hi = scan_hi)
    for i in eachindex(radii)
        @printf("  IR=%3d  sfmin=%.5e  ky=%.4f  w=%.4f\n", radii[i], cold.sfmin[i], cold.kyhat[i], cold.width[i])
    end
    @printf("  total eigensolves (cold) = %d\n", cold.total_evals)

    @printf("\n--- comparison ---\n")
    for i in eachindex(radii)
        d = abs(warm.sfmin[i] - cold.sfmin[i])
        rel = d / max(abs(cold.sfmin[i]), 1e-12)
        @printf("  IR=%3d  |Δsfmin|=%.3e  rel=%.3e\n", radii[i], d, rel)
    end
    @printf("  continuation speedup (eigensolves) = %.2fx  (%d warm vs %d cold)\n",
            cold.total_evals / max(warm.total_evals, 1), warm.total_evals, cold.total_evals)
    flush(stdout)
end

main()
