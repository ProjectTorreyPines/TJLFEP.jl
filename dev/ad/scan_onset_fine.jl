# Characterize the EP-AE onset on a FINE factor grid, to (a) see the band-entry
# discontinuity shape and (b) test whether the AD-Newton f*=0.0114 at
# (ky=0.25,width=1.571) is a root-finding artifact rather than real physics.
#
# At each factor we run gamma_dgamma_dfactor (IFLUX=false: γ, freq, dγ/dfactor per
# mode) and report the AE-band kept leading mode (freq < FREQ_AE_UPPER):
#   γ_keep, its mode index, dγ_keep/dfactor, and the mode's frequency.
#
# Then three onset estimates of f* (where kept-AE γ first reaches GAMMA_THRESH):
#   - fine first-crossing (linear interp between the two fine samples straddling
#     the band-entry)        → the grid-converged reference
#   - AD-Newton marginal_factor(ae_band=true)            → what the benchmark used
# A large gap means AD lands below the true band-entry onset (the artifact).
#
#   module load julia/1.11.7 && export JULIA_DEPOT_PATH=$PSCRATCH/.julia
#   julia -t 8 --project=. dev/ad/scan_onset_fine.jl

using TJLFEP
using Printf

const CASE = joinpath(@__DIR__, "..", "..", "examples", "DIIID_202017C42_500ms_v3.1")
const IR = 38
const THRESH = 1.0e-7

function keepγ(r, fu)
    cand = findall(<(fu), r.freq)
    isempty(cand) && return (0.0, 0, 0.0, NaN)
    j = cand[argmax(@view r.gamma[cand])]
    return (r.gamma[j], j, r.dgamma_dfactor[j], r.freq[j])
end

function scan(opts0, prof, fu, ky, w)
    @printf("\n===== kyhat=%.3f width=%.3f  (FREQ_AE_UPPER=%.4e, thresh=%.1e) =====\n", ky, w, fu, THRESH)
    println("   factor       AE_γ_keep   mode  dγ_keep/df    mode_freq")
    fs = exp.(range(log(4.0e-3), log(6.0e-2); length = 40))
    prevf = NaN; prevg = NaN; fcross = NaN
    for f in fs
        ep = deepcopy(opts0); ep.IR = IR; ep.KYHAT_IN = ky; ep.WIDTH_IN = w; ep.FACTOR_IN = f
        r = gamma_dgamma_dfactor(ep, prof)
        gk, j, dgk, frq = keepγ(r, fu)
        @printf("  %10.4e   %+.4e   %s   %+.4e   %s\n",
                f, gk, j == 0 ? "-" : string(j), dgk, isnan(frq) ? "-" : @sprintf("%+.3e", frq))
        if isnan(fcross) && !isnan(prevg) && prevg < THRESH && gk >= THRESH
            # linear interp in factor to γ=THRESH across the band-entry
            fcross = prevf + (THRESH - prevg) * (f - prevf) / (gk - prevg)
        end
        prevf = f; prevg = gk
    end

    ep = deepcopy(opts0); ep.IR = IR; ep.KYHAT_IN = ky; ep.WIDTH_IN = w
    rm = marginal_factor(ep, prof; gamma_thresh = THRESH, scan_lo = 1.0e-3, scan_hi = 10.0, ae_band = true)
    @printf("  fine first-crossing f* = %.5e\n", fcross)
    @printf("  AD-Newton          f* = %.5e   (evals=%d, converged=%s)\n", rm.factor, rm.evals, rm.converged)
end

function main()
    @printf("threads = %d\n", Threads.nthreads())
    opts0, prof, _ = preprocess_gacode_inputs(joinpath(CASE, "input.gacode"),
                                              joinpath(CASE, "input.TGLFEP"))
    opts0.IR = IR; opts0.N_BASIS = 6
    @printf("N_BASIS = %d\n", opts0.N_BASIS)
    fu = -abs(prof.omegaGAM[IR])
    scan(opts0, prof, fu, 0.25, 1.571)
    scan(opts0, prof, fu, 0.297, 1.571)
end

main()
