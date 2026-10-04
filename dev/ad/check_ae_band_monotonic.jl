# Diagnostic: is the AE-band-filtered leading growth rate a clean (monotonic,
# no false floor) instability indicator vs the raw max-over-modes γ?
#
# For a few (kyhat,width) grid points it scans FACTOR_IN over a geometric grid
# and prints, at each factor:
#   - raw γ_lead = max over ALL modes,  and the frequency of that argmax mode
#   - AE  γ_keep = max over modes with freq < FREQ_AE_UPPER (the primary keep
#                  filter in TJLFEP_ky), and which mode index it is
#   - per-mode (γ, freq) so the AE/background split is visible
#
# What we want to confirm:
#   (a) raw γ_lead is contaminated by background (non-AE-frequency) modes that
#       are unstable even at tiny factor → false "unstable at floor";
#   (b) AE γ_keep is ~0 below the EP-AE onset, rises through it, and stays above
#       threshold up to scan_hi → an endpoint (scan_hi) test classifies the point
#       correctly, enabling the 1-eval stable bail.
#
#   module load julia/1.11.7 && export JULIA_DEPOT_PATH=$PSCRATCH/.julia
#   julia -t 8 --project=. dev/ad/check_ae_band_monotonic.jl

using TJLFEP
using Printf

const CASE = joinpath(@__DIR__, "..", "..", "examples", "DIIID_202017C42_500ms_v3.1")
const IR = 38

# (kyhat, width) probes: the key low-onset point + a couple of neighbours.
const PROBES = [(0.25, 1.0), (0.50, 1.0), (0.75, 1.0), (0.25, 1.571)]

function scan_point(opts0, prof, fu, kyhat, width)
    @printf("\n===== kyhat=%.3f width=%.3f   (FREQ_AE_UPPER=%.4e) =====\n", kyhat, width, fu)
    println("   factor      raw_γ   raw_mode(freq)     AE_γ_keep  AE_mode   per-mode (γ | freq)")
    factors = exp.(range(log(1.0e-3), log(10.0); length = 12))
    prev_keep = -Inf
    monotonic = true
    for f in factors
        ep = deepcopy(opts0)
        ep.IR = IR; ep.KYHAT_IN = kyhat; ep.WIDTH_IN = width; ep.FACTOR_IN = f
        r = gamma_dgamma_dfactor(ep, prof)
        iraw = argmax(r.gamma)
        cand = findall(<(fu), r.freq)
        if isempty(cand)
            keepγ = 0.0; jkeep = 0
        else
            jkeep = cand[argmax(@view r.gamma[cand])]
            keepγ = r.gamma[jkeep]
        end
        permode = join([@sprintf("%d:%.3e|%+.3e", n, r.gamma[n], r.freq[n]) for n in eachindex(r.gamma)], "  ")
        @printf("  %9.3e  %+.3e  %d(%+.3e)   %+.3e  %s   %s\n",
                f, r.gamma[iraw], iraw, r.freq[iraw], keepγ, jkeep == 0 ? "-" : string(jkeep), permode)
        if keepγ + 1e-12 < prev_keep
            monotonic = false
        end
        prev_keep = max(prev_keep, keepγ)
    end
    @printf("  AE γ_keep monotonic non-decreasing? %s\n", monotonic ? "YES" : "NO")
end

function main()
    @printf("threads = %d\n", Threads.nthreads())
    opts0, prof, _ = preprocess_gacode_inputs(joinpath(CASE, "input.gacode"),
                                              joinpath(CASE, "input.TGLFEP"))
    opts0.IR = IR
    opts0.N_BASIS = 6
    @printf("N_BASIS = %d\n", opts0.N_BASIS)
    fu = TJLFEP._ae_band_upper(opts0, prof)
    for (ky, w) in PROBES
        scan_point(opts0, prof, fu, ky, w)
    end
end

main()
