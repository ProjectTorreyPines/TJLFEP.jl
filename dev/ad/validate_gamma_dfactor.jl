# Validate the single-point AD building block `gamma_dgamma_dfactor` against
# finite differences of the production single-ky path (`TJLFEP_ky`), on the
# DIII-D verification case.
#
# Three checks, in increasing strength:
#   (A) builder fidelity : the Float64 input built inside `gamma_dgamma_dfactor`
#       reproduces TJLFEP_ky's growth rate at the operating point (so the AD path
#       differentiates the SAME forward model the scan uses).
#   (B) AD value         : the Dual run's γ value matches the Float64 γ.
#   (C) AD derivative    : dγ/d(FACTOR_IN) from forward-mode AD matches a central
#       finite difference of TJLFEP_ky γ across FACTOR_IN ± h.
#
# Run (from the TJLFEP package dir):
#   module load julia/1.11.7 && export JULIA_DEPOT_PATH=$PSCRATCH/.julia
#   julia --project=. dev/ad/validate_gamma_dfactor.jl

using TJLFEP
using Printf

const CASE = joinpath(@__DIR__, "..", "..", "examples", "DIIID_202017C42_500ms_v3.1")

# ── Operating point (an unstable, well-separated point from out.scalefactor_r038)
const IR        = 38
const KYHAT_IN  = 0.25
const WIDTH_IN  = 1.0
const FACTOR_IN = 2.5

# Float64 growth rate at a given factor via the production single-ky path.
function gamma_f64(ep0, pr, factor)
    ep = deepcopy(ep0)
    ep.FACTOR_IN = factor
    gamma_out, freq_out, _, _, _, _, _ = TJLFEP.TJLFEP_ky(ep, pr, "", 0)
    nm = ep.NMODES
    return gamma_out[1:nm], freq_out[1:nm]
end

function main()
    gacode = joinpath(CASE, "input.gacode")
    tglfep = joinpath(CASE, "input.TGLFEP")
    opts, prof, _ = preprocess_gacode_inputs(gacode, tglfep)

    opts.IR        = IR
    opts.KYHAT_IN  = KYHAT_IN
    opts.WIDTH_IN  = WIDTH_IN
    opts.FACTOR_IN = FACTOR_IN

    @printf("DIII-D point: IR=%d  kyhat=%.4f  width=%.4f  factor=%.4f  SCAN_METHOD=%d  IS_EP=%d\n",
            IR, KYHAT_IN, WIDTH_IN, FACTOR_IN, opts.SCAN_METHOD, opts.IS_EP)

    # ── (A)+(B): AD pass returns γ and dγ/dfactor in one shot.
    res = gamma_dgamma_dfactor(opts, prof)
    nm  = opts.NMODES

    # Production Float64 reference at the same factor.
    g_ref, f_ref = gamma_f64(opts, prof, res.factor)

    println("\n(A) builder fidelity + (B) AD value  [mode: γ_TJLFEP_ky   γ_AD   |Δ|]")
    for n in 1:nm
        @printf("    mode %d:  %14.6e  %14.6e   %9.2e\n", n, g_ref[n], res.gamma[n], abs(g_ref[n]-res.gamma[n]))
    end

    # ── (C): central finite difference of TJLFEP_ky γ.
    println("\n(C) dγ/d(FACTOR_IN):  [mode:  AD        FD(h)       |Δ|       rel]")
    for h in (1e-2, 1e-3, 1e-4)
        gp, _ = gamma_f64(opts, prof, FACTOR_IN + h)
        gm, _ = gamma_f64(opts, prof, FACTOR_IN - h)
        @printf("  h=%.0e\n", h)
        for n in 1:nm
            fd  = (gp[n] - gm[n]) / (2h)
            ad  = res.dgamma_dfactor[n]
            rel = abs(ad - fd) / max(abs(fd), 1e-12)
            @printf("    mode %d:  %12.5e  %12.5e   %9.2e  %8.2e\n", n, ad, fd, abs(ad-fd), rel)
        end
    end

    println("\nγ      = ", res.gamma)
    println("dγ/df  = ", res.dgamma_dfactor)
    println("freq   = ", res.freq)
    println("dfreq/df = ", res.dfreq_dfactor)
end

main()
