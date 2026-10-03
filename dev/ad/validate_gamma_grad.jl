# Validate the multi-variable AD building block `gamma_grad` (Phase 1) on the
# DIII-D verification case. `gamma_grad` returns per-mode γ/freq and their exact
# partials w.r.t. (FACTOR_IN, KYHAT_IN, WIDTH_IN) in ONE forward-mode AD pass.
#
# Checks:
#   (A) value      : γ/freq values match the single-seed gamma_dgamma_dfactor.
#   (B) factor col : ∂γ/∂FACTOR_IN matches gamma_dgamma_dfactor exactly (same seed).
#   (C) kyhat col  : ∂γ/∂KYHAT_IN matches a central finite difference of γ.
#   (D) width col  : ∂γ/∂WIDTH_IN  matches a central finite difference of γ.
#
# Run (from the TJLFEP package dir):
#   module load julia/1.11.7 && export JULIA_DEPOT_PATH=$PSCRATCH/.julia
#   julia --project=. dev/ad/validate_gamma_grad.jl

using TJLFEP
using Printf

const CASE = joinpath(@__DIR__, "..", "..", "examples", "DIIID_202017C42_500ms_v3.1")

# Interior, well-separated unstable operating point (same as gamma_dgamma_dfactor validation).
const IR        = 38
const KYHAT_IN  = 0.25
const WIDTH_IN  = 1.0
const FACTOR_IN = 2.5

# Per-mode γ at a perturbed (kyhat, width, factor) via the AD building block's value.
function gamma_at(ep0, pr; kyhat = KYHAT_IN, width = WIDTH_IN, factor = FACTOR_IN)
    ep = deepcopy(ep0)
    ep.KYHAT_IN  = kyhat
    ep.WIDTH_IN  = width
    ep.FACTOR_IN = factor
    r = gamma_dgamma_dfactor(ep, pr)
    return r.gamma
end

function main()
    gacode = joinpath(CASE, "input.gacode")
    tglfep = joinpath(CASE, "input.TGLFEP")
    opts, prof, _ = preprocess_gacode_inputs(gacode, tglfep)

    opts.IR        = IR
    opts.KYHAT_IN  = KYHAT_IN
    opts.WIDTH_IN  = WIDTH_IN
    opts.FACTOR_IN = FACTOR_IN

    @printf("DIII-D point: IR=%d  kyhat=%.4f  width=%.4f  factor=%.4f\n",
            IR, KYHAT_IN, WIDTH_IN, FACTOR_IN)

    vars = (:FACTOR_IN, :KYHAT_IN, :WIDTH_IN)
    g = gamma_grad(opts, prof; vars = vars)
    nm = length(g.gamma)
    col = Dict(v => k for (k, v) in enumerate(vars))

    # Single-seed reference for the FACTOR_IN column + value cross-check.
    ref = gamma_dgamma_dfactor(opts, prof)

    println("\n(A) value match + (B) ∂γ/∂FACTOR_IN vs gamma_dgamma_dfactor")
    println("  [mode:  γ_grad        γ_ref       |Δγ|   |   ∂/∂f grad     ∂/∂f ref     |Δ|]")
    for n in 1:nm
        @printf("   %d:  %12.5e  %12.5e  %8.1e | %12.5e  %12.5e  %8.1e\n",
                n, g.gamma[n], ref.gamma[n], abs(g.gamma[n]-ref.gamma[n]),
                g.dgamma[n, col[:FACTOR_IN]], ref.dgamma_dfactor[n],
                abs(g.dgamma[n, col[:FACTOR_IN]] - ref.dgamma_dfactor[n]))
    end

    # ── (C),(D): finite-difference the kyhat / width columns.
    for (label, var, base, hs) in (
            ("KYHAT_IN", :KYHAT_IN, KYHAT_IN, (1e-3, 1e-4)),
            ("WIDTH_IN",  :WIDTH_IN,  WIDTH_IN,  (1e-3, 1e-4)))
        println("\n(", label == "KYHAT_IN" ? "C" : "D", ") ∂γ/∂", label, "  [mode:  AD          FD(h)        |Δ|        rel]")
        for h in hs
            gp = label == "KYHAT_IN" ? gamma_at(opts, prof; kyhat = base + h) : gamma_at(opts, prof; width = base + h)
            gm = label == "KYHAT_IN" ? gamma_at(opts, prof; kyhat = base - h) : gamma_at(opts, prof; width = base - h)
            @printf("  h=%.0e\n", h)
            for n in 1:nm
                fd  = (gp[n] - gm[n]) / (2h)
                ad  = g.dgamma[n, col[var]]
                rel = abs(ad - fd) / max(abs(fd), 1e-12)
                @printf("    mode %d:  %12.5e  %12.5e   %9.2e  %8.2e\n", n, ad, fd, abs(ad-fd), rel)
            end
        end
    end

    println("\nγ     = ", g.gamma)
    println("∂γ/∂(factor,kyhat,width):")
    for n in 1:nm
        @printf("  mode %d: % .5e  % .5e  % .5e\n", n, g.dgamma[n,1], g.dgamma[n,2], g.dgamma[n,3])
    end
end

main()
