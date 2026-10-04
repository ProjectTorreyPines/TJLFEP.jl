# Demo + validation of the AD-accelerated Newton marginal-factor finder.
#
#   1. Scans γ_lead(factor) on a coarse log grid (to show the instability curve
#      and the threshold crossing).
#   2. Runs `marginal_factor` (safeguarded Newton with exact dγ/dfactor).
#   3. Reports the iteration trace, the residual γ_lead(f★)-thresh, and the
#      eigensolve count — the speed-up vs the 8-point-per-round Fortran scan.
#
#   module load julia/1.11.7 && export JULIA_DEPOT_PATH=$PSCRATCH/.julia
#   julia --project=. dev/ad/demo_marginal_newton.jl

using TJLFEP
using Printf

const CASE = joinpath(@__DIR__, "..", "..", "examples", "DIIID_202017C42_500ms_v3.1")
const IR, KYHAT_IN, WIDTH_IN = 38, 0.25, 1.0
const GTHRESH = 1.0e-7

function main()
    opts, prof, _ = preprocess_gacode_inputs(joinpath(CASE, "input.gacode"),
                                             joinpath(CASE, "input.TGLFEP"))
    opts.IR = IR; opts.KYHAT_IN = KYHAT_IN; opts.WIDTH_IN = WIDTH_IN

    @printf("Marginal-factor Newton  |  IR=%d kyhat=%.3f width=%.3f  thresh=%.1e\n",
            IR, KYHAT_IN, WIDTH_IN, GTHRESH)

    # ── 1) coarse γ_lead(factor) curve ──
    println("\nγ_lead(factor) scan:")
    fs = exp.(range(log(1e-3), log(5.0); length = 8))
    for f in fs
        gl, dgl, _ = TJLFEP._gamma_lead_dfactor(opts, prof, f)
        @printf("    f=%9.4e   γ_lead=%12.5e   dγ/df=%12.5e   %s\n",
                f, gl, dgl, gl > GTHRESH ? "UNSTABLE" : "stable")
    end

    # ── 2) Newton root-find ──
    println("\nsafeguarded Newton:")
    res = marginal_factor(opts, prof; gamma_thresh = GTHRESH,
                          scan_lo = 1e-3, scan_hi = 5.0, nscan = 8, verbose = true)

    @printf("\nRESULT: f★ = %.8f   γ_lead(f★) = %.6e   residual = %.2e\n",
            res.factor, res.gamma_lead, res.gamma_lead - GTHRESH)
    @printf("        Newton iters = %d   total eigensolves = %d   converged = %s\n",
            res.iters, res.evals, res.converged)
    println("        final bracket = ", res.bracket)
    println("        γ modes at f★ = ", res.gamma)
end

main()
