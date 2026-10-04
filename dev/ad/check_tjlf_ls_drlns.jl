# Isolate the single-ky solver: does tjlf_LS propagate ∂γ/∂RLNS[is] correctly?
# Bypasses TJLF_map entirely. Builds the Float64 InputTJLF at the operating
# point, then compares AD (seed Dual on RLNS[is]) vs central FD (perturb
# RLNS[is], Float64 TJLF.run). A mismatch localizes the dropped partials to the
# matrix-assembly path inside tjlf_LS rather than the eigen rule.
#
#   module load julia/1.11.7 && export JULIA_DEPOT_PATH=$PSCRATCH/.julia
#   julia --project=. dev/ad/check_tjlf_ls_drlns.jl

using TJLFEP
using TJLF
using ForwardDiff
using Printf

const CASE = joinpath(@__DIR__, "..", "..", "examples", "DIIID_202017C42_500ms_v3.1")

function build_base(opts, prof)
    ep = deepcopy(opts)
    base = TJLFEP.TJLF_map(ep, prof)
    TJLFEP._configure_inputTJLF_for_ky!(base, ep)
    return ep, base
end

gamma_f64(inp) = TJLF.run(inp).eigenvalue[:, 1, 1]

function main()
    opts, prof, _ = preprocess_gacode_inputs(joinpath(CASE, "input.gacode"),
                                             joinpath(CASE, "input.TGLFEP"))
    opts.IR = 38; opts.KYHAT_IN = 0.25; opts.WIDTH_IN = 1.0; opts.FACTOR_IN = 2.5
    ep, base = build_base(opts, prof)
    is = ep.IS_EP + 1
    nm = ep.NMODES
    r0 = base.RLNS[is]
    @printf("RLNS[is=%d] = %.6f   (perturbing this directly)\n", is, r0)

    # ── AD: seed Dual on RLNS[is] (∂/∂RLNS[is] = 1) ──
    Tag = typeof(ForwardDiff.Tag(main, Float64))
    D   = ForwardDiff.Dual{Tag, Float64, 1}
    dual = TJLFEP._to_dual_input(base, D)
    RLNS = copy(dual.RLNS)
    RLNS[is] = ForwardDiff.Dual{Tag}(r0, 1.0)
    dual.RLNS = RLNS
    g = TJLFEP._tjlf_run_dual(dual, D).eigenvalue[:, 1, 1]
    ad = [ForwardDiff.partials(x, 1) for x in g]

    println("\n∂γ/∂RLNS[is]:  [mode:  AD          FD(h)        |Δ|        rel ]")
    for h in (1e-2, 1e-3, 1e-4)
        bp = deepcopy(base); bp.RLNS[is] = r0 + h
        bm = deepcopy(base); bm.RLNS[is] = r0 - h
        gp = gamma_f64(bp); gm = gamma_f64(bm)
        @printf("  h=%.0e\n", h)
        for n in 1:nm
            fd  = (gp[n] - gm[n]) / (2h)
            rel = abs(ad[n] - fd) / max(abs(fd), 1e-12)
            @printf("    mode %d:  %12.5e  %12.5e   %9.2e  %8.2e\n", n, ad[n], fd, abs(ad[n]-fd), rel)
        end
    end
end

main()
