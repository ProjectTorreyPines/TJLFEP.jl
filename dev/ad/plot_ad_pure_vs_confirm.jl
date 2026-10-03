# Overlay sfmin(IR) for the width-extended :ad variants vs robust_ad, nb=8, 20-radius DIII-D scan
# (data from the pure-vs-confirm GPU profile, job 54847106). Shows that the "pure AD" no-confirm
# variant (faithful_confirm=false) collapses to spuriously low / floor-pinned sfmin at the
# narrow-width edge (IR≈69–95), where the AE-band-unstable modes are rejected by the keep filters
# that the faithful confirm applies. confirm-ext (production) matches robust_ad bitwise there.
#   julia --project=. dev/ad/plot_ad_pure_vs_confirm.jl
using Pkg
Pkg.activate(normpath(@__DIR__, "..", ".."))
using Plots, Printf

# Experimental artifact (not a docs plot): write alongside the profile scripts.
const OUT_PNG = get(ENV, "OUT_PNG", joinpath(@__DIR__, "sfmin_ad_pure_vs_confirm_nb8.png"))
const FLOOR   = 10.0 / 512.0   # scan_lo = FACTOR_IN/512

ir      = [2, 7, 12, 17, 22, 28, 33, 38, 43, 48, 54, 59, 64, 69, 74, 80, 85, 90, 95, 101]
confirm = [0.881, 0.37194, 0.19165, 0.16687, 0.09135, 0.091357, 0.019531, 0.019531, 0.019531,
           0.019531, 0.051031, 0.05211, 0.027744, 0.12406, 0.20779, 0.24721, 0.23323, 0.24421,
           0.22529, 0.14607]
pure    = [0.5883, 0.29681, 0.18916, 0.16133, 0.081876, 0.091218, 0.019531, 0.019531, 0.019531,
           0.019531, 0.019531, 0.019531, 0.019531, 0.024205, 0.027036, 0.019531, 0.019531,
           0.019531, 0.030711, 0.11633]
robust  = [0.77859, 0.47823, 0.19129, 0.16649, 0.06296, 0.067384, 0.019531, 0.019531, 0.019531,
           0.019531, 0.032775, 0.05211, 0.027744, 0.12406, 0.20779, 0.24721, 0.23323, 0.24421,
           0.22529, 0.11341]

function main()
    default(legendfontsize=10, guidefontsize=11, tickfontsize=9, dpi=200, fontfamily="Computer Modern")

    plt = plot(ir, robust; marker=:star5, lw=2.6, ms=7, label="robust_ad (faithful ref)",
               color=:firebrick, yscale=:log10,
               xlabel="radial index IR", ylabel="critical scale factor  sfmin",
               title="DIII-D  N_BASIS=8  SCAN_N=20 :  width-extended :ad — pure vs confirm",
               legend=:bottomright, size=(900, 520))
    plot!(plt, ir, confirm; marker=:circle, lw=2, ms=5,
          label="confirm-ext :ad (faithful_confirm=true, production)", color=:seagreen)
    plot!(plt, ir, pure; marker=:diamond, lw=2, ms=5, linestyle=:dash,
          label="pure-ext :ad (faithful_confirm=false)", color=:darkorange)
    hline!(plt, [FLOOR]; lw=1, linestyle=:dot, color=:gray, label="scan floor = FACTOR_IN/512")

    # shade the narrow-width edge where pure collapses
    vspan!(plt, [66, 97]; color=:orange, alpha=0.08, label="")

    mkpath(dirname(OUT_PNG))
    savefig(plt, OUT_PNG)
    println("wrote ", OUT_PNG)
end

main()
