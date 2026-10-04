# Overlay sfmin(IR) for the "widen-the-box" single-descent :ad vs the original w≥1-only :ad and
# robust_ad (nb=8, 20-radius DIII-D scan; data from job 54904500). Shows that a SINGLE fast descent
# over w∈[0.05, WIDTH_MAX] with denser seeds + 1 faithful confirm recovers the narrow-width modes the
# w≥1 solver over-predicts or misses entirely (orig=Inf at IR=95,101), landing within ~1-2x of
# robust_ad at ~7x lower cost than the full locate-grid extension — and conservatively (≥ robust),
# unlike the no-confirm pure-ext which collapsed to the floor.
#   julia --project=. dev/ad/plot_ad_wide_box.jl
using Pkg
Pkg.activate(normpath(@__DIR__, "..", ".."))
using Plots, Printf

const OUT_PNG = get(ENV, "OUT_PNG", joinpath(@__DIR__, "sfmin_ad_wide_box_nb8.png"))
const FLOOR   = 10.0 / 512.0
const CEIL    = 10.0   # FACTOR_IN; orig railed here / returned Inf at the steep edge

ir     = [2, 7, 12, 17, 22, 28, 33, 38, 43, 48, 54, 59, 64, 69, 74, 80, 85, 90, 95, 101]
orig   = [0.70016, 0.91687, 1.8085, 0.15086, 0.081934, 0.080392, 1.5583, 1.949, 0.034778, 0.11987,
          0.16746, 0.18921, 0.12775, 1.443, 0.49991, 0.79179, 1.025, 2.7941, NaN, NaN]  # NaN = Inf (no mode)
wide   = [0.92891, 0.48951, 0.19068, 0.16629, 0.083654, 0.086003, 0.021234, 0.029437, 0.026629,
          0.02044, 0.083958, 0.051013, 0.043718, 0.21354, 0.31428, 0.37598, 0.26996, 0.51768,
          0.24878, 0.203]
robust = [0.77859, 0.47823, 0.19129, 0.16649, 0.06296, 0.067384, 0.019531, 0.019531, 0.019531,
          0.019531, 0.032775, 0.05211, 0.027744, 0.12406, 0.20779, 0.24721, 0.23323, 0.24421,
          0.22529, 0.11341]

function main()
    default(legendfontsize=10, guidefontsize=11, tickfontsize=9, dpi=200, fontfamily="Computer Modern")

    plt = plot(ir, robust; marker=:star5, lw=2.6, ms=7, label="robust_ad (faithful ref)",
               color=:firebrick, yscale=:log10, ylims=(0.015, 12),
               xlabel="radial index IR", ylabel="critical scale factor  sfmin",
               title="DIII-D  N_BASIS=8  SCAN_N=20 :  widen-the-box single-descent :ad",
               legend=:bottomright, size=(950, 540))
    plot!(plt, ir, orig; marker=:utriangle, lw=2, ms=6, linestyle=:dashdot,
          label="orig :ad (w≥1 only, single descent)", color=:steelblue)
    plot!(plt, ir, wide; marker=:diamond, lw=2.2, ms=6,
          label="wide :ad (w∈[0.05,2] single descent + 1 confirm)", color=:darkorange)
    hline!(plt, [FLOOR]; lw=1, linestyle=:dot, color=:gray, label="scan floor")

    # mark radii where orig found no mode at all (Inf), at the plot ceiling
    inf_ir = ir[isnan.(orig)]
    isempty(inf_ir) || scatter!(plt, inf_ir, fill(CEIL, length(inf_ir)); marker=:xcross, ms=8,
                                color=:steelblue, label="orig = Inf (no w≥1 mode)")

    mkpath(dirname(OUT_PNG))
    savefig(plt, OUT_PNG)
    println("wrote ", OUT_PNG)
end

main()
