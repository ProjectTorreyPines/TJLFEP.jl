# Overlay critical scale factor sfmin(radius) for the three solver tiers at DIII-D N_BASIS=32,
# SCAN_N=20: :grid (Fortran-equivalent kwscale_scan, w≥1 only), :robust_ad (production: w≥1 grid
# + autodiff narrow-width extension), and :ad (fast pure-autodiff onset). Shows why robust_ad is
# the recommended production model: it lowers sfmin into the narrow-width EP modes the grid misses,
# without the pure-:ad path's occasional cap/miss (e.g. IR=95).
#   julia --project=. dev/ad/plot_sfmin_grid_robust_ad_ad.jl
using Pkg
Pkg.activate(normpath(@__DIR__, "..", ".."))
using Plots
using Printf

const BUILD     = normpath(@__DIR__, "..")
const GRID_TXT  = get(ENV, "GRID_TXT",  joinpath(BUILD, "gacode_nb32_scan20_jgpu_53859063_tasks", "sfmin_scan.txt"))
const ROBUST_TXT= get(ENV, "ROBUST_TXT",joinpath(BUILD, "gacode_nb32_scan20_1node_robust_ad_54638544_tasks", "sfmin_scan.txt"))
const AD_TXT    = get(ENV, "AD_TXT",    joinpath(@__DIR__, "ad_threads_sfmin_nb32_ad.txt"))
const OUT_PNG   = joinpath(normpath(BUILD, ".."), "docs", "plots", "sfmin_grid_robust_ad_ad_nb32.png")

# read "idx ir sfmin" -> (ir, sfmin) vectors
function read_sfmin(path)
    ir = Int[]; sf = Float64[]
    for line in eachline(path)
        p = split(strip(line))
        length(p) >= 3 || continue
        push!(ir, parse(Int, p[2])); push!(sf, parse(Float64, p[3]))
    end
    return ir, sf
end

function relstats(ref_ir, ref_sf, ir, sf)
    m = Dict(zip(ref_ir, ref_sf))
    rels = Float64[]
    for (i, v) in zip(ir, sf)
        if haskey(m, i) && m[i] != 0
            push!(rels, abs(v - m[i]) / m[i])
        end
    end
    isempty(rels) ? (NaN, NaN) :
        (100 * sort(rels)[cld(length(rels), 2)], 100 * maximum(rels))
end

function main()
    gir, gsf = read_sfmin(GRID_TXT)
    rir, rsf = read_sfmin(ROBUST_TXT)
    air, asf = read_sfmin(AD_TXT)

    default(legendfontsize=10, guidefontsize=11, tickfontsize=9, dpi=200,
            fontfamily="Computer Modern")

    plt = plot(gir, gsf; marker=:circle, lw=2, ms=6, label="grid (Fortran-equiv., w≥1)",
               color=:dodgerblue, yscale=:log10,
               xlabel="radial index IR", ylabel="critical scale factor  sfmin",
               title="DIII-D  N_BASIS=32  SCAN_N=20 :  grid vs robust_ad vs ad",
               legend=:topleft, size=(900, 520))
    plot!(plt, air, asf; marker=:diamond, lw=2, ms=5, label="ad (width-extended)",
          color=:darkorange, linestyle=:dot)
    plot!(plt, rir, rsf; marker=:star5, lw=2.6, ms=7, label="robust_ad (production)",
          color=:firebrick)

    mg, xg = relstats(gir, gsf, rir, rsf)
    @printf("robust_ad vs grid: median |rel diff| = %.1f%%   max = %.1f%%\n", mg, xg)

    mkpath(dirname(OUT_PNG))
    savefig(plt, OUT_PNG)
    println("wrote ", OUT_PNG)
end

main()
