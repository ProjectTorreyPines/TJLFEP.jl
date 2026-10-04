# Plot grid (Fortran-equivalent kwscale_scan) vs AD critical scale factor sfmin(radius) for
# DIII-D N_BASIS=32, SCAN_N=20. Grid from the scan20 grid run; AD from the profile written by
# ad_threads_sfmin_profile.jl. By default plots the hardened robust path
# (ad_threads_sfmin_nb32_robust_ad_r1.txt); override the AD file/label with AD_TXT / AD_LABEL.
#   julia --project=. dev/ad/plot_grid_vs_ad_threads.jl
using Pkg
Pkg.activate(normpath(@__DIR__, "..", ".."))
using Plots
using Printf

const BUILD    = normpath(@__DIR__, "..")
const GRID_TXT = joinpath(BUILD, "gacode_nb32_scan20_jgpu_53859063_tasks", "sfmin_scan.txt")
const AD_TXT   = get(ENV, "AD_TXT", joinpath(@__DIR__, "ad_threads_sfmin_nb32_robust_ad_r1.txt"))
const AD_LABEL = get(ENV, "AD_LABEL", "robust_ad (threads, refine=1)")
const OUT_PNG  = joinpath(normpath(BUILD, ".."), "docs", "plots", "sfmin_grid_vs_ad_threads_nb32.png")

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

function main()
    gir, gsf = read_sfmin(GRID_TXT)
    air, asf = read_sfmin(AD_TXT)

    default(legendfontsize=10, guidefontsize=11, tickfontsize=9, dpi=200,
            fontfamily="Computer Modern")

    plt = plot(gir, gsf; marker=:circle, lw=2, ms=6, label="grid (Fortran-equiv.)",
               color=:dodgerblue, yscale=:log10,
               xlabel="radial index IR", ylabel="critical scale factor  sfmin",
               title="DIII-D  N_BASIS=32  SCAN_N=20 :  grid vs AD",
               legend=:topleft, size=(900, 520))
    plot!(plt, air, asf; marker=:diamond, lw=2, ms=6, label=AD_LABEL,
          color=:darkorange)

    # annotate per-radius relative difference where both exist (align by IR)
    gmap = Dict(zip(gir, gsf))
    rels = Float64[]
    for (ir, a) in zip(air, asf)
        if haskey(gmap, ir) && gmap[ir] != 0
            push!(rels, abs(a - gmap[ir]) / gmap[ir])
        end
    end
    if !isempty(rels)
        @printf("median |rel diff| = %.1f%%   max = %.1f%%   over %d radii\n",
                100*sort(rels)[cld(length(rels),2)], 100*maximum(rels), length(rels))
    end

    mkpath(dirname(OUT_PNG))
    savefig(plt, OUT_PNG)
    println("wrote ", OUT_PNG)
end

main()
