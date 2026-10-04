# Plot the grid → robust_ad → truth cost/fidelity ladder: critical scale factor sfmin(radius) for
# DIII-D N_BASIS=32, SCAN_N=20. The three lines are exactly the two components of the grid→truth gap.
#   - grid       : Fortran-equivalent kwscale_scan on the canonical w∈[1,2] box (the w≥1 reference)
#   - robust_ad  : critical_factor_robust with extend_width=true -- min(w≥1 grid-zoom, extended
#                  narrow-width locate) at nb=N_BASIS. grid→robust_ad is the WIDTH component (the 2–11×
#                  reduction from admitting the narrow EP-driven AE).
#   - truth      : critical_factor_truth = robust_ad + nbasis ladder ({32,40,48,56}) at the located
#                  optimum. robust_ad→truth is the (adverse) NBASIS component. Read directly from a
#                  full :truth scan20 run's sfmin_scan.txt.
# Set GRID_TXT / ROBUST_TXT / PROD_TXT to the three runs' sfmin_scan.txt (PROD_TXT = the :truth run).
#   module load julia/1.11.7 && export JULIA_DEPOT_PATH=$PSCRATCH/.julia
#   julia --project=. dev/ad/plot_grid_vs_truth.jl
using Pkg
Pkg.activate(normpath(@__DIR__, "..", ".."))
using Plots
using Printf

const BUILD      = normpath(@__DIR__, "..")
const GRID_TXT   = get(ENV, "GRID_TXT",
    joinpath(BUILD, "gacode_nb32_scan20_jgpu_53859063_tasks", "sfmin_scan.txt"))
const ROBUST_TXT = get(ENV, "ROBUST_TXT",
    joinpath(@__DIR__, "ad_threads_sfmin_nb32_robust_ad_r1.txt"))
const PROD_TXT   = get(ENV, "PROD_TXT",
    joinpath(BUILD, "gacode_nb32_scan20_jgpu_truth_54587438_tasks", "sfmin_scan.txt"))
const OUT_PNG    = get(ENV, "OUT_PNG",
    joinpath(normpath(BUILD, ".."), "docs", "plots", "sfmin_grid_vs_truth_nb32.png"))

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
    rir, rsf = read_sfmin(ROBUST_TXT)
    pir, psf = read_sfmin(PROD_TXT)

    gmap = Dict(zip(gir, gsf))
    rmap = Dict(zip(rir, rsf))
    pmap = Dict(zip(pir, psf))
    irs  = sort(collect(intersect(keys(rmap), keys(pmap), keys(gmap))))

    default(legendfontsize=10, guidefontsize=11, tickfontsize=9, dpi=200,
            fontfamily="Computer Modern")

    plt = plot(gir, gsf; marker=:circle, lw=2, ms=5, ls=:dot,
               label="grid (Fortran-faithful, w≥1)",
               color=:gray55, yscale=:log10,
               xlabel="radial index IR", ylabel="critical scale factor  sfmin",
               title="DIII-D  N_BASIS=32  SCAN_N=20 :  grid / robust_ad / production truth",
               legend=:topleft, size=(940, 540))
    plot!(plt, rir, rsf; marker=:circle, lw=2, ms=6,
          label="robust_ad (refined w≥1, floor)", color=:dodgerblue)
    plot!(plt, pir, psf; marker=:utriangle, lw=2.5, ms=6, ls=:dash,
          label="production  min(robust_ad, truth)  [:truth, nb→56]", color=:seagreen)

    # where the production min came from (floored vs truth-won) and how far below grid it sits
    floored = count(ir -> isapprox(pmap[ir], rmap[ir]; rtol=1e-3), irs)
    gratios = [gmap[ir] / pmap[ir] for ir in irs]      # grid / production
    order   = sortperm(gratios; rev=true)
    @printf("production floors on robust_ad at %d/%d radii; truth wins %d.\n",
            floored, length(irs), length(irs) - floored)
    @printf("grid/production ratio: median=%.2fx  max=%.2fx (IR=%d)\n",
            sort(gratios)[cld(length(gratios),2)], maximum(gratios), irs[argmax(gratios)])
    @printf("top-5 production reductions vs grid: %s\n",
            join([@sprintf("IR%d=%.1fx", irs[i], gratios[i]) for i in order[1:min(5,end)]], "  "))

    mkpath(dirname(OUT_PNG))
    savefig(plt, OUT_PNG)
    println("wrote ", OUT_PNG)
end

main()
