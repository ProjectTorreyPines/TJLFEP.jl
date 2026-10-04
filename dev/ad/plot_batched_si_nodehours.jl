using Pkg; Pkg.activate(normpath(@__DIR__, "..", ".."))
using Plots, Printf
using Plots.PlotMeasures: mm
gr()

# Node-hours vs N_BASIS for the 20-radius scan (README metric). Overlays the wired inner=:batched_si
# (4-GPU sharded, dev/ad/batched_si_nodehours.csv) on the reference tiers.
#
# Fortran (col10) and :grid GPU (col12) in docs/plots/scan20_timing.csv are current (they match the
# main TJLFEP/README node-hours: grid nb32 ~0.31). The :ad column in that CSV (col6/14) is STALE — it
# predates the confirm solves, so it *decreases* with N_BASIS. The current :ad :only *grows* with
# N_BASIS. To match the main README node-hours plot EXACTLY we parse :ad :only from the same PINNED
# 5-GPU-node run logs it uses (dev/timing/plot_scan20_timing_wide.jl's ONLY_JOBID), node-hours =
# nodes * scan_seconds / 3600. (:ad :only has NO backfill run — like :grid its per-radius cost is
# uniform so 5-node MPS is optimal; only :ad :locate/:wide use 1-node backfill.) Cross-check: nb32
# 5*132.2/3600 = 0.184 == README's "~0.18" and Fortran 4.29/0.184 ~= 23x.
BUILD = normpath(@__DIR__, "..")
REF   = joinpath(BUILD, "..", "docs", "plots", "scan20_timing.csv")
BSI   = get(ENV, "BSI", joinpath(@__DIR__, "batched_si_nodehours.csv"))
TIMING = joinpath(BUILD, "timing")   # pinned :only 5-node logs live in build/

# Pinned :ad :only 5-node run per N_BASIS (mirrors plot_scan20_timing_wide.jl so both plots agree).
const ONLY_JOBID = Dict(6 => 55326426, 8 => 55330924, 16 => 55330925, 32 => 55330927)
function ad_only_nh(nb)
    haskey(ONLY_JOBID, nb) || return NaN
    f = joinpath(BUILD, "time_scan20_nb$(nb)_julia_gpu_ad_ONLY_$(ONLY_JOBID[nb]).out")
    isfile(f) || return NaN
    for ln in eachline(f)
        occursin("phase=scan", ln) || continue
        m = match(r"seconds=\s*([0-9.]+)", ln); n = match(r"nodes=([0-9]+)", ln)
        m === nothing && continue
        nodes = n === nothing ? 5.0 : parse(Float64, n.captures[1])
        return parse(Float64, m.captures[1]) * nodes / 3600
    end
    return NaN
end

refnb=Int[]; fort=Float64[]; grid=Float64[]; ad=Float64[]
for (k, ln) in enumerate(eachline(REF))
    k == 1 && continue
    t = split(strip(ln), ','); length(t) >= 14 || continue
    nb = parse(Int,t[1])
    push!(refnb, nb); push!(fort, parse(Float64,t[10]))
    push!(grid, parse(Float64,t[12])); push!(ad, ad_only_nh(nb))
end
bnb=Int[]; bnh=Float64[]; bwall=Float64[]
for (k, ln) in enumerate(eachline(BSI))
    k == 1 && continue
    t = split(strip(ln), ','); length(t) >= 5 || continue
    push!(bnb, parse(Int,t[1])); push!(bwall, parse(Float64,t[4])); push!(bnh, parse(Float64,t[5]))
end

println(" nb   Fortran_nh  grid_nh   ad_nh   batched_si_nh  (bsi/grid)")
for i in eachindex(bnb)
    j = findfirst(==(bnb[i]), refnb)
    g = j===nothing ? NaN : grid[j]
    @printf("%3d  %9.3f  %7.3f  %6.3f  %11.3f    %5.2fx\n",
            bnb[i], j===nothing ? NaN : fort[j], g, j===nothing ? NaN : ad[j], bnh[i], bnh[i]/g)
end

default(legendfontsize=9, guidefontsize=11, tickfontsize=9, dpi=200)
p = plot(refnb, fort; marker=:circle, ms=6, lw=3, color=:steelblue, label="Fortran (CPU)",
         xlabel="N_BASIS", ylabel="node-hours (20-radius scan)", xscale=:log2, yscale=:log10,
         title="Node-hours vs N_BASIS (SCAN_N=20): batched_si vs grid/ad/Fortran",
         legend=:topleft, size=(950,560), left_margin=6mm, bottom_margin=5mm)
plot!(p, refnb, grid; marker=:diamond, ms=5, lw=2, ls=:dash, color=:gray55, label="Julia :grid (GPU)")
plot!(p, refnb, ad;   marker=:utriangle, ms=5, lw=2, color=:seagreen, label="Julia :ad :only (GPU)")
plot!(p, bnb, bnh;    marker=:star5, ms=8, lw=2, color=:royalblue, label="Julia :batched_si (4-GPU sharded)")
xticks!(p, Float64.(refnb), string.(refnb))
out = joinpath(@__DIR__, "batched_si_nodehours_vs_nbasis.png")
savefig(p, out); println("wrote ", out)
