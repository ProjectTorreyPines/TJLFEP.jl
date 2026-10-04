using Pkg; Pkg.activate(normpath(@__DIR__, "..", ".."))
using Plots, Printf
using Plots.PlotMeasures: mm
gr()

# Adaptive-union batched SI timing vs N_BASIS (README timing-vs-nbasis framing).
# CSV columns: nb,n,npencils,ngpu,ms_fixed,ms_adaptive,net_speedup,misses
CSV = get(ENV, "CSV", normpath(@__DIR__, "adaptive_si_vs_nbasis.csv"))
nb=Int[]; n=Int[]; msf=Float64[]; msa=Float64[]; spd=Float64[]
for (k, ln) in enumerate(eachline(CSV))
    k == 1 && continue
    t = split(strip(ln), ','); length(t) >= 7 || continue
    push!(nb, parse(Int,t[1])); push!(n, parse(Int,t[2]))
    push!(msf, parse(Float64,t[5])); push!(msa, parse(Float64,t[6])); push!(spd, parse(Float64,t[7]))
end
isempty(nb) && error("no rows in $CSV")
println(" nb    n   ms/pencil(adaptive)  net_speedup_vs_geev")
for i in eachindex(nb); @printf("%3d %5d   %8.2f            %6.1fx\n", nb[i], n[i], msa[i], spd[i]); end

default(legendfontsize=9, guidefontsize=11, tickfontsize=9, dpi=200)
p = plot(nb, msa; marker=:circle, ms=6, lw=2, color=:royalblue, label="adaptive-union batched SI",
         xlabel="N_BASIS", ylabel="ms / pencil (4× A100)", xscale=:log2, yscale=:log10,
         title="Adaptive-union batched SI: eigensolve time vs N_BASIS",
         legend=:topleft, size=(900,540), left_margin=5mm, bottom_margin=5mm)
plot!(p, nb, msf; marker=:diamond, ms=5, lw=2, ls=:dash, color=:gray55, label="fixed-shift batched SI")
xticks!(p, Float64.(nb), string.(nb))
for i in eachindex(nb)
    annotate!(p, nb[i], msa[i]*1.18, text(@sprintf("%.0f×", spd[i]), 7, :seagreen, :center))
end
out = joinpath(@__DIR__, "adaptive_si_timing_vs_nbasis.png")
savefig(p, out); println("wrote ", out, "  (green = net speedup vs full geev, incl. calibration)")
