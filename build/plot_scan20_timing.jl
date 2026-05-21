#!/usr/bin/env julia
# Bar plot: SCAN_N=20 wall time vs N_BASIS for Fortran / Julia CPU / Julia GPU.
using Pkg
Pkg.activate(normpath(@__DIR__, ".."))

include(joinpath(@__DIR__, "collect_scan20_timing.jl"))
using Plots

const OUTDIR = joinpath(@__DIR__, "timing_runs")

function read_timing_csv(path::String)
    lines = filter(!isempty, strip.(readlines(path)))
    basis = Int[]
    fort, cpu, gpu = Float64[], Float64[], Float64[]
    for line in lines[2:end]
        p = split(line, ",")
        length(p) < 4 && continue
        push!(basis, parse(Int, p[1]))
        push!(fort, isempty(strip(p[2])) ? NaN : parse(Float64, p[2]))
        push!(cpu, isempty(strip(p[3])) ? NaN : parse(Float64, p[3]))
        push!(gpu, isempty(strip(p[4])) ? NaN : parse(Float64, p[4]))
    end
    return basis, fort, cpu, gpu
end

csv_path, _ = collect_scan20_timing!()
basis, fort, cpu, gpu = read_timing_csv(csv_path)

default(legendfontsize=10, guidefontsize=11, tickfontsize=9, dpi=200)
mkpath(OUTDIR)

ymax = maximum(filter(!isnan, vcat(fort, cpu, gpu)); init=0.0)
ymax = ymax > 0 ? ymax * 1.12 : 250.0

bw = 0.22
x = Float64.(basis)
p = plot(
    xlim=(minimum(x) - 0.6, maximum(x) + 0.6),
    ylim=(0, ymax),
    xlabel="N_BASIS (basis functions)",
    ylabel="Wall time (s)",
    title="SCAN_N=20 radial scan — DIII-D 202017C42_500ms_v3.1",
    legend=:topright,
    xticks=(x, string.(basis)),
)
bar!(p, x .- bw, fort; bar_width=bw, label="Fortran (10n, 20 MPI)", color=:steelblue)
bar!(p, x, cpu; bar_width=bw, label="Julia CPU (10n, 20 workers)", color=:darkorange)
bar!(p, x .+ bw, gpu; bar_width=bw, label="Julia GPU (5n, 20 tasks)", color=:seagreen)

out_png = joinpath(OUTDIR, "scan20_timing_by_nbasis.png")
savefig(p, out_png)
println("Wrote ", out_png)

p2 = plot(
    basis, fort; label="Fortran", marker=:circle, linewidth=2, markersize=6,
    xlabel="N_BASIS", ylabel="Time (s)", title="SCAN_N=20 timing vs N_BASIS",
    xticks=basis, legend=:topleft, ylim=(0, ymax),
)
plot!(p2, basis, cpu; label="Julia CPU", marker=:square, linewidth=2, markersize=6)
plot!(p2, basis, gpu; label="Julia GPU", marker=:diamond, linewidth=2, markersize=6)
out2 = joinpath(OUTDIR, "scan20_timing_lines.png")
savefig(p2, out2)
println("Wrote ", out2)

println("\nData ($csv_path):")
for i in eachindex(basis)
    println("  N_BASIS=$(basis[i]): Fortran=$(fort[i])  JuliaCPU=$(cpu[i])  JuliaGPU=$(gpu[i])")
end
