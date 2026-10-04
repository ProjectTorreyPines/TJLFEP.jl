#!/usr/bin/env julia
# Node-hours vs N_BASIS for the SCAN_N=20 DIII-D scan. Plots the 2x2 solver-selection
# options (+ the Fortran baseline) so the figure maps onto the README decision matrix:
#   * Fortran (CPU, reference baseline; its own node count from the CSV)
#   * Julia :grid       (julia_gpu_nh)   -- matches Fortran (w>=1 box)
#   * Julia :ad :only   -- approximates :grid (bare w>=1 AD), parsed from the ONLY logs
#   * Julia :ad :locate (julia_gpu_ad_nh)-- extends Fortran with narrow w<1 AE modes (default)
#   * Julia :ad :wide   (kdesc=2)        -- fast approximation of :locate, parsed from WIDE logs
# The robust_ad/truth reference tiers are intentionally omitted: :ad :locate matches their
# accuracy at lower cost, so the plot guides users along the four matrix options.
# Node-hours (nodes x wallclock) is the fair metric regardless of node count (Fortran CPU on
# 10 nodes, the GPU tiers on 5), so the legend just tags (CPU)/(GPU).
#   julia --project=. dev/timing/plot_scan20_timing_wide.jl
using Pkg
Pkg.activate(normpath(@__DIR__, "..", ".."))
using Plots, Printf
using Plots.PlotMeasures: mm

const BUILD = normpath(@__DIR__, "..")
const CSV   = joinpath(BUILD, "timing_runs", "scan20_timing.csv")
const OUT   = joinpath(BUILD, "timing_runs", "scan20_timing_wide_lines.png")

# --- established series from the node-hours CSV (col 1 = n_basis; node-hours block at cols 10..17) ---
num(p, i) = (length(p) >= i && !isempty(strip(p[i]))) ? parse(Float64, p[i]) : NaN
basis = Int[]; fort = Float64[]; grid = Float64[]; locate = Float64[]; robust = Float64[]; truth = Float64[]
for ln in readlines(CSV)[2:end]
    isempty(strip(ln)) && continue
    p = split(ln, ",")
    length(p) < 17 && continue
    push!(basis, parse(Int, p[1]))
    push!(fort,   num(p, 10))   # fortran_nh
    push!(grid,   num(p, 12))   # julia_gpu_nh (grid MPS)
    push!(locate, num(p, 14))   # julia_gpu_ad_nh (:ad :locate, 1-node backfill)
    push!(truth,  num(p, 16))   # julia_gpu_truth_nh
    push!(robust, num(p, 17))   # julia_gpu_robust_ad_nh
end

# --- :wide series: parse node-hours from the 1-node-backfill WIDE logs (scan-phase seconds, nodes=1) ---
function wide_nh(nb::Int)
    pat = Regex("time_scan20_nb$(nb)_julia_gpu_ad_WIDE_[0-9]+\\.out")
    files = sort(filter(f -> occursin(pat, f), readdir(BUILD)))
    isempty(files) && return NaN
    for f in reverse(files)                       # newest job id wins
        for ln in eachline(joinpath(BUILD, f))
            occursin("phase=scan", ln) || continue
            m = match(r"seconds=\s*([0-9.]+)", ln)
            n = match(r"nodes=([0-9]+)", ln)
            m === nothing && continue
            nodes = n === nothing ? 1.0 : parse(Float64, n.captures[1])
            return parse(Float64, m.captures[1]) * nodes / 3600
        end
    end
    return NaN
end
wide = [wide_nh(nb) for nb in basis]

# --- :ad :only (bare w>=1, no confirm) series ---
# 5-node MPS layout (5 nodes x 4 GPU = 20 GPUs, one radius/task, wall = slowest task),
# node-matched to grid/:wide. Node-hours = nodes x full-scan wallclock. Sources are
# PINNED per N_BASIS (job id) for reproducibility -- "newest-wins" globbing previously
# grabbed stray reruns. nb8/16/32 use the FILE-ONLY GPU sysimage runs (leaner 1.1 GB
# image -> faster load + identical compute, so faster than the generic image). nb6 uses
# the GENERIC-image run: its file-only scan wallclock is dominated by 20-way cold 1.1 GB
# sysimage-load I/O (~150-170 s, reproducible over two runs) rather than compute
# (~13 s/radius), which would make the line non-monotone. Mixed-image by design: :only is
# file-only where clean; grid/:locate/:wide stay on the generic image.
const ONLY_JOBID = Dict(6 => 55326426, 8 => 55330924, 16 => 55330925, 32 => 55330927)
function only_nh(nb::Int)
    haskey(ONLY_JOBID, nb) || return NaN
    f = joinpath(BUILD, "time_scan20_nb$(nb)_julia_gpu_ad_ONLY_$(ONLY_JOBID[nb]).out")
    isfile(f) || (@warn "missing pinned :only log" file=f; return NaN)
    for ln in eachline(f)
        occursin("phase=scan", ln) || continue
        m = match(r"seconds=\s*([0-9.]+)", ln)
        n = match(r"nodes=([0-9]+)", ln)
        m === nothing && continue
        nodes = n === nothing ? 5.0 : parse(Float64, n.captures[1])
        return parse(Float64, m.captures[1]) * nodes / 3600
    end
    return NaN
end
only = [only_nh(nb) for nb in basis]

default(legendfontsize=9, guidefontsize=11, tickfontsize=9, dpi=200, fontfamily="Computer Modern")
# Plot only the 2x2 solver-selection options (+ Fortran baseline): the reference tiers
# robust_ad/truth are omitted here because :ad :locate matches their accuracy at lower cost.
ys = filter(!isnan, vcat(fort, grid, locate, wide, only))
ymax = (isempty(ys) ? 1.0 : maximum(ys)) * 1.12
p = plot(xlabel="N_BASIS", ylabel="Node-hours (nodes × wallclock)",
         title="SCAN_N=20 node-hours vs N_BASIS: DIII-D",
         xticks=basis, legend=:topleft, ylim=(0, ymax), size=(900, 560),
         left_margin=4mm, bottom_margin=4mm)
plot!(p, basis, fort;   label="Fortran (CPU)",          marker=:circle,   linewidth=2, markersize=6, color=:steelblue)
plot!(p, basis, grid;   label="Julia :grid (GPU)",      marker=:diamond,  linewidth=2, markersize=6, color=:gray55, linestyle=:dash)
any(!isnan, only) && plot!(p, basis, only; label="Julia :ad :only (GPU)", marker=:rect, linewidth=2, markersize=5, color=:darkorange, linestyle=:dot)
plot!(p, basis, locate; label="Julia :ad :locate (GPU)",marker=:utriangle,linewidth=2, markersize=6, color=:seagreen,   linestyle=:solid)
plot!(p, basis, wide;   label="Julia :ad :wide (GPU)", marker=:star5, linewidth=2.5, markersize=8, color=:firebrick)
savefig(p, OUT)
println("Wrote ", OUT)

println("\n N_BASIS  Fortran   :grid    :only    :locate   :wide   (locate/wide)")
for i in eachindex(basis)
    c(x) = isnan(x) ? "  —   " : @sprintf("%6.3f", x)
    rat = (!isnan(locate[i]) && !isnan(wide[i]) && wide[i] > 0) ? @sprintf("%5.2fx", locate[i]/wide[i]) : "  —  "
    @printf("  %4d   %s  %s  %s  %s  %s   %s\n", basis[i], c(fort[i]), c(grid[i]), c(only[i]), c(locate[i]), c(wide[i]), rat)
end
