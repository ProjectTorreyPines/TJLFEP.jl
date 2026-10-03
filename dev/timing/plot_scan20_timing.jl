#!/usr/bin/env julia
# Bar plot: SCAN_N=20 wall time vs N_BASIS for Fortran / Julia CPU / Julia GPU.
using Pkg
Pkg.activate(normpath(@__DIR__, "..", ".."))

include(joinpath(@__DIR__, "collect_scan20_timing.jl"))
using Plots
using Printf

const OUTDIR = joinpath(normpath(@__DIR__, ".."), "timing_runs")

# Reads the CSV into (seconds NamedTuple, node-hours NamedTuple). Current schema (18 cols):
#   1=n_basis | 2..9 = wallclock seconds | 10..17 = node-hours (nodes×s/3600) | 18 = notes.
# Older schemas (no node-hours block, fewer series) are tolerated: every field is guarded by count
# and missing values become NaN. node-hours falls back to NaN for pre-node-hours CSVs.
function read_timing_csv(path::String)
    lines = filter(!isempty, strip.(readlines(path)))
    mk() = (fort=Float64[], cpu=Float64[], gpu=Float64[], cpu_ad=Float64[],
            gpu_ad=Float64[], gpu_ad_mps=Float64[], gpu_truth=Float64[], gpu_robust=Float64[])
    basis = Int[]; s = mk(); nh = mk()
    num(p, i) = (length(p) >= i && !isempty(strip(p[i]))) ? parse(Float64, p[i]) : NaN
    for line in lines[2:end]
        p = split(line, ",")
        length(p) < 4 && continue
        push!(basis, parse(Int, p[1]))
        # seconds @ 2..9
        push!(s.fort, num(p, 2)); push!(s.cpu, num(p, 3)); push!(s.gpu, num(p, 4))
        push!(s.cpu_ad,     length(p) >= 7  ? num(p, 5) : NaN)
        push!(s.gpu_ad,     length(p) >= 7  ? num(p, 6) : NaN)
        push!(s.gpu_ad_mps, length(p) >= 8  ? num(p, 7) : NaN)
        push!(s.gpu_truth,  length(p) >= 9  ? num(p, 8) : NaN)
        push!(s.gpu_robust, length(p) >= 10 ? num(p, 9) : NaN)
        # node-hours @ 10..17 (only present in the 18-col schema)
        nh18 = length(p) >= 18
        push!(nh.fort,       nh18 ? num(p, 10) : NaN)
        push!(nh.cpu,        nh18 ? num(p, 11) : NaN)
        push!(nh.gpu,        nh18 ? num(p, 12) : NaN)
        push!(nh.cpu_ad,     nh18 ? num(p, 13) : NaN)
        push!(nh.gpu_ad,     nh18 ? num(p, 14) : NaN)
        push!(nh.gpu_ad_mps, nh18 ? num(p, 15) : NaN)
        push!(nh.gpu_truth,  nh18 ? num(p, 16) : NaN)
        push!(nh.gpu_robust, nh18 ? num(p, 17) : NaN)
    end
    return basis, s, nh
end

csv_path, _ = collect_scan20_timing!()
basis, S, NH = read_timing_csv(csv_path)
# Wallclock minutes (for the bar plot).
fort = S.fort ./ 60; cpu = S.cpu ./ 60; gpu = S.gpu ./ 60
cpu_ad = S.cpu_ad ./ 60; gpu_ad = S.gpu_ad ./ 60; gpu_ad_mps = S.gpu_ad_mps ./ 60
gpu_truth = S.gpu_truth ./ 60; gpu_robust = S.gpu_robust ./ 60
# Node-hours (for the timing-vs-nbasis line plot). Already in node-hours; no unit conversion.
nh_fort = NH.fort; nh_cpu = NH.cpu; nh_gpu = NH.gpu
nh_cpu_ad = NH.cpu_ad; nh_gpu_ad = NH.gpu_ad; nh_gpu_ad_mps = NH.gpu_ad_mps
nh_gpu_truth = NH.gpu_truth; nh_gpu_robust = NH.gpu_robust

# The Julia CPU sysimage runs may still be pending; the collector then backfills STALE
# pre-optimization numbers. Set SKIP_JULIA_CPU=1 to omit the Julia CPU series entirely
# (Fortran + Julia GPU only) until the fresh CPU runs complete.
const SHOW_CPU = get(ENV, "SKIP_JULIA_CPU", "0") != "1"

default(legendfontsize=10, guidefontsize=11, tickfontsize=9, dpi=200,
        fontfamily="Computer Modern")
mkpath(OUTDIR)

# AD series (line plot only). Simplified to show ONLY the fastest AD route (GPU in-process
# threads); AD-CPU and AD-GPU+MPS are measured but hidden to keep the plot readable.
const SHOW_AD_CPU = false
const SHOW_AD_GPU = any(!isnan, gpu_ad)
const SHOW_AD_GPU_MPS = false
# Physical-truth (solver=:truth) GPU+MPS series: the production min(robust_ad, truth) profile.
const SHOW_TRUTH = any(!isnan, gpu_truth)
# Robust_ad (solver=:robust_ad, width-extended) GPU+MPS series: the WIDTH tier (no nbasis ladder).
const SHOW_ROBUST = any(!isnan, gpu_robust)

ymax_series = SHOW_CPU ? vcat(fort, cpu, gpu) : vcat(fort, gpu)
SHOW_AD_CPU && (ymax_series = vcat(ymax_series, cpu_ad))
SHOW_AD_GPU && (ymax_series = vcat(ymax_series, gpu_ad))
SHOW_AD_GPU_MPS && (ymax_series = vcat(ymax_series, gpu_ad_mps))
SHOW_TRUTH && (ymax_series = vcat(ymax_series, gpu_truth))
SHOW_ROBUST && (ymax_series = vcat(ymax_series, gpu_robust))
ymax = maximum(filter(!isnan, ymax_series); init=0.0)
ymax = ymax > 0 ? ymax * 1.12 : 30.0

# Separate y-axis ceiling for the node-hours line plot.
ymax_nh_series = SHOW_CPU ? vcat(nh_fort, nh_cpu, nh_gpu) : vcat(nh_fort, nh_gpu)
SHOW_AD_GPU && (ymax_nh_series = vcat(ymax_nh_series, nh_gpu_ad))
SHOW_TRUTH && (ymax_nh_series = vcat(ymax_nh_series, nh_gpu_truth))
SHOW_ROBUST && (ymax_nh_series = vcat(ymax_nh_series, nh_gpu_robust))
ymax_nh = maximum(filter(!isnan, ymax_nh_series); init=0.0)
ymax_nh = ymax_nh > 0 ? ymax_nh * 1.12 : 10.0

bw = 0.22
x = Float64.(basis)
p = plot(
    xlim=(minimum(x) - 0.6, maximum(x) + 0.6),
    ylim=(0, ymax),
    xlabel="N_BASIS (basis functions)",
    ylabel="Total wallclock (min)",
    title="SCAN_N=20 radial scan — DIII-D 202017C42_500ms_v3.1",
    legend=:topright,
    xticks=(x, string.(basis)),
)
if SHOW_CPU
    bar!(p, x .- bw, fort; bar_width=bw, label="Fortran (10 CPU nodes)", color=:steelblue)
    bar!(p, x, cpu; bar_width=bw, label="Julia (10 CPU nodes)", color=:darkorange)
    bar!(p, x .+ bw, gpu; bar_width=bw, label="Julia (5 GPU nodes)", color=:seagreen)
else
    # Recenter the two remaining series so there's no empty gap where CPU would be.
    bar!(p, x .- bw/2, fort; bar_width=bw, label="Fortran (10 CPU nodes)", color=:steelblue)
    bar!(p, x .+ bw/2, gpu; bar_width=bw, label="Julia (5 GPU nodes)", color=:seagreen)
end

out_png = joinpath(OUTDIR, "scan20_timing_by_nbasis.png")
savefig(p, out_png)
println("Wrote ", out_png)

# Timing-vs-nbasis plot, now in NODE-HOURS (nodes × wallclock): the resource-cost metric for
# NN-database generation. node-hours auto-uses each series' real node count (5-node wave runs vs
# 1-node backfill runs), so a single-node backfill run plots at its true (much lower) cost.
p2 = plot(
    basis, nh_fort; label="Fortran (10 CPU nodes)", marker=:circle, linewidth=2, markersize=6,
    xlabel="N_BASIS", ylabel="Node-hours (nodes × wallclock)", title="SCAN_N=20 node-hours vs N_BASIS",
    xticks=basis, legend=:topleft, ylim=(0, ymax_nh),
)
SHOW_CPU && plot!(p2, basis, nh_cpu; label="Julia grid (10 CPU nodes)", marker=:square, linewidth=2, markersize=6)
plot!(p2, basis, nh_gpu; label="Julia grid MPS (GPU)", marker=:diamond, linewidth=2, markersize=6)
SHOW_AD_CPU && plot!(p2, basis, nh_cpu_ad; label="Julia AD (10 CPU nodes)", marker=:utriangle, linewidth=2, markersize=6, linestyle=:dash)
SHOW_AD_GPU && plot!(p2, basis, nh_gpu_ad; label="Julia AD threads (GPU)", marker=:star5, linewidth=2, markersize=7, linestyle=:dash)
SHOW_AD_GPU_MPS && plot!(p2, basis, nh_gpu_ad_mps; label="Julia AD GPU+MPS", marker=:hexagon, linewidth=2, markersize=6, linestyle=:dot)
SHOW_ROBUST && plot!(p2, basis, nh_gpu_robust; label="Julia robust_ad MPS (GPU)", marker=:star6, linewidth=2, markersize=7, linestyle=:dash, color=:firebrick)
SHOW_TRUTH && plot!(p2, basis, nh_gpu_truth; label="Julia truth MPS (GPU)", marker=:pentagon, linewidth=2, markersize=7, linestyle=:dashdot, color=:purple)
out2 = joinpath(OUTDIR, "scan20_timing_lines.png")
savefig(p2, out2)
println("Wrote ", out2)

println("\nData ($csv_path), node-hours", SHOW_CPU ? "" : " (Julia CPU omitted: SKIP_JULIA_CPU=1)", ":")
for i in eachindex(basis)
    @printf("  N_BASIS=%d: Fortran=%s  JuliaCPU=%s  JuliaGPU=%s  AD-GPU=%s  Robust_ad=%s  Truth=%s\n", basis[i],
        isnan(nh_fort[i]) ? "—" : @sprintf("%.2f", nh_fort[i]),
        (!SHOW_CPU || isnan(nh_cpu[i])) ? "—" : @sprintf("%.2f", nh_cpu[i]),
        isnan(nh_gpu[i]) ? "—" : @sprintf("%.2f", nh_gpu[i]),
        isnan(nh_gpu_ad[i]) ? "—" : @sprintf("%.2f", nh_gpu_ad[i]),
        isnan(nh_gpu_robust[i]) ? "—" : @sprintf("%.2f", nh_gpu_robust[i]),
        isnan(nh_gpu_truth[i]) ? "—" : @sprintf("%.2f", nh_gpu_truth[i]))
end
