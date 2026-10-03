#!/usr/bin/env julia
# UCP_complete accuracy at N_BASIS=32: sfmin (critical EP factor) vs radial grid index IR,
# for the reactor-relevant case, in the same legend order / per-solver colors as the DIII-D
# figure (dev/ad/plot_ad_wide_nb32.jl):
#   * Fortran (CPU -n 1280; SFmin block from the newest ucp_fortran_scan20_nb32_* run dir)
#   * Julia :grid       (5N MPS)      ucp_nb32_scan20_jgpu_*_tasks/sfmin_scan.txt
#   * Julia :ad :only   (5N threads)  ucp_nb32_scan20_jgpu_ad_only_*_tasks/sfmin_scan.txt
#   * Julia :ad :locate (1N backfill) ucp_nb32_scan20_1node_ad_locate_*_tasks/sfmin_scan.txt
#   * Julia :ad :wide   (1N backfill) ucp_nb32_scan20_1node_ad_wide_*_tasks/sfmin_scan.txt
# Sentinel sfmin (>=9000, stable/rejected edge radii) are masked to NaN so the log axis stays
# on the physical range.
#   cd dev && julia --project=.. ad/plot_ucp_accuracy_nb32.jl
using Pkg
Pkg.activate(normpath(@__DIR__, "..", ".."))
using Plots, Printf
using Plots.PlotMeasures: mm

const BUILD = normpath(@__DIR__, "..")
const SENTINEL = 9000.0

_jobid(f) = (m = match(r"_(\d+)(?:_tasks)?$", basename(f)); m === nothing ? -1 : parse(Int, m.captures[1]))

# newest directory in BUILD matching a prefix..suffix (jobid embedded), or nothing.
function newest_dir(rx::Regex)
    dirs = [joinpath(BUILD, d) for d in readdir(BUILD) if isdir(joinpath(BUILD, d)) && occursin(rx, d)]
    isempty(dirs) && return nothing
    return first(sort(dirs; by=_jobid, rev=true))
end

read_sfmin(p) = begin
    ir = Int[]; f = Float64[]
    for ln in eachline(p)
        t = split(strip(ln)); length(t) >= 3 || continue
        push!(ir, parse(Int, t[2]))
        v = parse(Float64, t[3]); push!(f, v >= SENTINEL ? NaN : v)
    end
    ir, f
end

read_fortran_sfmin(p) = begin
    vals = Float64[]; inblk = false
    for ln in eachline(p)
        s = strip(ln)
        s == "SFmin" && (inblk = true; continue)
        if inblk
            startswith(s, "-") && break
            v = tryparse(Float64, s)
            v === nothing || push!(vals, v >= SENTINEL ? NaN : v)
        end
    end
    vals
end

sfmin_path(rx) = (d = newest_dir(rx); d === nothing ? nothing : joinpath(d, "sfmin_scan.txt"))

series = [
    ("Julia :grid (GPU)",              sfmin_path(r"^ucp_nb32_scan20_jgpu_\d+_tasks$"),          :gray55,    :dash,  :diamond),
    ("Julia :ad :only (GPU)",          sfmin_path(r"^ucp_nb32_scan20_jgpu_ad_only_\d+_tasks$"),  :darkorange,:dot,   :rect),
    ("Julia :ad :locate (GPU, default)", sfmin_path(r"^ucp_nb32_scan20_1node_ad_locate_\d+_tasks$"), :seagreen, :solid, :utriangle),
    ("Julia :ad :wide (GPU)",          sfmin_path(r"^ucp_nb32_scan20_1node_ad_wide_\d+_tasks$"), :firebrick, :solid, :star5),
]

fdir = newest_dir(r"^ucp_fortran_scan20_nb32_\d+$")
# fortran run dirs live under timing_runs/, not build/ root:
if fdir === nothing
    tr = joinpath(BUILD, "timing_runs")
    if isdir(tr)
        cands = [joinpath(tr, d) for d in readdir(tr) if occursin(r"^ucp_fortran_scan20_nb32_\d+$", d)]
        fdir = isempty(cands) ? nothing : first(sort(cands; by=_jobid, rev=true))
    end
end
FORT = fdir === nothing ? nothing : joinpath(fdir, "out.TGLFEP")

default(legendfontsize=9, guidefontsize=11, tickfontsize=9, dpi=200)
p = plot(xlabel="radial grid index IR", ylabel="sfmin (critical factor)",
         title="UCP_complete SCAN_N=20 accuracy @ N_BASIS=32",
         yscale=:log10, legend=:topleft, size=(950, 560),
         left_margin=4mm, bottom_margin=4mm)

# reference IR axis from the grid series (fall back to any available series).
ir_ref = Int[]
for (_, path, _, _, _) in series
    global ir_ref
    path !== nothing && isfile(path) && (ir_ref = first(read_sfmin(path)); break)
end

if FORT !== nothing && isfile(FORT) && !isempty(ir_ref)
    fort = read_fortran_sfmin(FORT)
    if length(fort) == length(ir_ref)
        plot!(p, ir_ref, fort; label="Fortran (CPU -n 1280)", color=:steelblue, linestyle=:solid,
              marker=:circle, markersize=6, linewidth=3)
    else
        @warn "Fortran SFmin length mismatch; skipping" nfort=length(fort) nir=length(ir_ref)
    end
else
    @warn "Fortran out.TGLFEP not found; plotting Julia series only" FORT
end

for (lab, path, col, ls, mk) in series
    (path === nothing || !isfile(path)) && (@warn "missing $lab ($path)"; continue)
    ir, f = read_sfmin(path)
    plot!(p, ir, f; label=lab, color=col, linestyle=ls, marker=mk, markersize=5, linewidth=2)
end

out = joinpath(@__DIR__, "ucp_accuracy_nb32.png")
savefig(p, out)
println("Wrote ", out)

# numeric table (Fortran vs grid vs locate vs wide)
if !isempty(ir_ref)
    fort = (FORT !== nothing && isfile(FORT)) ? read_fortran_sfmin(FORT) : Float64[]
    _, grd = series[1][2] !== nothing && isfile(series[1][2]) ? read_sfmin(series[1][2]) : (Int[], Float64[])
    _, aon = series[2][2] !== nothing && isfile(series[2][2]) ? read_sfmin(series[2][2]) : (Int[], Float64[])
    _, loc = series[3][2] !== nothing && isfile(series[3][2]) ? read_sfmin(series[3][2]) : (Int[], Float64[])
    _, wid = series[4][2] !== nothing && isfile(series[4][2]) ? read_sfmin(series[4][2]) : (Int[], Float64[])
    println("\n IR    Fortran      grid   ad:only     locate      wide")
    for i in eachindex(ir_ref)
        g(v) = i <= length(v) ? v[i] : NaN
        @printf("%4d  %9.4f  %9.4f  %9.4f  %9.4f  %9.4f\n", ir_ref[i], g(fort), g(grd), g(aon), g(loc), g(wid))
    end
end
