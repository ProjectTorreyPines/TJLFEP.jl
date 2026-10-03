using Pkg; Pkg.activate(normpath(@__DIR__, "..", ".."))
using Plots, Printf
using Plots.PlotMeasures: mm
gr()

# Overlay the wired inner=:batched_si (fixed-shift hybrid) sfmin(IR) on the stored N_BASIS=32
# reference: Fortran (out.TGLFEP SFmin block) and Julia :grid (reproduces Fortran bit-for-bit).
BUILD = normpath(@__DIR__, "..")
GRIDF = joinpath(BUILD, "gacode_nb32_scan20_jgpu_53859063_tasks", "sfmin_scan.txt")   # idx IR sfmin
FORT  = joinpath(BUILD, "fortran_runs", "53155032", "out.TGLFEP")
BSIF  = get(ENV, "BSI", joinpath(@__DIR__, "batched_si_sfmin_nb32.txt"))              # idx IR sfmin

read_sfmin(p) = begin
    ir = Int[]; f = Float64[]
    for ln in eachline(p)
        t = split(strip(ln)); length(t) >= 3 || continue
        push!(ir, parse(Int, t[2])); push!(f, parse(Float64, t[3]))
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
            v = tryparse(Float64, s); v === nothing || push!(vals, v)
        end
    end
    vals
end

irg, grd = read_sfmin(GRIDF)
fort = read_fortran_sfmin(FORT)
irb, bsi = read_sfmin(BSIF)
@assert irb == irg "batched_si radii $(irb) != grid radii $(irg)"
dsf = abs.(bsi .- grd) ./ max.(abs.(grd), 1e-9)
mis = dsf .> 0.02
@printf("nb32: %d radii, %d mismatch >2%% vs grid golden (max %.1f%%)\n",
        length(irg), count(mis), 100*maximum(dsf))

default(legendfontsize=9, guidefontsize=11, tickfontsize=9, dpi=200)
p = plot(xlabel="radial grid index IR", ylabel="sfmin (critical factor)",
         title="SCAN_N=20 accuracy @ N_BASIS=32: :grid (golden) vs :batched_si",
         yscale=:log10, legend=:topleft, size=(950, 560), left_margin=5mm, bottom_margin=5mm)
length(fort) == length(irg) && plot!(p, irg, fort; label="Fortran (CPU)", color=:steelblue,
                                     marker=:circle, ms=6, lw=3)
plot!(p, irg, grd; label="Julia :grid (golden)", color=:gray55, ls=:dash, marker=:diamond, ms=5, lw=2)
plot!(p, irb, bsi; label="Julia :batched_si (fixed-shift hybrid)", color=:royalblue,
      marker=:utriangle, ms=6, lw=2)
any(mis) && scatter!(p, irb[mis], bsi[mis]; marker=:xcross, ms=10, color=:red, label="mismatch (>2%)")
for i in eachindex(irb)
    mis[i] && annotate!(p, irb[i], bsi[i]*1.28, text(@sprintf("%.0f%%", 100*dsf[i]), 7, :red, :left))
end
out = joinpath(@__DIR__, "batched_si_accuracy_nb32.png")
savefig(p, out); println("wrote ", out)
