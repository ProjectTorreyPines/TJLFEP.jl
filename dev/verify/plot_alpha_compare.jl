#!/usr/bin/env julia
# Overlay Fortran vs Julia α(dn/dr) and α(dp/dr) from saved alpha_*.input files.

using Pkg
Pkg.activate(normpath(@__DIR__, "..", ".."))

using TJLFEP
using Plots

const BUILD = normpath(@__DIR__, "..")
const CASE = normpath(@__DIR__, "..", "..", "examples", "DIIID_202017C42_500ms_v3.1")
const FORTRAN_DIR = get(ENV, "FORTRAN_DIR", joinpath(BUILD, "fortran_runs", "debug_nb6_local"))
const JULIA_DIR = get(ENV, "JULIA_DIR", joinpath(BUILD, "debug_out_nb6_local"))
const OUTDIR = get(ENV, "PLOT_OUTDIR", joinpath(BUILD, "compare_alpha_plots"))
const FORTRAN_LABEL = get(ENV, "FORTRAN_LABEL", "Fortran")
const JULIA_LABEL = get(ENV, "JULIA_LABEL", "Julia")
const PLOT_TITLE = get(ENV, "PLOT_TITLE", "N_BASIS=6, SCAN_N=1")

function read_alpha_file(path::String)
    lines = readlines(path)
    header = strip(lines[1])
    body = strip(lines[2])
    if startswith(body, "[")
        arr = split(strip(body, ['[', ']']), ",")
        vals = parse.(Float64, strip.(arr))
    else
        vals = Float64[]
        for line in lines[2:end]
            s = strip(line)
            isempty(s) && continue
            x = tryparse(Float64, s)
            x === nothing && continue
            push!(vals, x)
        end
    end
    return header, vals
end

prof = read_input_profile(joinpath(CASE, "dump.profile"))
rho = prof.RMIN ./ prof.RMIN[end]

mkpath(OUTDIR)

_, dndr_f = read_alpha_file(joinpath(FORTRAN_DIR, "alpha_dndr_crit.input"))
_, dpdr_f = read_alpha_file(joinpath(FORTRAN_DIR, "alpha_dpdr_crit.input"))
_, dndr_j = read_alpha_file(joinpath(JULIA_DIR, "alpha_dndr_crit.input"))
_, dpdr_j = read_alpha_file(joinpath(JULIA_DIR, "alpha_dpdr_crit.input"))

pairs = [
    (dndr_f, dndr_j, "dn/dr_crit", "10¹⁹ m⁻⁴", "alpha_dndr_crit"),
    (dpdr_f, dpdr_j, "dp/dr_crit", "10 kPa/m", "alpha_dpdr_crit"),
]

for (vf, vj, ylab_short, units, stem) in pairs
    nf = min(length(rho), length(vf))
    nj = min(length(rho), length(vj))
    p = plot(
        rho[1:nf], vf[1:nf];
        label = FORTRAN_LABEL,
        lw = 2,
        color = :blue,
        xlabel = "ρ (r / r_edge)",
        ylabel = "$ylab_short ($units)",
        title = "$PLOT_TITLE — $(replace(stem, "_" => " "))",
        legend = :best,
        size = (900, 550),
    )
    plot!(p, rho[1:nj], vj[1:nj]; label = JULIA_LABEL, lw = 2, color = :red, linestyle = :dash)
    savefig(p, joinpath(OUTDIR, "$(stem)_compare.png"))
    println("Wrote ", joinpath(OUTDIR, "$(stem)_compare.png"))
end

p = plot(layout = (2, 1), size = (900, 800), legend = :best)
nf = min(length(rho), length(dndr_f), length(dndr_j))
plot!(p[1], rho[1:nf], dndr_f[1:nf]; label = FORTRAN_LABEL, lw = 2, color = :blue)
plot!(p[1], rho[1:nf], dndr_j[1:nf]; label = JULIA_LABEL, lw = 2, color = :red, linestyle = :dash)
plot!(p[1]; xlabel = "ρ", ylabel = "dn/dr_crit (10¹⁹ m⁻⁴)", title = "Critical density gradient")

np = min(length(rho), length(dpdr_f), length(dpdr_j))
plot!(p[2], rho[1:np], dpdr_f[1:np]; label = FORTRAN_LABEL, lw = 2, color = :blue)
plot!(p[2], rho[1:np], dpdr_j[1:np]; label = JULIA_LABEL, lw = 2, color = :red, linestyle = :dash)
plot!(p[2]; xlabel = "ρ", ylabel = "dp/dr_crit (10 kPa/m)", title = "Critical pressure gradient")

combined = joinpath(OUTDIR, "alpha_crit_grads_compare.png")
savefig(p, combined)
println("Wrote ", combined)
