#!/usr/bin/env julia
# Overlay Fortran vs Julia α profiles for nb6 SCAN_N=20 (subsample at IR_EXP).

using Pkg
Pkg.activate(normpath(@__DIR__, "..", ".."))

using TJLFEP
using Plots

const BUILD = normpath(@__DIR__, "..")
const CASE = normpath(@__DIR__, "..", "..", "examples", "DIIID_202017C42_500ms_v3.1")
const FORTRAN_DIR = get(ENV, "FORTRAN_DIR", joinpath(BUILD, "fortran_runs", "debug_nb6_scan20_10n_53171364"))
const JULIA_DIR = get(ENV, "JULIA_DIR", joinpath(BUILD, "debug_out_nb6_scan20_53171385_dist"))
const FILE_DIR = get(ENV, "FILE_DIR", joinpath(BUILD, "fileInput_nb6_scan20_53171385"))
const OUTDIR = get(ENV, "PLOT_OUTDIR", joinpath(BUILD, "compare_nb6_scan20_plots"))
const FORTRAN_LABEL = get(ENV, "FORTRAN_LABEL", "Fortran")
const JULIA_LABEL = get(ENV, "JULIA_LABEL", "Julia")
const PLOT_TITLE = get(ENV, "PLOT_TITLE", "N_BASIS=6, SCAN_N=20 (10 nodes)")

function read_sfmin(dir::String)
    p = joinpath(dir, "out.TGLFEP")
    if isfile(p)
        lines = readlines(p)
        i = findfirst(l -> strip(l) == "SFmin", lines)
        out = Float64[]
        for line in lines[i+1:end]
            s = strip(line)
            isempty(s) && break
            startswith(s, "The ") && break
            if startswith(s, "[")
                m = match(r"\[(.*)\]", s)
                return parse.(Float64, strip.(split(m.captures[1], ",")))
            end
            x = tryparse(Float64, split(s, " ")[1])
            x === nothing && break
            push!(out, x)
        end
        return out
    end
    p = joinpath(dir, "sfmin_scan.txt")
    isfile(p) || error("need out.TGLFEP or sfmin_scan.txt in $dir")
    return [parse(Float64, split(strip(line))[3]) for line in readlines(p) if !isempty(strip(line))]
end

function read_alpha_file(path::String)
    lines = readlines(path)
    body = strip(lines[2])
    if startswith(body, "[")
        arr = split(strip(body, ['[', ']']), ",")
        return strip(lines[1]), parse.(Float64, strip.(arr))
    end
    vals = Float64[]
    for line in lines[2:end]
        s = strip(line)
        isempty(s) && continue
        x = tryparse(Float64, s)
        x === nothing && continue
        push!(vals, x)
    end
    return strip(lines[1]), vals
end

prof = read_input_profile(joinpath(CASE, "dump.profile"))
rho_full = prof.RMIN ./ prof.RMIN[end]

_, ir_exp = readMTGLF(joinpath(FILE_DIR, "input.MTGLF"))
ir_exp = Int.(ir_exp)
@assert length(ir_exp) == 20

rho = rho_full[ir_exp]

_, dndr_f = read_alpha_file(joinpath(FORTRAN_DIR, "alpha_dndr_crit.input"))
_, dpdr_f = read_alpha_file(joinpath(FORTRAN_DIR, "alpha_dpdr_crit.input"))
_, dndr_j = read_alpha_file(joinpath(JULIA_DIR, "alpha_dndr_crit.input"))
_, dpdr_j = read_alpha_file(joinpath(JULIA_DIR, "alpha_dpdr_crit.input"))

@assert length(dndr_f) >= maximum(ir_exp)
@assert length(dndr_j) >= maximum(ir_exp)

dndr_f_s = dndr_f[ir_exp]
dpdr_f_s = dpdr_f[ir_exp]
dndr_j_s = dndr_j[ir_exp]
dpdr_j_s = dpdr_j[ir_exp]

mkpath(OUTDIR)

pairs = [
    (dndr_f_s, dndr_j_s, "dn/dr_crit", "10¹⁹ m⁻⁴", "alpha_dndr_crit"),
    (dpdr_f_s, dpdr_j_s, "dp/dr_crit", "10 kPa/m", "alpha_dpdr_crit"),
]

for (vf, vj, ylab_short, units, stem) in pairs
    p = plot(
        rho, vf;
        label = FORTRAN_LABEL,
        lw = 2,
        marker = :circle,
        color = :blue,
        xlabel = "ρ (r / r_edge)",
        ylabel = "$ylab_short ($units)",
        title = "$PLOT_TITLE — $(replace(stem, "_" => " "))",
        legend = :best,
        size = (900, 550),
    )
    plot!(p, rho, vj; label = JULIA_LABEL, lw = 2, marker = :diamond, color = :red, linestyle = :dash)
    savefig(p, joinpath(OUTDIR, "$(stem)_compare.png"))
    rel = abs.(vj .- vf) ./ max.(abs.(vf), 1e-30)
    println(stem, ": max rel err = ", maximum(rel), " mean = ", sum(rel) / length(rel))
    println("  Wrote ", joinpath(OUTDIR, "$(stem)_compare.png"))
end

p = plot(layout = (2, 1), size = (900, 800), legend = :best)
plot!(p[1], rho, dndr_f_s; label = FORTRAN_LABEL, lw = 2, marker = :circle, color = :blue)
plot!(p[1], rho, dndr_j_s; label = JULIA_LABEL, lw = 2, marker = :diamond, color = :red, linestyle = :dash)
plot!(p[1]; xlabel = "ρ", ylabel = "dn/dr_crit (10¹⁹ m⁻⁴)", title = "Critical density gradient")

plot!(p[2], rho, dpdr_f_s; label = FORTRAN_LABEL, lw = 2, marker = :circle, color = :blue)
plot!(p[2], rho, dpdr_j_s; label = JULIA_LABEL, lw = 2, marker = :diamond, color = :red, linestyle = :dash)
plot!(p[2]; xlabel = "ρ", ylabel = "dp/dr_crit (10 kPa/m)", title = "Critical pressure gradient")

combined = joinpath(OUTDIR, "alpha_crit_grads_compare.png")
savefig(p, combined)
println("Wrote ", combined)

sf_f = read_sfmin(FORTRAN_DIR)
sf_j = read_sfmin(JULIA_DIR)
n = min(length(sf_f), length(sf_j), length(ir_exp))
p_sf = plot(
    rho[1:n], sf_f[1:n];
    label = FORTRAN_LABEL,
    lw = 2,
    marker = :circle,
    color = :blue,
    xlabel = "ρ",
    ylabel = "SFmin",
    title = "$PLOT_TITLE — scale factor",
    legend = :best,
    size = (900, 550),
)
plot!(p_sf, rho[1:n], sf_j[1:n]; label = JULIA_LABEL, lw = 2, marker = :diamond, color = :red, linestyle = :dash)
sf_path = joinpath(OUTDIR, "sfmin_compare.png")
savefig(p_sf, sf_path)
rel_sf = abs.(sf_j[1:n] .- sf_f[1:n]) ./ max.(abs.(sf_f[1:n]), 1e-30)
println("SFmin: max rel err = ", maximum(rel_sf), " mean = ", sum(rel_sf) / n)
println("Wrote ", sf_path)
