#!/usr/bin/env julia
# Compare Julia validation outputs to Fortran reference (202017C42_500ms_v3.1).
# Usage:
#   julia --project=../.. compare_fortran_julia.jl [julia_outdir]
# Default julia_outdir: most recent GPU_* or CPU_* directory in pwd.

using Printf

# Archived reference outputs (read-only; do not run jobs in this directory).
const FORTRAN_REF_DIR = get(ENV, "FORTRAN_REF_DIR", joinpath(@__DIR__, "202017C42_500ms_v3.1"))
const RHO_SCAN = [0.01, 0.06, 0.11, 0.16, 0.21, 0.27, 0.32, 0.37, 0.42, 0.47,
    0.53, 0.58, 0.63, 0.68, 0.73, 0.79, 0.84, 0.89, 0.94, 1.0]

function latest_julia_outdir()
    candidates = filter(d -> startswith(d, "GPU_") || startswith(d, "CPU_"), readdir(@__DIR__))
    isempty(candidates) && return nothing
    sort!(candidates, by = d -> stat(joinpath(@__DIR__, d)).mtime, rev = true)
    return joinpath(@__DIR__, candidates[1])
end

"""Parse final converged scale factor from out.scalefactor_r### (last non-empty data line)."""
function parse_sfmin(path::String)
    isfile(path) || return missing
    lines = readlines(path)
    for line in reverse(lines)
        s = strip(line)
        isempty(s) && continue
        startswith(s, "factor") && continue
        startswith(s, "---") && continue
        startswith(s, "omega") && continue
        startswith(s, "Frequencies") && continue
        parts = split(s)
        length(parts) < 1 && continue
        x = tryparse(Float64, parts[1])
        x !== nothing && return x
    end
    return missing
end

function parse_crit_file(path::String)
    isfile(path) || return Float64[]
    vals = Float64[]
    for line in readlines(path)
        s = strip(line)
        isempty(s) && continue
        x = tryparse(Float64, s)
        x === nothing && continue
        push!(vals, x)
    end
    return vals
end

function compare_sfmin(f_dir::String, j_dir::String)
    println("\n=== out.scalefactor_r### (final factor per radius) ===")
    println(@sprintf("%-4s %-22s %-12s %-12s %-10s", "k", "file", "F", "J", "rel_err"))
    max_rel = 0.0
    n = 0
    f_files = sort(filter(f -> startswith(f, "out.scalefactor_r"), readdir(f_dir)))
    j_files = sort(filter(f -> startswith(f, "out.scalefactor_r"), readdir(j_dir)))
    ncomp = min(length(f_files), length(j_files))
    for k in 1:ncomp
        ff = joinpath(f_dir, f_files[k])
        jf = joinpath(j_dir, j_files[k])
        sf_f = parse_sfmin(ff)
        sf_j = parse_sfmin(jf)
        rel = (sf_f, sf_j) isa Tuple{Float64, Float64} ? abs(sf_j - sf_f) / max(abs(sf_f), 1e-30) : NaN
        if sf_f !== missing && sf_j !== missing
            max_rel = max(max_rel, rel)
            n += 1
        end
        @printf("%3d  %-20s  F=%.6g  J=%.6g  rel=%.4g\n", k, f_files[k], something(sf_f, NaN), something(sf_j, NaN), rel)
    end
    println(@sprintf("Compared %d radius files; max relative |SF_J-SF_F|/|SF_F| = %.4g", n, max_rel))
    return max_rel
end

function compare_crit(f_dir::String, j_dir::String, name::String)
    f_path = joinpath(f_dir, name)
    j_path = joinpath(j_dir, name)
    vf = parse_crit_file(f_path)
    vj = parse_crit_file(j_path)
    println("\n=== $name (length F=$(length(vf)), J=$(length(vj))) ===")
    n = min(length(vf), length(vj))
    n == 0 && (println("  missing or empty"); return)
    rel = [abs(vj[i] - vf[i]) / max(abs(vf[i]), 1e-30) for i in 1:n]
    @printf("  max rel err = %.4g  mean rel err = %.4g  at i=%d\n", maximum(rel), sum(rel) / n, argmax(rel))
end

function main()
    j_dir = length(ARGS) >= 1 ? ARGS[1] : latest_julia_outdir()
    if j_dir === nothing || !isdir(j_dir)
        println("No Julia output directory found. Pass path as first argument.")
        println("Fortran reference: ", FORTRAN_REF_DIR)
        exit(1)
    end
    j_dir = abspath(j_dir)
    f_dir = get(ENV, "FORTRAN_DIR", FORTRAN_REF_DIR)
    println("Fortran reference: ", f_dir)
    if abspath(f_dir) == abspath(FORTRAN_REF_DIR)
        println("(archived ref in 202017C42_500ms_v3.1; set FORTRAN_DIR=.../build/fortran_runs/<job> for a fresh Fortran run)")
    end
    println("Julia output:      ", j_dir)
    compare_sfmin(f_dir, j_dir)
    compare_crit(f_dir, j_dir, "alpha_dndr_crit.input")
    compare_crit(f_dir, j_dir, "alpha_dpdr_crit.input")
end

main()
