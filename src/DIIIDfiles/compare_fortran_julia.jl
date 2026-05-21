#!/usr/bin/env julia
# Compare Julia validation outputs to Fortran reference (202017C42_500ms_v3.1).
# Usage:
#   julia --project=../.. compare_fortran_julia.jl [julia_outdir]
# Default julia_outdir: most recent GPU_* or CPU_* directory in pwd.

using Printf

const FORTRAN_REF_DIR = get(ENV, "FORTRAN_REF_DIR", joinpath(@__DIR__, "202017C42_500ms_v3.1"))

function latest_julia_outdir()
    candidates = filter(d -> startswith(d, "GPU_") || startswith(d, "CPU_"), readdir(@__DIR__))
    isempty(candidates) && return nothing
    sort!(candidates, by = d -> stat(joinpath(@__DIR__, d)).mtime, rev = true)
    return joinpath(@__DIR__, candidates[1])
end

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
    println("\n=== SFmin (out.TGLFEP or sfmin_scan.txt) ===")
    sf_f = read_sfmin(f_dir)
    sf_j = read_sfmin(j_dir)
    n = min(length(sf_f), length(sf_j))
    println(@sprintf("%-4s %-12s %-12s %-10s", "i", "F", "J", "rel_err"))
    max_rel = 0.0
    for i in 1:n
        rel = abs(sf_j[i] - sf_f[i]) / max(abs(sf_f[i]), 1e-30)
        max_rel = max(max_rel, rel)
        @printf("%3d  F=%.6g  J=%.6g  rel=%.4g\n", i, sf_f[i], sf_j[i], rel)
    end
    println(@sprintf("Compared %d radii; max relative |SF_J-SF_F|/|SF_F| = %.4g", n, max_rel))
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
