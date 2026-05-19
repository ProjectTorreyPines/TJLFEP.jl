#!/usr/bin/env julia
# Decompose critical-gradient differences: SFmin vs profile factors at IR_EXP.
# Usage: julia --project=../.. diagnose_crit_grad.jl <julia_outdir> [scan_index]

using Printf

const FORTRAN_DIR = joinpath(@__DIR__, "202017C42_500ms_v3.1")
const RHO = [0.01, 0.06, 0.11, 0.16, 0.21, 0.27, 0.32, 0.37, 0.42, 0.47,
    0.53, 0.58, 0.63, 0.68, 0.73, 0.79, 0.84, 0.89, 0.94, 1.0]

function parse_sfmin(path)
    for line in reverse(readlines(path))
        parts = split(strip(line))
        isempty(parts) && continue
        x = tryparse(Float64, parts[1])
        x !== nothing && return x
    end
    return NaN
end

function main()
    length(ARGS) >= 1 || error("usage: diagnose_crit_grad.jl <julia_outdir> [iscan]")
    j_dir = abspath(ARGS[1])
    iscan = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 4
    f_files = sort(filter(f -> startswith(f, "out.scalefactor_r"), readdir(FORTRAN_DIR)))
    j_files = sort(filter(f -> startswith(f, "out.scalefactor_r"), readdir(j_dir)))
    ff = joinpath(FORTRAN_DIR, f_files[iscan])
    jf = joinpath(j_dir, j_files[iscan])
    sf_f = parse_sfmin(ff)
    sf_j = parse_sfmin(jf)
    @printf("scan i=%d  rho=%.2f\n", iscan, RHO[iscan])
    @printf("  SFmin Fortran = %.6g\n", sf_f)
    @printf("  SFmin Julia   = %.6g\n", sf_j)
    @printf("  ratio SF_J/SF_F = %.6g\n", isfinite(sf_f) && sf_f != 0 ? sf_j / sf_f : NaN)
    IRS, NR, SCAN_N = 2, 101, 20
    ir_f = IRS + floor(Int, (iscan - 1) * (NR - IRS) / (SCAN_N - 1))
    println("  IR_EXP (Fortran formula) = ", ir_f)
    println("\nIf dndr_J/dndr_F ≈ (SF_J/SF_F)*(ni*dlnn)_J/(ni*dlnn)_F, residual is TJLF physics; else check profiles.")
end

main()
