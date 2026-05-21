#!/usr/bin/env julia
# Validate profile_from_gacode against Fortran dump.profile on the same case.
#
#   julia --startup-file=no --project=.. build/compare_gacode_profile_vs_dump.jl

using Pkg
Pkg.activate(normpath(@__DIR__, ".."))

using TJLFEP
using TJLFEP: profile_from_gacode, read_input_profile
using Printf
using Statistics

const CASE_DIR = get(ENV, "CASE_DIR",
    normpath(@__DIR__, "..", "src", "DIIIDfiles", "202017C42_500ms_v3.1"))
const GACODE_FILE = get(ENV, "GACODE_FILE", joinpath(CASE_DIR, "dump.gacode"))
const PROFILE_REF = joinpath(CASE_DIR, "dump.profile")
const IS_EP = parse(Int, get(ENV, "IS_EP", "2"))
const N_ION = parse(Int, get(ENV, "N_ION", "2"))

function rel_err(a, b)
    denom = max(abs(a), abs(b), 1e-30)
    return abs(b - a) / denom
end

function compare_vec(name, ref::AbstractVector, test::AbstractVector)
    errs = rel_err.(ref, test)
    @printf("  %-18s max=%.3e mean=%.3e\n", name, maximum(errs), mean(errs))
    return maximum(errs)
end

function compare_mat(name, ref::AbstractMatrix, test::AbstractMatrix)
    errs = rel_err.(ref, test)
    @printf("  %-18s max=%.3e mean=%.3e\n", name, maximum(errs), mean(errs))
    return maximum(errs)
end

prof_ref = read_input_profile(PROFILE_REF)
prof_test = profile_from_gacode(GACODE_FILE; is_ep=IS_EP, tglfep_nion=N_ION)

@assert prof_ref.NR == prof_test.NR
@assert prof_ref.NS == prof_test.NS

println("=== profile_from_gacode vs dump.profile ===")
println("GACODE: ", GACODE_FILE)
println("REF:    ", PROFILE_REF)
println("NR=$(prof_ref.NR) NS=$(prof_ref.NS) is_ep=$IS_EP")

println("\n--- species ---")
compare_mat("AS", prof_ref.AS, prof_test.AS)
compare_mat("TAUS", prof_ref.TAUS, prof_test.TAUS)
compare_mat("RLNS", prof_ref.RLNS, prof_test.RLNS)
compare_mat("RLTS", prof_ref.RLTS, prof_test.RLTS)

println("\n--- geometry ---")
compare_vec("RMIN", prof_ref.RMIN, prof_test.RMIN)
compare_vec("RMAJ", prof_ref.RMAJ, prof_test.RMAJ)
compare_vec("Q", prof_ref.Q, prof_test.Q)
compare_vec("SHEAR", prof_ref.SHEAR, prof_test.SHEAR)
compare_vec("Q_PRIME", prof_ref.Q_PRIME, prof_test.Q_PRIME)
compare_vec("P_PRIME", prof_ref.P_PRIME, prof_test.P_PRIME)
compare_vec("KAPPA", prof_ref.KAPPA, prof_test.KAPPA)
compare_vec("S_KAPPA", prof_ref.S_KAPPA, prof_test.S_KAPPA)
compare_vec("DELTA", prof_ref.DELTA, prof_test.DELTA)
compare_vec("BETAE", prof_ref.BETAE, prof_test.BETAE)
compare_vec("ZEFF", prof_ref.ZEFF, prof_test.ZEFF)

println("\n--- sample RLNS_EP at scan-like indices ---")
for ir in [2, 22, 48, 74, 101]
    @printf("  ir=%3d  RLNS_EP ref=%.5f test=%.5f  rel=%.3e\n",
        ir, prof_ref.RLNS[ir, 3], prof_test.RLNS[ir, 3],
        rel_err(prof_ref.RLNS[ir, 3], prof_test.RLNS[ir, 3]))
end

println("\n=== done ===")
