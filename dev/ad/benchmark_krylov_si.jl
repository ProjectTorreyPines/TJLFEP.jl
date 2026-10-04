#!/usr/bin/env julia
# Phase B1b: generalized shift-invert Arnoldi for TJLF's dominant modes.
#
# Naive KrylovKit :LR on M=B^-1 A FAILED on real pencils (48/48 unconverged, 0.9x speed): the
# spectrum bulk has |lambda|~40 (damped) while the physically relevant unstable modes sit in a
# dense cluster |lambda| <~ 1.5 near the origin, so power-type iterations lock onto the wrong end.
#
# Shift-invert fixes this: for a shift sigma, the operator
#     S x = (A - sigma B)^{-1} (B x)      (one LU of A - sigma B, then gemv + 2 triangular solves)
# has eigenvalues mu = 1/(lambda - sigma), so modes NEAREST sigma become DOMINANT (:LM), which
# Arnoldi converges to quickly. We never form M = B^{-1}A (saves the dense B\A solve).
# lambda = sigma + 1/mu.
#
# Strategy tested here ("SI union"): a few shifts placed to cover the unstable window
# (gamma>0, |freq| <~ 1.5 from the spectra inspected): sigma in {0.05, 0.05+0.5im, 0.05-0.5im,
# 0.05+1.0im, ...}; take the union of Ritz values, dedup, and compare the TJLF-consumed outputs
# (top-NMODES by Re, IBRANCH=0 branch leaders) against full geev. Also reports how many exact
# eigenvalues lie in the candidate disk (context for required nev).
#
# Usage: julia --project=. dev/ad/benchmark_krylov_si.jl
# Env: PENCILS, NEV (10), KRYLOVDIM (30), TOL (1e-8), NMODES (4), EPS1 (1e-12),
#      SHIFTS ("0.05,0.05+0.5im,0.05-0.5im,0.05+1.0im,0.05-1.0im")
using Pkg
Pkg.activate(normpath(@__DIR__, "..", ".."))
using LinearAlgebra, Serialization, Printf, Statistics
import LinearAlgebra.LAPACK: gesv!, geev!
import KrylovKit

const PDIR   = get(ENV, "PENCILS", normpath(@__DIR__, "pencils_nb16"))
const NEV    = parse(Int, get(ENV, "NEV", "10"))
const KDIM   = parse(Int, get(ENV, "KRYLOVDIM", "30"))
const TOL    = parse(Float64, get(ENV, "TOL", "1e-8"))
const NMODES = parse(Int, get(ENV, "NMODES", "4"))
const EPS1   = parse(Float64, get(ENV, "EPS1", "1e-12"))
const SHIFTS = [parse(ComplexF64, s) for s in
                split(get(ENV, "SHIFTS", "0.05,0.05+0.5im,0.05-0.5im,0.05+1.0im,0.05-1.0im"), ",")]

BLAS.set_num_threads(1)

topk(vals, k) = sort(vals; by=real, rev=true)[1:min(k, length(vals))]

function branch_leaders(vals)
    ion = filter(v -> real(v) > EPS1 && imag(v) > 0, vals)
    ele = filter(v -> real(v) > EPS1 && imag(v) <= 0, vals)
    (ion = isempty(ion) ? nothing : ion[argmax(real.(ion))],
     ele = isempty(ele) ? nothing : ele[argmax(real.(ele))])
end

function run_ref(A, B)
    A2, B2 = copy(A), copy(B)
    t = @elapsed begin
        (A3, _, _) = gesv!(B2, A2)
        vals = geev!('N', 'N', A3)[1]
    end
    return vals, t
end

# Shift-invert Arnoldi at one shift: returns back-transformed Ritz values + stats.
function si_arnoldi(A, B, sigma; nev=NEV, kdim=KDIM, tol=TOL)
    n = size(A, 1)
    F = lu(A - sigma * B)
    op = x -> F \ (B * x)
    x0 = ones(ComplexF64, n) ./ sqrt(n)
    mus, _, info = KrylovKit.eigsolve(op, x0, nev, :LM;
                                      krylovdim=kdim, tol=tol, maxiter=30,
                                      ishermitian=false, verbosity=0)
    lams = [sigma + 1 / mu for mu in mus if abs(mu) > 1e-300]
    return lams, info.converged, info.numops
end

function run_si_union(A, B)
    lams = ComplexF64[]
    nconv = 0; nops = 0
    t = @elapsed for sigma in SHIFTS
        l, c, o = si_arnoldi(A, B, sigma)
        append!(lams, l); nconv += c; nops += o
    end
    # dedup (different shifts re-find the same modes)
    sort!(lams; by=x->(real(x), imag(x)))
    ded = ComplexF64[]
    for l in lams
        if isempty(ded) || abs(l - ded[end]) > 1e-8 * max(1.0, abs(l))
            push!(ded, l)
        end
    end
    return ded, t, nconv, nops
end

function main()
    files = sort(filter(f -> endswith(f, ".jls"), readdir(PDIR; join=true)))
    isempty(files) && error("no pencils in $PDIR")
    p1 = Serialization.deserialize(files[1])
    println("pencils: ", length(files), "   n = ", size(p1.A, 1))
    println("shifts: ", SHIFTS, "   NEV=", NEV, " KRYLOVDIM=", KDIM, " TOL=", TOL)

    t_ref = Float64[]; t_si = Float64[]; nops_all = Int[]
    err_top = Float64[]; err_ion = Float64[]; err_ele = Float64[]
    n_cluster = Int[]   # exact modes with Re > EPS1 (unstable count)
    n_ionmiss = 0; n_elemiss = 0; n_ion = 0; n_ele = 0
    # gamma errors in the physically-used scalar (growth rate of leaders)
    gerr_ion = Float64[]; gerr_ele = Float64[]

    for f in files
        p = Serialization.deserialize(f)
        vals_ref, tr = run_ref(p.A, p.B)
        cand, ts, _, no = run_si_union(p.A, p.B)
        push!(t_ref, tr); push!(t_si, ts); push!(nops_all, no)
        push!(n_cluster, count(v -> real(v) > EPS1, vals_ref))

        tr_ref = topk(vals_ref, NMODES); tr_si = topk(cand, NMODES)
        m = min(length(tr_ref), length(tr_si))
        m > 0 && push!(err_top, maximum(abs.(tr_ref[1:m] .- tr_si[1:m])))

        bl  = branch_leaders(vals_ref)
        bls = branch_leaders(cand)
        if bl.ion !== nothing
            n_ion += 1
            if bls.ion === nothing; n_ionmiss += 1
            else push!(err_ion, abs(bl.ion - bls.ion)); push!(gerr_ion, abs(real(bl.ion) - real(bls.ion))) end
        end
        if bl.ele !== nothing
            n_ele += 1
            if bls.ele === nothing; n_elemiss += 1
            else push!(err_ele, abs(bl.ele - bls.ele)); push!(gerr_ele, abs(real(bl.ele) - real(bls.ele))) end
        end
    end

    q(v, p) = isempty(v) ? NaN : quantile(v, p)
    println("\n================ timings (ms/pencil, 1 BLAS thread) ================")
    @printf("  full geev (ref)   : median %7.1f   p90 %7.1f\n", 1e3*q(t_ref,0.5), 1e3*q(t_ref,0.9))
    @printf("  SI union (%d LUs) : median %7.1f   p90 %7.1f   speedup %.1fx (median)\n",
            length(SHIFTS), 1e3*q(t_si,0.5), 1e3*q(t_si,0.9), q(t_ref,0.5)/q(t_si,0.5))
    @printf("  total matvecs     : median %d   max %d\n",
            round(Int, q(Float64.(nops_all),0.5)), maximum(nops_all))

    println("\n================ accuracy vs full geev ================")
    @printf("  unstable modes per pencil (Re>eps): median %d   max %d\n",
            round(Int, q(Float64.(n_cluster),0.5)), maximum(n_cluster))
    @printf("  top-%d |dLambda|  : median %.3g   p90 %.3g   max %.3g\n",
            NMODES, q(err_top,0.5), q(err_top,0.9), isempty(err_top) ? NaN : maximum(err_top))
    @printf("  ion leader  : missed %d/%d   |dLambda| median %.3g max %.3g   |dgamma| max %.3g\n",
            n_ionmiss, n_ion, q(err_ion,0.5), isempty(err_ion) ? NaN : maximum(err_ion),
            isempty(gerr_ion) ? NaN : maximum(gerr_ion))
    @printf("  ele leader  : missed %d/%d   |dLambda| median %.3g max %.3g   |dgamma| max %.3g\n",
            n_elemiss, n_ele, q(err_ele,0.5), isempty(err_ele) ? NaN : maximum(err_ele),
            isempty(gerr_ele) ? NaN : maximum(gerr_ele))
end

main()
