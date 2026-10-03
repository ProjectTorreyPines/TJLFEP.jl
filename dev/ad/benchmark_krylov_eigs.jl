#!/usr/bin/env julia
# Phase B1: can a Krylov (Arnoldi) dominant-mode solver replace the full dense eigensolve
# (LAPACK geev / cuSOLVER Xgeev) in TJLF's hot path?
#
# Context (measured earlier): Xgeev is 98% of each GPU solve, cannot be batched, and is NOT
# concurrency-safe within a CUDA context; CPU geev is ~40-60 ms/solve at n~240. But TJLF only
# USES the few most-unstable modes (tjlf_LINEAR_SOLUTION sorts by real part and keeps
# NMODES; IBRANCH=0 keeps the most unstable electron- and ion-branch modes). So computing the
# full spectrum is mostly waste — IF a Krylov solver targeting :LR (largest real part = growth
# rate) converges reliably on REAL TJLF pencils.
#
# This script benchmarks on pencils harvested from the DIII-D nb=16 grid scan via
# TJLF_DUMP_PENCILS (see _maybe_dump_pencil in tjlf_eigensolver.jl):
#   ref    : gesv!(B,A) + geev! full spectrum (what TJLF does today on CPU)
#   krylov : M = lu(B) \ A once (O(n^3/3), cheap), then KrylovKit.eigsolve(:LR, nev)
# and reports, per pencil: wall times, Arnoldi matvec count, and the errors in the quantities
# TJLF consumes: gamma/freq of the top-NMODES modes and of the per-branch (ri sign) leaders.
# It also reports the real-part RANK of the ion/electron branch leaders in the full spectrum —
# this sets the nev needed for IBRANCH=0 correctness.
#
# Usage: julia --project=. -t N dev/ad/benchmark_krylov_eigs.jl
# Env: PENCILS (dir, default dev/ad/pencils_nb16), NEV (8), KRYLOVDIM (40), TOL (1e-9),
#      NMODES (4), EPS1 (1e-12: unstable-mode cutoff, epsilon1 in tjlf_LS)
using Pkg
Pkg.activate(normpath(@__DIR__, "..", ".."))
using LinearAlgebra, Serialization, Printf, Statistics
import LinearAlgebra.LAPACK: gesv!, geev!
import KrylovKit

const PDIR   = get(ENV, "PENCILS", normpath(@__DIR__, "pencils_nb16"))
const NEV    = parse(Int, get(ENV, "NEV", "8"))
const KDIM   = parse(Int, get(ENV, "KRYLOVDIM", "40"))
const TOL    = parse(Float64, get(ENV, "TOL", "1e-9"))
const NMODES = parse(Int, get(ENV, "NMODES", "4"))
const EPS1   = parse(Float64, get(ENV, "EPS1", "1e-12"))

BLAS.set_num_threads(1)

# TJLF keeps modes sorted by descending real part (growth rate). Return the top-k eigenvalues.
topk(vals, k) = sort(vals; by=real, rev=true)[1:min(k, length(vals))]

# IBRANCH=0 branch leaders: most unstable with ri>0 ("ion") and ri<=0 ("electron").
function branch_leaders(vals)
    ion = filter(v -> real(v) > EPS1 && imag(v) > 0, vals)
    ele = filter(v -> real(v) > EPS1 && imag(v) <= 0, vals)
    (ion = isempty(ion) ? nothing : ion[argmax(real.(ion))],
     ele = isempty(ele) ? nothing : ele[argmax(real.(ele))])
end

# Rank (1-based) of eigenvalue v in the full spectrum ordered by descending real part.
rank_of(v, vals) = findfirst(==(v), sort(vals; by=real, rev=true))

function run_ref(A, B)
    A2, B2 = copy(A), copy(B)
    t = @elapsed begin
        (A3, _, _) = gesv!(B2, A2)
        vals = geev!('N', 'N', A3)[1]
    end
    return vals, t
end

function run_krylov(A, B; nev=NEV, kdim=KDIM, tol=TOL)
    local vals, info, tlu, tar
    tlu = @elapsed M = lu(B) \ A          # dense M = B^-1 A, one LU + n triangular solves
    x0 = ones(ComplexF64, size(M, 1))
    tar = @elapsed vals, _, info = KrylovKit.eigsolve(M, x0, nev, :LR;
                                                      krylovdim=kdim, tol=tol, maxiter=100,
                                                      ishermitian=false, verbosity=0)
    return vals, info, tlu, tar
end

function main()
    files = sort(filter(f -> endswith(f, ".jls"), readdir(PDIR; join=true)))
    isempty(files) && error("no pencils in $PDIR — run a grid radius with TJLF_DUMP_PENCILS set")
    println("pencils: ", length(files), " from ", PDIR)
    p1 = Serialization.deserialize(files[1])
    println("matrix size n = ", size(p1.A, 1), "   NEV=", NEV, " KRYLOVDIM=", KDIM, " TOL=", TOL)

    t_ref = Float64[]; t_lu = Float64[]; t_ar = Float64[]; nops = Int[]
    err_top = Float64[]           # worst |Δλ| over top-NMODES modes (matched by order)
    err_ion = Float64[]; err_ele = Float64[]
    rank_ion = Int[]; rank_ele = Int[]
    n_unconv = 0; n_topmiss = 0

    for (i, f) in enumerate(files)
        p = Serialization.deserialize(f)
        vals_ref, tr = run_ref(p.A, p.B)
        vals_kr, info, tl, ta = run_krylov(p.A, p.B)
        push!(t_ref, tr); push!(t_lu, tl); push!(t_ar, ta); push!(nops, info.numops)
        info.converged < min(NEV, length(vals_kr)) && (n_unconv += 1)

        # top-NMODES comparison (the gamma/freq TJLF returns for IBRANCH=-1)
        tr_ref = topk(vals_ref, NMODES)
        tr_kr  = topk(vals_kr,  NMODES)
        if length(tr_kr) < length(tr_ref)
            n_topmiss += 1
        else
            push!(err_top, maximum(abs.(tr_ref .- tr_kr[1:length(tr_ref)])))
        end

        # IBRANCH=0 branch leaders + their spectral rank
        bl = branch_leaders(vals_ref)
        if bl.ion !== nothing
            push!(rank_ion, rank_of(bl.ion, vals_ref))
            blk = branch_leaders(vals_kr)
            push!(err_ion, blk.ion === nothing ? NaN : abs(bl.ion - blk.ion))
        end
        if bl.ele !== nothing
            push!(rank_ele, rank_of(bl.ele, vals_ref))
            blk = branch_leaders(vals_kr)
            push!(err_ele, blk.ele === nothing ? NaN : abs(bl.ele - blk.ele))
        end
    end

    q(v, p) = isempty(v) ? NaN : quantile(v, p)
    println("\n================ timings (ms/pencil) ================")
    @printf("  full geev (ref)     : median %7.2f   p90 %7.2f\n", 1e3*q(t_ref,0.5), 1e3*q(t_ref,0.9))
    @printf("  M=B\\A LU form       : median %7.2f   p90 %7.2f\n", 1e3*q(t_lu,0.5),  1e3*q(t_lu,0.9))
    @printf("  Arnoldi eigsolve    : median %7.2f   p90 %7.2f\n", 1e3*q(t_ar,0.5),  1e3*q(t_ar,0.9))
    tk = t_lu .+ t_ar
    @printf("  krylov total        : median %7.2f   p90 %7.2f   speedup vs geev: %.1fx (median)\n",
            1e3*q(tk,0.5), 1e3*q(tk,0.9), q(t_ref,0.5)/q(tk,0.5))
    @printf("  Arnoldi matvecs     : median %d   p90 %d   max %d\n",
            round(Int,q(Float64.(nops),0.5)), round(Int,q(Float64.(nops),0.9)), maximum(nops))

    println("\n================ accuracy vs full geev ================")
    @printf("  top-%d modes |dLambda| : median %.3g   max %.3g   (missed sets: %d/%d)\n",
            NMODES, q(err_top,0.5), isempty(err_top) ? NaN : maximum(err_top), n_topmiss, length(files))
    ei = filter(!isnan, err_ion); ee = filter(!isnan, err_ele)
    @printf("  ion-branch leader     : found %d/%d   |dLambda| median %.3g max %.3g\n",
            length(ei), length(err_ion), q(ei,0.5), isempty(ei) ? NaN : maximum(ei))
    @printf("  ele-branch leader     : found %d/%d   |dLambda| median %.3g max %.3g\n",
            length(ee), length(err_ele), q(ee,0.5), isempty(ee) ? NaN : maximum(ee))
    @printf("  unconverged eigsolves : %d/%d\n", n_unconv, length(files))

    println("\n================ branch-leader spectral rank (sets required NEV) ================")
    for (name, r) in (("ion", rank_ion), ("ele", rank_ele))
        isempty(r) && continue
        @printf("  %s: median %d   p90 %d   max %d   (pencils with a leader: %d)\n",
                name, round(Int,q(Float64.(r),0.5)), round(Int,q(Float64.(r),0.9)), maximum(r), length(r))
    end
end

main()
