#!/usr/bin/env julia
# Phase B1c: BATCHED shift-invert subspace iteration on GPU — the concurrency-safe way to fill
# the A100 that Xgeev cannot provide.
#
# Idea: TJLF only consumes the few most-unstable modes, which live in a small cluster
# (gamma in (0, ~0.2], |freq| <~ 1.5) while the spectrum bulk is strongly damped (|lambda|~40).
# CPU prototyping (benchmark_krylov_si.jl) showed shift-invert Arnoldi with a union of shifts
# covering that window recovers the branch leaders to ~1e-14 — but sequential per-pencil Krylov
# is 2x SLOWER than geev on CPU. The GPU changes the economics: EVERY step of shift-invert
# subspace iteration exists as a batched cuBLAS primitive, so P pencils are processed in ONE
# kernel launch sequence per shift:
#   G  = A - sigma*B          (broadcast)
#   G  = LU(G)                (getrf_strided_batched!)          NOTE: getrs_batched segfaults in
#   Gi = G^-1                 (getri_strided_batched!)          this cublas build, hence explicit
#   S  = Gi * B               (gemm_strided_batched!)           inverse; shifts keep G well-cond.
#   q x subspace iteration:  X <- orth(S*X)   (gemm + thin host QR)
#   Rayleigh-Ritz: T = X^H S X (m x m, host eig), mu -> lambda = sigma + 1/mu
# Union over shifts, dedup, then compare TJLF-consumed outputs vs full CPU geev.
#
# Usage: julia --project=. dev/ad/benchmark_batched_si_gpu.jl
# Env: PENCILS, M (16 subspace), Q (10 iters), NMODES (4), EPS1 (1e-12),
#      SHIFTS ("0.02,0.02+0.15im,0.02-0.15im,0.05+0.5im,0.05-0.5im,0.05+1.0im,0.05-1.0im")
using Pkg
Pkg.activate(normpath(@__DIR__, "..", ".."))
import CUDA
using CUDA: CuArray
using LinearAlgebra, Serialization, Printf, Statistics
import LinearAlgebra.LAPACK: gesv!, geev!

const PDIR   = get(ENV, "PENCILS", normpath(@__DIR__, "pencils_nb16"))
const M      = parse(Int, get(ENV, "M", "16"))
const Q      = parse(Int, get(ENV, "Q", "10"))
const NMODES = parse(Int, get(ENV, "NMODES", "4"))
const EPS1   = parse(Float64, get(ENV, "EPS1", "1e-12"))
const METHOD = Symbol(get(ENV, "METHOD", "trsm"))   # :trsm (optimized) or :inv (getri baseline)
const SHIFTS = [parse(ComplexF64, s) for s in
                split(get(ENV, "SHIFTS",
                          "0.02,0.02+0.05im,0.02-0.05im,0.02+0.12im,0.02-0.12im,0.02+0.25im,0.02-0.25im,0.05+0.6im,0.05-0.6im,0.05+1.1im,0.05-1.1im,0.05+1.5im,0.05-1.5im"), ",")]

# GPU gather of the M-column block by per-pencil LU pivots: dst[i,c,p] = src[perm[i,p],c,p].
# 3D launch (row / col / pencil) avoids div/rem index decoding.
function _gather_rows_kernel!(dst, src, perm, n)
    i = (CUDA.blockIdx().x - 1) * CUDA.blockDim().x + CUDA.threadIdx().x   # row
    c = CUDA.blockIdx().y                                                 # column
    p = CUDA.blockIdx().z                                                 # pencil
    if i <= n
        @inbounds dst[i, c, p] = src[perm[i, p], c, p]
    end
    return
end
function _gather_rows!(dst, src, perm)
    n, mm, P = size(dst)
    CUDA.@cuda threads=256 blocks=(cld(n, 256), mm, P) _gather_rows_kernel!(dst, src, perm, n)
    return dst
end

# Convert LAPACK-style ipiv (sequential row swaps from getrf) to a permutation `perm` with
# (P·b)[i] = b[perm[i]]. Host, O(n) per pencil.
function _ipiv_to_perm(ipiv::AbstractMatrix{<:Integer})
    n, P = size(ipiv)
    perm = Matrix{Int32}(undef, n, P)
    idx = Vector{Int}(undef, n)
    for p in 1:P
        @inbounds for i in 1:n; idx[i] = i; end
        @inbounds for i in 1:n
            j = ipiv[i, p]
            if j != i; idx[i], idx[j] = idx[j], idx[i]; end
        end
        @inbounds for i in 1:n; perm[i, p] = idx[i]; end
    end
    return perm
end

BLAS.set_num_threads(8)

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

const ORTH = Symbol(get(ENV, "ORTH", "cholqr"))   # :cholqr (GPU-resident) or :hostqr

# Host QR (thin): transfers the full n*M*P block both ways every call (PCIe-heavy at large n).
function _orth_hostqr!(X, scratch)
    Xh = copyto!(scratch, X)
    n, mm, P = size(Xh)
    Threads.@threads for p in 1:P
        @views Xh[:, :, p] .= Matrix(qr!(Xh[:, :, p]).Q)
    end
    copyto!(X, Xh)
    return X
end

# Cholesky-QR: X = QR with R from chol(XᴴX). Only the M×M Gram/Rinv cross PCIe (tiny), the
# n×M block stays on the GPU. gram_d/rinv_d/tmp_d are reusable device buffers.
function _orth_cholqr!(X, gram_d, rinv_d, tmp_d)
    n, mm, P = size(X)
    CUDA.CUBLAS.gemm_strided_batched!('C','N', ComplexF64(1), X, X, ComplexF64(0), gram_d)
    Gh = Array(gram_d)                      # M×M×P, tiny
    Rinv = similar(Gh)
    Threads.@threads for p in 1:P
        @views begin
            H = Hermitian((Gh[:, :, p] .+ Gh[:, :, p]') ./ 2)
            R = cholesky(H).U
            Rinv[:, :, p] .= inv(R)
        end
    end
    copyto!(rinv_d, Rinv)
    CUDA.CUBLAS.gemm_strided_batched!('N','N', ComplexF64(1), X, rinv_d, ComplexF64(0), tmp_d)
    copyto!(X, tmp_d)
    return X
end

_orth!(X, bufs) = ORTH === :cholqr ?
    _orth_cholqr!(X, bufs.gram, bufs.rinv, bufs.tmp) : _orth_hostqr!(X, bufs.scratch)

# :inv method — baseline: explicit inverse + dense S = G^-1 B (2.5 n^3 per shift).
function _si_shift_inv!(cands, A3g, B3g, sigma, bufs)
    G = A3g .- sigma .* B3g
    piv, _, _ = CUDA.CUBLAS.getrf_strided_batched!(G, true)
    Gi = similar(G); CUDA.CUBLAS.getri_strided_batched!(G, Gi, piv)
    S = similar(G); CUDA.CUBLAS.gemm_strided_batched!('N','N', ComplexF64(1), Gi, B3g, ComplexF64(0), S)
    CUDA.unsafe_free!(G); CUDA.unsafe_free!(Gi)
    apply_S! = (dst, src) -> CUDA.CUBLAS.gemm_strided_batched!('N','N', ComplexF64(1), S, src, ComplexF64(0), dst)
    _ritz!(cands, apply_S!, sigma, bufs)
    CUDA.unsafe_free!(S)
    return
end

# :trsm method — optimized: factor G once, apply S x = U^-1 L^-1 P (B x) to the M-block only.
# No explicit inverse, no dense S: ~1 n^3 per shift. Pivots applied to the n x M block via gather.
function _si_shift_trsm!(cands, A3g, B3g, sigma, bufs)
    n, _, P = size(A3g)
    G = A3g .- sigma .* B3g
    piv, _, _ = CUDA.CUBLAS.getrf_strided_batched!(G, true)   # G holds LU in place
    perm = CuArray(_ipiv_to_perm(Array(piv)))
    Lv = [view(G, :, :, p) for p in 1:P]
    permbuf = bufs.permbuf
    # S x = block-solve: W = B x (gemm) -> PW = P W (gather) -> L\PW -> U\ : result in dst
    apply_S! = (dst, src) -> begin
        CUDA.CUBLAS.gemm_strided_batched!('N','N', ComplexF64(1), B3g, src, ComplexF64(0), permbuf)
        _gather_rows!(dst, permbuf, perm)
        Dv = [view(dst, :, :, p) for p in 1:P]
        CUDA.CUBLAS.trsm_batched!('L','L','N','U', ComplexF64(1), Lv, Dv)
        CUDA.CUBLAS.trsm_batched!('L','U','N','N', ComplexF64(1), Lv, Dv)
        return dst
    end
    _ritz!(cands, apply_S!, sigma, bufs)
    CUDA.unsafe_free!(G)
    return
end

# Shared subspace iteration + Rayleigh-Ritz given an operator apply_S!(dst, src) = S*src.
function _ritz!(cands, apply_S!, sigma, bufs)
    X, Y = bufs.X, bufs.Y
    n, mm, P = size(X)
    for _ in 1:Q
        apply_S!(Y, X)
        _orth!(Y, bufs)
        X, Y = Y, X            # swap: X holds latest orthonormal block
    end
    apply_S!(Y, X)             # Y = S X
    T = bufs.gram
    CUDA.CUBLAS.gemm_strided_batched!('C','N', ComplexF64(1), X, Y, ComplexF64(0), T)
    Th = Array(T)
    Threads.@threads for p in 1:P
        mus = eigvals(Th[:, :, p])
        append!(cands[p], [sigma + 1 / mu for mu in mus if abs(mu) > 1e-10])
    end
    return
end

# One batched shift; buffers reused across shifts. X0h is the shared (n,M) host start block.
function batched_si_shift!(cands, A3g, B3g, sigma, X0h, bufs)
    n, _, P = size(A3g)
    copyto!(bufs.X, repeat(reshape(X0h, n, M, 1), 1, 1, P))
    METHOD === :trsm ? _si_shift_trsm!(cands, A3g, B3g, sigma, bufs) :
                       _si_shift_inv!(cands, A3g, B3g, sigma, bufs)
    return
end

function main()
    files = sort(filter(f -> endswith(f, ".jls"), readdir(PDIR; join=true)))
    isempty(files) && error("no pencils in $PDIR")
    allpens = [Serialization.deserialize(f) for f in files]
    # A grid scan emits pencils at several sizes (coarse width-search vs full-nbasis). Keep only
    # the modal size so the batched stack is uniform.
    szs = [size(p.A, 1) for p in allpens]
    usz = unique(szs); modal = usz[argmax([count(==(s), szs) for s in usz])]
    length(usz) > 1 && println("pencil sizes present: ",
        join(["n=$s×$(count(==(s),szs))" for s in sort(usz)], " "), " -> using n=", modal)
    pens = [p for p in allpens if size(p.A, 1) == modal]
    P = length(pens); n = modal
    println("pencils: ", P, "   n = ", n, "   M=", M, " Q=", Q, "   shifts=", length(SHIFTS),
            "   method=", METHOD)
    CUDA.functional() || error("no functional GPU")

    # reference spectra (CPU, threaded across pencils)
    t_ref = zeros(P); refs = Vector{Vector{ComplexF64}}(undef, P)
    Threads.@threads for p in 1:P
        refs[p], t_ref[p] = run_ref(pens[p].A, pens[p].B)
    end

    # stack on device
    A3h = zeros(ComplexF64, n, n, P); B3h = zeros(ComplexF64, n, n, P)
    for p in 1:P
        A3h[:, :, p] = pens[p].A; B3h[:, :, p] = pens[p].B
    end
    t_up = @elapsed begin
        A3g = CuArray(A3h); B3g = CuArray(B3h); CUDA.synchronize()
    end

    X0h = Matrix(qr(randn(ComplexF64, n, M)).Q)
    bufs = (X = CUDA.zeros(ComplexF64, n, M, P), Y = CUDA.zeros(ComplexF64, n, M, P),
            permbuf = CUDA.zeros(ComplexF64, n, M, P), tmp = CUDA.zeros(ComplexF64, n, M, P),
            gram = CUDA.zeros(ComplexF64, M, M, P), rinv = CUDA.zeros(ComplexF64, M, M, P),
            scratch = Array{ComplexF64}(undef, n, M, P))
    println("orth=", ORTH)
    cands = [ComplexF64[] for _ in 1:P]
    # warm one shift pass to exclude compile time
    batched_si_shift!([ComplexF64[] for _ in 1:P], A3g, B3g, SHIFTS[1], X0h, bufs)
    t_gpu = @elapsed begin
        for sigma in SHIFTS
            batched_si_shift!(cands, A3g, B3g, sigma, X0h, bufs)
        end
        CUDA.synchronize()
    end

    # dedup per pencil
    for p in 1:P
        sort!(cands[p]; by=x->(real(x), imag(x)))
        ded = ComplexF64[]
        for l in cands[p]
            (isempty(ded) || abs(l - ded[end]) > 1e-8 * max(1.0, abs(l))) && push!(ded, l)
        end
        cands[p] = ded
    end

    # accuracy
    err_top = Float64[]; gerr_ion = Float64[]; gerr_ele = Float64[]
    n_ion = 0; n_ele = 0; n_ionmiss = 0; n_elemiss = 0
    for p in 1:P
        tr_ref = topk(refs[p], NMODES); tr_si = topk(cands[p], NMODES)
        m = min(length(tr_ref), length(tr_si))
        m > 0 && push!(err_top, maximum(abs.(tr_ref[1:m] .- tr_si[1:m])))
        bl = branch_leaders(refs[p]); bls = branch_leaders(cands[p])
        if bl.ion !== nothing
            n_ion += 1
            bls.ion === nothing ? (n_ionmiss += 1) : push!(gerr_ion, abs(real(bl.ion) - real(bls.ion)))
        end
        if bl.ele !== nothing
            n_ele += 1
            bls.ele === nothing ? (n_elemiss += 1) : push!(gerr_ele, abs(real(bl.ele) - real(bls.ele)))
        end
    end

    q(v, pr) = isempty(v) ? NaN : quantile(v, pr)
    println("\n================ timings ================")
    @printf("  CPU geev ref        : median %6.1f ms/pencil (threaded ref pass total %.1f s)\n",
            1e3*q(t_ref,0.5), sum(t_ref))
    @printf("  H2D upload          : %6.1f ms total\n", 1e3*t_up)
    @printf("  GPU batched SI      : %6.1f ms TOTAL for %d pencils x %d shifts = %5.2f ms/pencil\n",
            1e3*t_gpu, P, length(SHIFTS), 1e3*t_gpu/P)
    @printf("  speedup vs 1-thread CPU geev: %.1fx per pencil\n", q(t_ref,0.5)/(t_gpu/P))

    println("\n================ accuracy vs full geev ================")
    @printf("  top-%d |dLambda|: median %.3g   p90 %.3g   max %.3g\n",
            NMODES, q(err_top,0.5), q(err_top,0.9), isempty(err_top) ? NaN : maximum(err_top))
    @printf("  ion leader : missed %d/%d   |dgamma| median %.3g max %.3g\n",
            n_ionmiss, n_ion, q(gerr_ion,0.5), isempty(gerr_ion) ? NaN : maximum(gerr_ion))
    @printf("  ele leader : missed %d/%d   |dgamma| median %.3g max %.3g\n",
            n_elemiss, n_ele, q(gerr_ele,0.5), isempty(gerr_ele) ? NaN : maximum(gerr_ele))
end

main()
