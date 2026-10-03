#!/usr/bin/env julia
# Adaptive-shift salvage of batched shift-invert.
#
# The fixed 13-shift set misses ion-branch leaders at high-drive radii (IR101, factor=10: ion
# missed 333/503) because the leaders sit where no shift was placed. Fix: DERIVE the shift set
# from the actual spectrum — geev a cheap CALIBRATION subset of the batch, collect the observed
# branch leaders (both ion imag>0 and electron imag<=0), cluster them, and place a shift at each
# cluster. Coverage then follows the drive level automatically. The bulk of the batch is solved
# with batched SI using those data-driven shifts. Optionally shard the batch over NGPU A100s.
#
# Usage: NGPU=1 CALIB_FRAC=0.1 julia --project=. dev/ad/benchmark_adaptive_si_gpu.jl
# Env: PENCILS, M(16), Q(12), NMODES(4), EPS1(1e-12), CALIB_FRAC(0.1), CLUSTER_TOL(0.12),
#      MAXSHIFTS(28), PER_BRANCH(2), OFFSET(0.02), NGPU(1)
using Pkg
Pkg.activate(normpath(@__DIR__, "..", ".."))
import CUDA
using CUDA: CuArray
using LinearAlgebra, Serialization, Printf, Statistics
import LinearAlgebra.LAPACK: gesv!, geev!

const PDIR       = get(ENV, "PENCILS", normpath(@__DIR__, "pencils_nb16"))
const M          = parse(Int, get(ENV, "M", "16"))
const Q          = parse(Int, get(ENV, "Q", "12"))
const NMODES     = parse(Int, get(ENV, "NMODES", "4"))
const EPS1       = parse(Float64, get(ENV, "EPS1", "1e-12"))
const CALIB_FRAC = parse(Float64, get(ENV, "CALIB_FRAC", "0.1"))
const CLUSTER_TOL= parse(Float64, get(ENV, "CLUSTER_TOL", "0.12"))
const MAXSHIFTS  = parse(Int, get(ENV, "MAXSHIFTS", "28"))
const PER_BRANCH = parse(Int, get(ENV, "PER_BRANCH", "2"))
const OFFSET     = parse(Float64, get(ENV, "OFFSET", "0.02"))
const NGPU       = parse(Int, get(ENV, "NGPU", "1"))
# UNION_FIXED=1: keep the full fixed set (guarantees the dense near-axis electron band) and ADD
# adaptive ion-branch shifts on top, rather than redistributing a fixed budget. Prevents the
# electron branch from degrading when the ion cloud is spread out (high drive).
const UNION_FIXED = get(ENV, "UNION_FIXED", "0") == "1"
const FIXED_SHIFTS = ComplexF64[0.02, 0.02+0.05im,0.02-0.05im, 0.02+0.12im,0.02-0.12im,
                                0.02+0.25im,0.02-0.25im, 0.05+0.6im,0.05-0.6im,
                                0.05+1.1im,0.05-1.1im, 0.05+1.5im,0.05-1.5im]

# ---------- batched SI primitives (identical math to TJLFCUDAExt._si_*) --------------------------
function _gather_rows_kernel!(dst, src, perm, n)
    i = (CUDA.blockIdx().x - 1) * CUDA.blockDim().x + CUDA.threadIdx().x
    c = CUDA.blockIdx().y; p = CUDA.blockIdx().z
    if i <= n; @inbounds dst[i, c, p] = src[perm[i, p], c, p]; end
    return
end
function _gather_rows!(dst, src, perm)
    n, mm, P = size(dst)
    CUDA.@cuda threads=256 blocks=(cld(n,256), mm, P) _gather_rows_kernel!(dst, src, perm, n)
    return dst
end
function _ipiv_to_perm(ipiv::AbstractMatrix{<:Integer})
    n, P = size(ipiv); perm = Matrix{Int32}(undef, n, P); idx = Vector{Int}(undef, n)
    for p in 1:P
        @inbounds for i in 1:n; idx[i] = i; end
        @inbounds for i in 1:n
            j = ipiv[i, p]; j != i && ((idx[i], idx[j]) = (idx[j], idx[i]))
        end
        @inbounds for i in 1:n; perm[i, p] = idx[i]; end
    end
    return perm
end
function _cholqr!(X, gram, rinv, tmp)
    _, mm, P = size(X)
    CUDA.CUBLAS.gemm_strided_batched!('C','N', ComplexF64(1), X, X, ComplexF64(0), gram)
    Gh = Array(gram); Rinv = similar(Gh)
    Threads.@threads for p in 1:P
        @views Rinv[:, :, p] .= inv(cholesky(Hermitian((Gh[:,:,p] .+ Gh[:,:,p]') ./ 2)).U)
    end
    copyto!(rinv, Rinv)
    CUDA.CUBLAS.gemm_strided_batched!('N','N', ComplexF64(1), X, rinv, ComplexF64(0), tmp)
    copyto!(X, tmp); return X
end
# One shift σ over the whole (already-on-device) batch; appends Ritz values to cands.
function _si_shift!(cands, A3, B3, σ, X0h, bufs, mu_tol)
    n, _, P = size(A3); X, Y, permbuf, tmp, gram, rinv = bufs
    copyto!(X, repeat(reshape(X0h, n, M, 1), 1, 1, P))
    G = A3 .- σ .* B3
    piv, _, _ = CUDA.CUBLAS.getrf_strided_batched!(G, true)
    perm = CuArray(_ipiv_to_perm(Array(piv)))
    Lv = [view(G, :, :, p) for p in 1:P]
    apply_S! = (dst, src) -> begin
        CUDA.CUBLAS.gemm_strided_batched!('N','N', ComplexF64(1), B3, src, ComplexF64(0), permbuf)
        _gather_rows!(dst, permbuf, perm)
        Dv = [view(dst, :, :, p) for p in 1:P]
        CUDA.CUBLAS.trsm_batched!('L','L','N','U', ComplexF64(1), Lv, Dv)
        CUDA.CUBLAS.trsm_batched!('L','U','N','N', ComplexF64(1), Lv, Dv)
        return dst
    end
    for _ in 1:Q
        apply_S!(Y, X); _cholqr!(Y, gram, rinv, tmp); X, Y = Y, X
    end
    apply_S!(Y, X)
    T = gram; CUDA.CUBLAS.gemm_strided_batched!('C','N', ComplexF64(1), X, Y, ComplexF64(0), T)
    Th = Array(T); CUDA.unsafe_free!(G)
    Threads.@threads for p in 1:P
        μs = eigvals(@view Th[:, :, p])
        append!(cands[p], (σ + 1/μ for μ in μs if abs(μ) > mu_tol))
    end
    return
end
# Solve a shard of pencils on the CURRENT device with a given shift list. Returns candidate lists.
function si_solve_shard(As, Bs, shifts)
    P = length(As); n = size(As[1], 1)
    A3h = Array{ComplexF64}(undef, n, n, P); B3h = similar(A3h)
    @inbounds for p in 1:P; A3h[:,:,p] .= As[p]; B3h[:,:,p] .= Bs[p]; end
    A3 = CuArray(A3h); B3 = CuArray(B3h)
    bufs = (CUDA.zeros(ComplexF64,n,M,P), CUDA.zeros(ComplexF64,n,M,P), CUDA.zeros(ComplexF64,n,M,P),
            CUDA.zeros(ComplexF64,n,M,P), CUDA.zeros(ComplexF64,M,M,P), CUDA.zeros(ComplexF64,M,M,P))
    X0h = Matrix(qr(randn(ComplexF64, n, M)).Q)
    cands = [ComplexF64[] for _ in 1:P]
    for σ in shifts; _si_shift!(cands, A3, B3, σ, X0h, bufs, 1e-10); end
    CUDA.synchronize()
    return [_dedup(cands[p], 1e-8) for p in 1:P]
end
# Multi-GPU: split pencils into NGPU contiguous shards, one task per device.
function si_solve_multigpu(As, Bs, shifts, ngpu)
    P = length(As); ngpu = min(ngpu, length(CUDA.devices()), P)
    ngpu <= 1 && (CUDA.device!(0); return si_solve_shard(As, Bs, shifts))
    bnd = round.(Int, range(0, P; length = ngpu + 1))
    parts = Vector{Vector{Vector{ComplexF64}}}(undef, ngpu)
    @sync for g in 1:ngpu
        rng = (bnd[g]+1):bnd[g+1]
        Threads.@spawn begin
            CUDA.device!(g - 1)
            parts[g] = si_solve_shard(As[rng], Bs[rng], shifts)
        end
    end
    return reduce(vcat, parts)
end

# ---------- accuracy helpers ---------------------------------------------------------------------
_dedup(lams, tol) = (isempty(lams) ? ComplexF64[] : begin
    s = sort(lams; by = x -> (real(x), imag(x))); out = ComplexF64[s[1]]
    for l in @view s[2:end]; abs(l - out[end]) > tol*max(1.0,abs(l)) && push!(out, l); end; out end)
topk(vals, k) = sort(vals; by=real, rev=true)[1:min(k, length(vals))]
function branch_leaders(vals)
    ion = filter(v -> real(v) > EPS1 && imag(v) > 0, vals)
    ele = filter(v -> real(v) > EPS1 && imag(v) <= 0, vals)
    (ion = isempty(ion) ? nothing : ion[argmax(real.(ion))],
     ele = isempty(ele) ? nothing : ele[argmax(real.(ele))])
end
function run_ref(A, B)
    A2, B2 = copy(A), copy(B); (A3,_,_) = gesv!(B2, A2); return geev!('N','N', A3)[1]
end

# Build adaptive shifts from calibration spectra: take the top PER_BRANCH unstable modes of EACH
# branch per calib pencil, cluster the union in the complex plane, place a shift at each cluster
# centroid (nudged +OFFSET in real to stay off the eigenvalue). Coverage tracks the real spectrum.
function build_adaptive_shifts(refs_calib)
    pts = ComplexF64[]
    for vals in refs_calib
        ion = sort(filter(v->real(v)>EPS1 && imag(v)>0,  vals); by=v->-real(v))
        ele = sort(filter(v->real(v)>EPS1 && imag(v)<=0, vals); by=v->-real(v))
        append!(pts, ion[1:min(PER_BRANCH,end)]); append!(pts, ele[1:min(PER_BRANCH,end)])
    end
    isempty(pts) && return copy(FIXED_SHIFTS)
    sort!(pts; by=v->(real(v), imag(v)))
    clusters = Vector{Vector{ComplexF64}}()
    for z in pts
        placed = false
        for c in clusters
            if abs(z - sum(c)/length(c)) <= CLUSTER_TOL; push!(c, z); placed = true; break; end
        end
        placed || push!(clusters, [z])
    end
    sort!(clusters; by = c -> -length(c))
    cents = ComplexF64[complex(real(sum(c)/length(c)) + OFFSET, imag(sum(c)/length(c))) for c in clusters]
    if UNION_FIXED
        # Start from the guaranteed fixed set, then append adaptive centroids that aren't already
        # covered (farther than CLUSTER_TOL from an existing shift), up to MAXSHIFTS.
        out = copy(FIXED_SHIFTS)
        for z in cents
            length(out) >= MAXSHIFTS && break
            minimum(abs(z - s) for s in out) > CLUSTER_TOL && push!(out, z)
        end
        return out
    end
    return cents[1:min(MAXSHIFTS, end)]
end

function score(cands, refs, label, shifts, t_gpu, P)
    gerr_ion = Float64[]; gerr_ele = Float64[]; err_top = Float64[]
    n_ion=0; n_ele=0; n_ionmiss=0; n_elemiss=0
    for p in 1:P
        tr_ref = topk(refs[p], NMODES); tr_si = topk(cands[p], NMODES)
        m = min(length(tr_ref), length(tr_si)); m>0 && push!(err_top, maximum(abs.(tr_ref[1:m] .- tr_si[1:m])))
        bl = branch_leaders(refs[p]); bls = branch_leaders(cands[p])
        if bl.ion !== nothing; n_ion += 1
            bls.ion === nothing ? (n_ionmiss += 1) : push!(gerr_ion, abs(real(bl.ion)-real(bls.ion))); end
        if bl.ele !== nothing; n_ele += 1
            bls.ele === nothing ? (n_elemiss += 1) : push!(gerr_ele, abs(real(bl.ele)-real(bls.ele))); end
    end
    q(v,pr) = isempty(v) ? NaN : quantile(v,pr)
    @printf("\n---- %s : %d shifts, %.1f ms/pencil ----\n", label, length(shifts), 1e3*t_gpu/P)
    @printf("  top-%d |dLambda|: median %.3g  p90 %.3g  max %.3g\n",
            NMODES, q(err_top,0.5), q(err_top,0.9), isempty(err_top) ? NaN : maximum(err_top))
    @printf("  ion leader : MISSED %d/%d   |dgamma| median %.3g max %.3g\n",
            n_ionmiss, n_ion, q(gerr_ion,0.5), isempty(gerr_ion) ? NaN : maximum(gerr_ion))
    @printf("  ele leader : MISSED %d/%d   |dgamma| median %.3g max %.3g\n",
            n_elemiss, n_ele, q(gerr_ele,0.5), isempty(gerr_ele) ? NaN : maximum(gerr_ele))
    return n_ionmiss + n_elemiss
end

function main()
    files = sort(filter(f -> endswith(f, ".jls"), readdir(PDIR; join=true)))
    isempty(files) && error("no pencils in $PDIR")
    allpens = [Serialization.deserialize(f) for f in files]
    szs = [size(p.A,1) for p in allpens]; usz = unique(szs)
    modal = usz[argmax([count(==(s), szs) for s in usz])]
    pens = [p for p in allpens if size(p.A,1) == modal]
    As = [ComplexF64.(p.A) for p in pens]; Bs = [ComplexF64.(p.B) for p in pens]
    P = length(pens); n = modal
    CUDA.functional() || error("no functional GPU")
    @printf("pencils=%d  n=%d  M=%d Q=%d  NGPU=%d  devices=%d\n", P, n, M, Q, NGPU, length(CUDA.devices()))

    # geev reference for ALL pencils (threaded)
    refs = Vector{Vector{ComplexF64}}(undef, P)
    t_ref = @elapsed (Threads.@threads for p in 1:P; refs[p] = run_ref(As[p], Bs[p]); end)
    @printf("geev ref (threaded, all %d): %.1f s total\n", P, t_ref)

    # calibration subset (evenly strided) -> adaptive shifts
    stride = max(1, round(Int, 1/CALIB_FRAC)); calib_idx = collect(1:stride:P)
    adaptive = build_adaptive_shifts(refs[calib_idx])
    @printf("calibration: %d/%d pencils (frac~%.2f) -> %d adaptive shifts\n",
            length(calib_idx), P, length(calib_idx)/P, length(adaptive))

    # warm
    si_solve_multigpu(As[1:min(P,8)], Bs[1:min(P,8)], FIXED_SHIFTS[1:1], NGPU)

    t_fix = @elapsed cf = si_solve_multigpu(As, Bs, FIXED_SHIFTS, NGPU)
    miss_fix = score(cf, refs, "FIXED (13-shift)", FIXED_SHIFTS, t_fix, P)
    t_ada = @elapsed ca = si_solve_multigpu(As, Bs, adaptive, NGPU)
    miss_ada = score(ca, refs, "ADAPTIVE (geev-calibrated)", adaptive, t_ada, P)

    @printf("\n================ verdict ================\n")
    @printf("  total leader misses:  fixed=%d   adaptive=%d\n", miss_fix, miss_ada)
    @printf("  batched-SI wall:      fixed=%.1fs  adaptive=%.1fs  (NGPU=%d)\n", t_fix, t_ada, NGPU)
    @printf("  effective speedup vs full geev (adaptive, incl calib): %.1fx\n",
            t_ref / (t_ada + t_ref*length(calib_idx)/P))

    # Optional machine-readable row for the timing-vs-nbasis scan.
    csv = get(ENV, "CSV_OUT", "")
    if !isempty(csv)
        nb = get(ENV, "NB", "?")
        open(csv, "a") do io
            @printf(io, "%s,%d,%d,%d,%.3f,%.3f,%.3f,%d\n", nb, n, P, NGPU,
                    1e3*t_fix/P, 1e3*t_ada/P, t_ref/(t_ada + t_ref*length(calib_idx)/P), miss_ada)
        end
    end
end
main()
