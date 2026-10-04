#!/usr/bin/env julia
# Smoke test: the GPU contour-moment path (TJLFCUDAExt._contour_moments — batched
# getrf/gather/trsm/gemm, no Xgeev) must agree with the CPU reference path on planted spectra.
# Prints CONTOUR_GPU_OK or CONTOUR_GPU_MISMATCH.
using Pkg
Pkg.activate(normpath(@__DIR__, "..", ".."))
import CUDA
using TJLF
using LinearAlgebra, Random, Printf

CUDA.functional() || error("no functional GPU")
Random.seed!(11)

n, P = 96, 6
cfg = TJLF.ContourConfig(re_lo = -0.05, re_hi = 0.4, im_max = 1.2,
                         n_long = 12, n_short = 6, L = 12, K = 2)

As = Vector{Matrix{ComplexF64}}(undef, P)
Bs = Vector{Matrix{ComplexF64}}(undef, P)
inside = Vector{Vector{ComplexF64}}(undef, P)
for p in 1:P
    λin = [complex(-0.03 + 0.4 * rand(), -1.0 + 2.0 * rand()) for _ in 1:4]
    λin = [λ for λ in λin if TJLF._contour_inside(λ, cfg) && abs(imag(λ)) < cfg.im_max - 0.05]
    bulk = ComplexF64[(20 + 30 * rand()) * cis(π * (0.55 + 0.9 * rand()))
                      for _ in 1:(n - length(λin))]
    Q = Matrix(qr(randn(ComplexF64, n, n)).Q)
    M = Q * Diagonal(vcat(λin, bulk)) * Q'
    Bs[p] = Matrix{ComplexF64}(I, n, n) + 0.1 * randn(ComplexF64, n, n)
    As[p] = Bs[p] * M
    inside[p] = λin
end

vg, flg, rkg = TJLF.contour_eigvals_batch(As, Bs; cfg, use_gpu = true,  dense_fallback = :none)
vc, flc, rkc = TJLF.contour_eigvals_batch(As, Bs; cfg, use_gpu = false, dense_fallback = :none)

worst_gpu = worst_cpu = 0.0
spurious = 0
for p in 1:P
    for λ in inside[p]
        global worst_gpu = max(worst_gpu, minimum(abs(λ - g) for g in vg[p]; init = Inf))
        global worst_cpu = max(worst_cpu, minimum(abs(λ - g) for g in vc[p]; init = Inf))
    end
    global spurious += count(g -> minimum(abs(λ - g) for λ in inside[p]; init = Inf) > 1e-6, vg[p])
end
@printf("planted-mode recovery: worst |dLambda| gpu=%.3g cpu=%.3g   spurious(gpu)=%d\n",
        worst_gpu, worst_cpu, spurious)
@printf("flagged: gpu=%d cpu=%d   ranks gpu=%s cpu=%s\n", count(flg), count(flc), rkg, rkc)
ok = worst_gpu < 1e-7 && spurious == 0 && count(flg) == 0
println(ok ? "CONTOUR_GPU_OK" : "CONTOUR_GPU_MISMATCH")
