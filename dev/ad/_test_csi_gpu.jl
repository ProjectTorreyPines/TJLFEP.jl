#!/usr/bin/env julia
# GPU smoke test for the coverage-certified adaptive SI solver:
#  (a) the GPU handle's certified Ritz values match the CPU handle's on planted spectra
#      (both per-shift batched sweeps with per-pencil shift vectors);
#  (b) the full certified_si_eigvals_batch(use_gpu=true) recovers every planted unstable mode.
using Pkg
Pkg.activate(normpath(@__DIR__, "..", ".."))
import CUDA
using TJLF
using LinearAlgebra, Random, Printf

CUDA.functional() || error("no functional GPU")

Random.seed!(33)
n, P = 96, 6
cfg = TJLF.CertifiedSIConfig(re_hi = 0.6, im_max = 1.5, row_dy = 0.5, M = 16, Q = 16,
                             trust = 8, max_rounds = 60)
As = Vector{Matrix{ComplexF64}}(undef, P); Bs = similar(As)
unstable = Vector{Vector{ComplexF64}}(undef, P)
for p in 1:P
    λu = [complex(0.02 + 0.55 * rand(), -1.4 + 2.8 * rand()) for _ in 1:(2 + p % 3)]
    crowd = [complex(-0.06 * rand(), -1.4 + 2.8 * rand()) for _ in 1:30]
    bulk = ComplexF64[(15 + 25 * rand()) * cis(π * (0.55 + 0.9 * rand()))
                      for _ in 1:(n - length(λu) - 30)]
    Q0 = Matrix(qr(randn(ComplexF64, n, n)).Q)
    M0 = Q0 * Diagonal(vcat(λu, crowd, bulk)) * Q0'
    Bs[p] = Matrix{ComplexF64}(I, n, n) + 0.1 * randn(ComplexF64, n, n)
    As[p] = Bs[p] * M0
    unstable[p] = λu
end

# (a) one batched sweep, per-pencil shifts: GPU handle vs CPU handle
σs = [complex(0.02 + 0.1 * p, 0.3 * p - 1.0) for p in 1:P]
hg = TJLF._CUDA_CSI_OPEN[](As, Bs, cfg.M, cfg.Q)
λg, rg = hg.solve(σs)
hg.close()
hc = TJLF._csi_open_cpu(As, Bs, cfg.M, cfg.Q)
λc, rc = hc.solve(σs)
nbad = 0
for p in 1:P
    certg = sort([λg[i, p] for i in 1:cfg.M if rg[i, p] <= cfg.resid_tol]; by = x -> (real(x), imag(x)))
    certc = sort([λc[i, p] for i in 1:cfg.M if rc[i, p] <= cfg.resid_tol]; by = x -> (real(x), imag(x)))
    dev = (isempty(certg) || length(certg) != length(certc)) ?
          (isempty(certg) && isempty(certc) ? 0.0 : Inf) :
          maximum(abs.(certg .- certc))
    @printf("p%d: ncert gpu=%d cpu=%d  max|dl|=%.3g\n", p, length(certg), length(certc), dev)
    dev < 1e-8 || (global nbad += 1)
end

# (b) full solver on GPU: every planted unstable mode recovered (or pencil flagged)
vals, fl, ns = TJLF.certified_si_eigvals_batch(As, Bs; cfg, use_gpu = true,
                                               dense_fallback = :cpu)
for p in 1:P
    miss = [λ for λ in unstable[p] if minimum(abs.(vals[p] .- λ); init = Inf) > 1e-6]
    @printf("p%d: fl=%d shifts=%d found=%d missed=%d\n", p, fl[p], ns[p], length(vals[p]),
            length(miss))
    isempty(miss) || (global nbad += 1)
end

println(nbad == 0 ? "CSI_GPU_OK" : "CSI_GPU_BAD n=$nbad")
