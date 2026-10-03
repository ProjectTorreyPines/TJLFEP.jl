#!/usr/bin/env julia
# CPU validation of the SS-RR contour solver on REAL harvested pencils (no GPU needed):
# recovery of all unstable in-window modes + branch leaders vs dense geev, spurious count,
# rank/flag behavior. Env: PENCILS (default pencils_nb16), NP (default 8, pencil count).
using Pkg
Pkg.activate(normpath(@__DIR__, "..", ".."))
using TJLF
using LinearAlgebra, Serialization, Printf
import LinearAlgebra.LAPACK: gesv!, geev!

const PDIR = get(ENV, "PENCILS", normpath(@__DIR__, "pencils_nb16"))
const NP   = parse(Int, get(ENV, "NP", "8"))
const EPS1 = 1e-12

files = sort(filter(f -> endswith(f, ".jls"), readdir(PDIR; join=true)))
pens = [Serialization.deserialize(f) for f in files]
szs = [size(p.A, 1) for p in pens]; usz = unique(szs)
modal = usz[argmax([count(==(s), szs) for s in usz])]
pens = [p for p in pens if size(p.A, 1) == modal][1:min(NP, end)]
As = [ComplexF64.(p.A) for p in pens]; Bs = [ComplexF64.(p.B) for p in pens]
P = length(As); n = modal
println("pencils=", P, " n=", n)

refs = Vector{Vector{ComplexF64}}(undef, P)
t_ref = @elapsed Threads.@threads for p in 1:P
    A2, B2 = copy(As[p]), copy(Bs[p]); (A3, _, _) = gesv!(B2, A2)
    refs[p] = geev!('N', 'N', A3)[1]
end
@printf("geev refs: %.1f s\n", t_ref)

cfg = TJLF.ContourConfig()
nin = [count(v -> TJLF._contour_inside(v, cfg), r) for r in refs]
println("in-window counts: ", nin, "  capacity=", 2 * cfg.K * cfg.L)

t = @elapsed vals, fl, rk = TJLF.contour_eigvals_batch(As, Bs; cfg, use_gpu=false,
                                                       dense_fallback=:none)
@printf("cpu contour: %.1f s total (%.1f s/pencil)\nranks=%s\nflagged=%s\n", t, t / P, rk, fl)

branch(vals) = begin
    ion = filter(v -> real(v) > EPS1 && imag(v) > 0, vals)
    ele = filter(v -> real(v) > EPS1 && imag(v) <= 0, vals)
    (isempty(ion) ? nothing : ion[argmax(real.(ion))],
     isempty(ele) ? nothing : ele[argmax(real.(ele))])
end

nbad = 0
for p in 1:P
    inw = filter(v -> TJLF._contour_inside(v, cfg), refs[p])
    unst = filter(v -> real(v) > EPS1, inw)
    worst = isempty(unst) ? 0.0 :
            maximum(minimum(abs(l - g) for g in vals[p]; init=Inf) for l in unst)
    spur = count(g -> minimum(abs(l - g) for l in refs[p]; init=Inf) > 1e-6, vals[p])
    bi, be = branch(refs[p]); ci, ce = branch(vals[p])
    di = bi === nothing ? NaN : (ci === nothing ? Inf : abs(bi - ci))
    de = be === nothing ? NaN : (ce === nothing ? Inf : abs(be - ce))
    ok = (worst < 1e-6 || fl[p]) && spur == 0 &&
         (!(di isa Float64) || isnan(di) || di < 1e-6 || fl[p]) &&
         (isnan(de) || de < 1e-6 || fl[p])
    ok || (global nbad += 1)
    @printf("p%d: inwin=%3d found=%3d rank=%3d fl=%d  worst_unstable|dl|=%.3g spur=%d  dlead_ion=%.3g ele=%.3g %s\n",
            p, length(inw), length(vals[p]), rk[p], fl[p], worst, spur, di, de,
            ok ? "" : "<-- BAD")
end
println(nbad == 0 ? "CONTOUR_CPU_REAL_OK" : "CONTOUR_CPU_REAL_BAD n=$nbad")
