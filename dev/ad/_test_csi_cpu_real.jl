#!/usr/bin/env julia
# CPU validation of the coverage-certified adaptive SI solver on REAL harvested pencils:
# every unstable in-window mode must be recovered (or the pencil flagged), no spurious values,
# branch leaders must match dense geev. Env: PENCILS (default pencils_nb16), NP (default 8).
using Pkg
Pkg.activate(normpath(@__DIR__, "..", ".."))
using TJLF
using LinearAlgebra, Serialization, Printf
import LinearAlgebra.LAPACK: gesv!, geev!

# Julia threads drive the pencil-parallel loops; keep BLAS single-threaded to avoid
# oversubscribing the cores (32 Julia threads × an OpenBLAS pool busy-wait otherwise).
LinearAlgebra.BLAS.set_num_threads(1)

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

cfg = TJLF.CertifiedSIConfig()
inwin(v) = 0.0 ≤ real(v) ≤ cfg.re_hi && abs(imag(v)) ≤ cfg.im_max
println("unstable in-window counts: ",
        [count(v -> inwin(v) && real(v) > EPS1, r) for r in refs])

t = @elapsed vals, fl, ns, rs = TJLF.certified_si_eigvals_batch(As, Bs; cfg, use_gpu=false,
                                                                dense_fallback=:none)
@printf("cpu certified-SI: %.1f s total (%.1f s/pencil)\nnshifts=%s\nflagged=%s\nreasons=%s\n",
        t, t / P, ns, fl, rs)

branch(vals) = begin
    ion = filter(v -> real(v) > EPS1 && imag(v) > 0, vals)
    ele = filter(v -> real(v) > EPS1 && imag(v) <= 0, vals)
    (isempty(ion) ? nothing : ion[argmax(real.(ion))],
     isempty(ele) ? nothing : ele[argmax(real.(ele))])
end

nbad = 0
for p in 1:P
    unst = filter(v -> inwin(v) && real(v) > EPS1, refs[p])
    worst = isempty(unst) ? 0.0 :
            maximum(minimum(abs(l - g) for g in vals[p]; init=Inf) for l in unst)
    spur = count(g -> minimum(abs(l - g) for l in refs[p]; init=Inf) > 1e-6, vals[p])
    bi, be = branch(refs[p]); ci, ce = branch(vals[p])
    di = bi === nothing ? NaN : (ci === nothing ? Inf : abs(bi - ci))
    de = be === nothing ? NaN : (ce === nothing ? Inf : abs(be - ce))
    ok = (worst < 1e-6 || fl[p]) && spur == 0 &&
         (isnan(di) || di < 1e-6 || fl[p]) && (isnan(de) || de < 1e-6 || fl[p])
    ok || (global nbad += 1)
    @printf("p%d: unst=%2d found=%3d shifts=%2d fl=%d  worst|dl|=%.3g spur=%d  dlead_ion=%.3g ele=%.3g %s\n",
            p, length(unst), length(vals[p]), ns[p], fl[p], worst, spur, di, de,
            ok ? "" : "<-- BAD")
end
println(nbad == 0 ? "CSI_CPU_REAL_OK" : "CSI_CPU_REAL_BAD n=$nbad")
