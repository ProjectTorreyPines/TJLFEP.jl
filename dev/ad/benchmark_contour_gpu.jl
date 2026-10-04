#!/usr/bin/env julia
# Contour-integral (Beyn) batched eigensolver vs fixed-shift SI on REAL harvested pencils.
#
# The fixed 13-shift SI missed 333/503 IR101 ion leaders (coverage-by-sampling); the adaptive
# union salvage still left ~8 misses (modes absent from the calibration sample). The contour
# solver covers the whole unstable window GEOMETRICALLY (every eigenvalue inside the ellipse),
# needs no calibration geev, and flags its own failures (rank saturation / gray residuals) for a
# per-pencil dense fallback. The metric that matters for grid exactness is therefore
#   "leader misses on UNFLAGGED pencils"  (silent failures — must be 0),
# with the flagged fraction being a cost (dense redo), not a correctness, number.
#
# Usage: julia --project=. -t 32 dev/ad/benchmark_contour_gpu.jl
# Env: PENCILS (dir), NGPU(1), NMODES(4), EPS1(1e-12), RUN_SI(1: also run fixed-shift SI),
#      window: RE_LO(-0.02) RE_HI(0.8) IM_MAX(2.6),
#      solver: N_LONG(18) N_SHORT(6) L(64) K(3) RANK_TOL(1e-8) RESID_TOL(1e-7)
#              REFINE(2) CERT_RE(-0.005) MAX_MOVE(0.05) SAT_FRAC(0.9)
using Pkg
Pkg.activate(normpath(@__DIR__, "..", ".."))
import CUDA
using TJLF
using LinearAlgebra, Serialization, Printf, Statistics
import LinearAlgebra.LAPACK: gesv!, geev!

const PDIR   = get(ENV, "PENCILS", normpath(@__DIR__, "pencils_nb16"))
const NGPU   = parse(Int, get(ENV, "NGPU", "1"))
const NMODES = parse(Int, get(ENV, "NMODES", "4"))
const EPS1   = parse(Float64, get(ENV, "EPS1", "1e-12"))
const RUN_SI = parse(Int, get(ENV, "RUN_SI", "1")) != 0

const CFG = TJLF.ContourConfig(
    re_lo     = parse(Float64, get(ENV, "RE_LO", "-0.02")),
    re_hi     = parse(Float64, get(ENV, "RE_HI", "0.8")),
    im_max    = parse(Float64, get(ENV, "IM_MAX", "2.6")),
    n_long    = parse(Int, get(ENV, "N_LONG", "18")),
    n_short   = parse(Int, get(ENV, "N_SHORT", "6")),
    L         = parse(Int, get(ENV, "L", "64")),
    K         = parse(Int, get(ENV, "K", "3")),
    rank_tol  = parse(Float64, get(ENV, "RANK_TOL", "1e-8")),
    resid_tol = parse(Float64, get(ENV, "RESID_TOL", "1e-7")),
    refine    = parse(Int, get(ENV, "REFINE", "2")),
    cert_re   = parse(Float64, get(ENV, "CERT_RE", "-0.005")),
    max_move  = parse(Float64, get(ENV, "MAX_MOVE", "0.05")),
    sat_frac  = parse(Float64, get(ENV, "SAT_FRAC", "0.9")))

topk(vals, k) = sort(vals; by=real, rev=true)[1:min(k, length(vals))]
function branch_leaders(vals)
    ion = filter(v -> real(v) > EPS1 && imag(v) > 0, vals)
    ele = filter(v -> real(v) > EPS1 && imag(v) <= 0, vals)
    (ion = isempty(ion) ? nothing : ion[argmax(real.(ion))],
     ele = isempty(ele) ? nothing : ele[argmax(real.(ele))])
end
run_ref(A, B) = (A2 = copy(A); B2 = copy(B); (A3,_,_) = gesv!(B2, A2); geev!('N','N', A3)[1])

# Shard a batched-solver call over NGPU devices (each call round-robins onto its own device via
# TJLF's _with_device_slot; concurrent calls land on distinct GPUs).
function sharded(f, As, Bs, ngpu)
    P = length(As)
    ngpu = min(ngpu, P)
    ngpu <= 1 && return f(As, Bs)
    bnd = round.(Int, range(0, P; length = ngpu + 1))
    parts = Vector{Any}(undef, ngpu)
    @sync for g in 1:ngpu
        rng = (bnd[g]+1):bnd[g+1]
        Threads.@spawn (parts[g] = f(As[rng], Bs[rng]))
    end
    return parts
end

function score(cands, refs, flagged, label, t, P)
    gerr_ion = Float64[]; gerr_ele = Float64[]; err_top = Float64[]
    n_ion=0; n_ele=0; n_ionmiss=0; n_elemiss=0; n_ionmiss_unfl=0; n_elemiss_unfl=0
    for p in 1:P
        tr_ref = topk(refs[p], NMODES); tr_c = topk(cands[p], NMODES)
        m = min(length(tr_ref), length(tr_c))
        m > 0 && push!(err_top, maximum(abs.(tr_ref[1:m] .- tr_c[1:m])))
        bl = branch_leaders(refs[p]); blc = branch_leaders(cands[p])
        if bl.ion !== nothing
            n_ion += 1
            if blc.ion === nothing
                n_ionmiss += 1; flagged[p] || (n_ionmiss_unfl += 1)
            else
                push!(gerr_ion, abs(real(bl.ion) - real(blc.ion)))
            end
        end
        if bl.ele !== nothing
            n_ele += 1
            if blc.ele === nothing
                n_elemiss += 1; flagged[p] || (n_elemiss_unfl += 1)
            else
                push!(gerr_ele, abs(real(bl.ele) - real(blc.ele)))
            end
        end
    end
    q(v,pr) = isempty(v) ? NaN : quantile(v,pr)
    @printf("\n---- %s : %.1f ms/pencil (wall %.1f s) ----\n", label, 1e3*t/P, t)
    @printf("  top-%d |dLambda|: median %.3g  p90 %.3g  max %.3g\n",
            NMODES, q(err_top,0.5), q(err_top,0.9), isempty(err_top) ? NaN : maximum(err_top))
    @printf("  ion leader : MISSED %d/%d (unflagged: %d)   |dgamma| median %.3g max %.3g\n",
            n_ionmiss, n_ion, n_ionmiss_unfl, q(gerr_ion,0.5),
            isempty(gerr_ion) ? NaN : maximum(gerr_ion))
    @printf("  ele leader : MISSED %d/%d (unflagged: %d)   |dgamma| median %.3g max %.3g\n",
            n_elemiss, n_ele, n_elemiss_unfl, q(gerr_ele,0.5),
            isempty(gerr_ele) ? NaN : maximum(gerr_ele))
    @printf("  flagged pencils (dense fallback): %d/%d\n", count(flagged), P)
    return (n_ionmiss + n_elemiss, n_ionmiss_unfl + n_elemiss_unfl)
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
    @printf("pencils=%d  n=%d  NGPU=%d  devices=%d\n", P, n, NGPU, length(CUDA.devices()))
    nq = 2 * (CFG.n_long + CFG.n_short)
    @printf("window: gamma in [%.3g, %.3g], |freq| <= %.3g   quad nodes=%d  L=%d K=%d (capacity %d)\n",
            CFG.re_lo, CFG.re_hi, CFG.im_max, nq, CFG.L, CFG.K, 2*CFG.K*CFG.L)

    # geev reference for ALL pencils (threaded)
    refs = Vector{Vector{ComplexF64}}(undef, P)
    t_ref = @elapsed (Threads.@threads for p in 1:P; refs[p] = run_ref(As[p], Bs[p]); end)
    @printf("geev ref (threaded, all %d): %.1f s total\n", P, t_ref)

    # Window audit: every mode TJLF consumes (top-NMODES + branch leaders) must be INSIDE the
    # ellipse, otherwise the contour cannot find it BY DESIGN (a window, not solver, failure).
    consumed = ComplexF64[]
    for p in 1:P
        append!(consumed, topk(refs[p], NMODES))
        bl = branch_leaders(refs[p])
        bl.ion === nothing || push!(consumed, bl.ion)
        bl.ele === nothing || push!(consumed, bl.ele)
    end
    unstable = filter(v -> real(v) > EPS1, consumed)
    nout = count(v -> !TJLF._contour_inside(v, CFG), unstable)
    @printf("window audit: consumed unstable modes: gamma in [%.3g, %.3g], |freq| max %.3g -> %d/%d OUTSIDE window\n",
            minimum(real, unstable), maximum(real, unstable), maximum(abs ∘ imag, unstable),
            nout, length(unstable))
    # In-window eigenvalue counts (capacity check against 2K*L).
    nin = [count(v -> TJLF._contour_inside(v, CFG), refs[p]) for p in 1:P]
    @printf("in-window eigenvalue count: median %d  p90 %d  max %d  (capacity %d)\n",
            round(Int, quantile(Float64.(nin), 0.5)), round(Int, quantile(Float64.(nin), 0.9)),
            maximum(nin), 2*CFG.K*CFG.L)

    # warm (JIT + cuBLAS handles) on a small slice
    TJLF.contour_eigvals_batch(As[1:min(P,4)], Bs[1:min(P,4)];
                               cfg = CFG, use_gpu = true, dense_fallback = :none)

    # contour, raw (no fallback): exposes silent misses on unflagged pencils
    t_raw = @elapsed raw = sharded(As, Bs, NGPU) do a, b
        TJLF.contour_eigvals_batch(a, b; cfg = CFG, use_gpu = true, dense_fallback = :none)
    end
    vals_raw = NGPU <= 1 ? raw[1] : reduce(vcat, (r[1] for r in raw))
    flags    = NGPU <= 1 ? raw[2] : reduce(vcat, (r[2] for r in raw))
    ranks    = NGPU <= 1 ? raw[3] : reduce(vcat, (r[3] for r in raw))
    @printf("\nranks: median %d  p90 %d  max %d\n", round(Int, quantile(Float64.(ranks), 0.5)),
            round(Int, quantile(Float64.(ranks), 0.9)), maximum(ranks))
    _, miss_unfl_raw = score(vals_raw, refs, flags, "CONTOUR raw (no fallback)", t_raw, P)

    # contour + per-pencil dense CPU fallback (the production configuration)
    t_fb = @elapsed fb = sharded(As, Bs, NGPU) do a, b
        TJLF.contour_eigvals_batch(a, b; cfg = CFG, use_gpu = true, dense_fallback = :cpu)
    end
    vals_fb = NGPU <= 1 ? fb[1] : reduce(vcat, (r[1] for r in fb))
    flags_fb = NGPU <= 1 ? fb[2] : reduce(vcat, (r[2] for r in fb))
    miss_fb, miss_unfl_fb = score(vals_fb, refs, flags_fb, "CONTOUR + dense fallback", t_fb, P)

    miss_si = -1
    if RUN_SI
        t_si = @elapsed si = sharded(As, Bs, NGPU) do a, b
            TJLF.si_eigvals_batch(a, b; cfg = TJLF.SIConfig(), use_gpu = true)
        end
        vals_si = NGPU <= 1 ? si : reduce(vcat, si)
        miss_si, _ = score(vals_si, refs, falses(P), "FIXED-SHIFT SI (13 shifts)", t_si, P)
    end

    @printf("\n================ verdict ================\n")
    @printf("  leader misses: contour+fallback=%d (silent/unflagged: raw=%d, fb=%d)",
            miss_fb, miss_unfl_raw, miss_unfl_fb)
    RUN_SI && @printf("   fixed-SI=%d", miss_si)
    @printf("\n  wall: contour+fallback=%.1fs (%.1f ms/pencil, NGPU=%d)   geev ref=%.1fs\n",
            t_fb, 1e3*t_fb/P, NGPU, t_ref)
    @printf("  speedup vs threaded geev: %.1fx\n", t_ref / t_fb)
    println(miss_unfl_fb == 0 ? "CONTOUR_SILENT_MISSES_ZERO" : "CONTOUR_HAS_SILENT_MISSES")
end

main()
