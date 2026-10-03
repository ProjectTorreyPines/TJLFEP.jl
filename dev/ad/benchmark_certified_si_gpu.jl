#!/usr/bin/env julia
# Coverage-certified adaptive SI vs fixed-shift SI on REAL harvested pencils.
#
# History on IR101 (nb16): fixed 13-shift SI missed 333/503 ion leaders; geev-calibrated
# adaptive-union still missed ~8/1024; the contour (SS-RR) solver was safe but saturated on the
# axis-hugging mode crowd (everything flagged -> all-dense). The certified solver keeps the SI
# kernel but (a) residual-certifies every returned value against the original pencil and
# (b) requires a per-pencil geometric COVERAGE certificate of the unstable window, adding
# per-pencil shifts (batched rounds) at uncovered spots. The correctness metric is
#   "leader misses on UNFLAGGED pencils"  (silent failures — must be 0);
# flagged pencils are a cost (dense redo), not a correctness issue.
#
# Usage: julia --project=. -t 32 dev/ad/benchmark_certified_si_gpu.jl
# Env: PENCILS (dir), NGPU(1), NMODES(4), EPS1(1e-12), RUN_SI(1: also run fixed-shift SI);
#      solver overrides (defaults = CertifiedSIConfig() struct defaults): RE_HI IM_MAX ROW_RE
#      ROW_DY ROW_DENSE_W M Q RESID_TOL ENUM_TOL TRUST EMPTY_FRAC REFINE MAX_ROUNDS
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

const DEF = TJLF.CertifiedSIConfig()
envf(k, d) = parse(Float64, get(ENV, k, string(d)))
envi(k, d) = parse(Int, get(ENV, k, string(d)))
const CFG = TJLF.CertifiedSIConfig(
    re_hi       = envf("RE_HI", DEF.re_hi),
    im_max      = envf("IM_MAX", DEF.im_max),
    row_re      = envf("ROW_RE", DEF.row_re),
    row_dy      = envf("ROW_DY", DEF.row_dy),
    row_dense_w = envf("ROW_DENSE_W", DEF.row_dense_w),
    M           = envi("M", DEF.M),
    Q           = envi("Q", DEF.Q),
    resid_tol   = envf("RESID_TOL", DEF.resid_tol),
    enum_tol    = envf("ENUM_TOL", DEF.enum_tol),
    trust       = envi("TRUST", DEF.trust),
    empty_frac  = envf("EMPTY_FRAC", DEF.empty_frac),
    refine      = envi("REFINE", DEF.refine),
    max_rounds  = envi("MAX_ROUNDS", DEF.max_rounds))

topk(vals, k) = sort(vals; by=real, rev=true)[1:min(k, length(vals))]
function branch_leaders(vals)
    ion = filter(v -> real(v) > EPS1 && imag(v) > 0, vals)
    ele = filter(v -> real(v) > EPS1 && imag(v) <= 0, vals)
    (ion = isempty(ion) ? nothing : ion[argmax(real.(ion))],
     ele = isempty(ele) ? nothing : ele[argmax(real.(ele))])
end
run_ref(A, B) = (A2 = copy(A); B2 = copy(B); (A3,_,_) = gesv!(B2, A2); geev!('N','N', A3)[1])

# Shard a batched-solver call over NGPU devices (each concurrent call pins its own device).
function sharded(f, As, Bs, ngpu)
    P = length(As)
    ngpu = min(ngpu, P)
    ngpu <= 1 && return [f(As, Bs)]
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
            if blc.ion === nothing || abs(bl.ion - blc.ion) > 1e-6
                n_ionmiss += 1; flagged[p] || (n_ionmiss_unfl += 1)
            else
                push!(gerr_ion, abs(real(bl.ion) - real(blc.ion)))
            end
        end
        if bl.ele !== nothing
            n_ele += 1
            if blc.ele === nothing || abs(bl.ele - blc.ele) > 1e-6
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
    @printf("  ion leader : MISSED/WRONG %d/%d (unflagged: %d)   |dgamma| median %.3g max %.3g\n",
            n_ionmiss, n_ion, n_ionmiss_unfl, q(gerr_ion,0.5),
            isempty(gerr_ion) ? NaN : maximum(gerr_ion))
    @printf("  ele leader : MISSED/WRONG %d/%d (unflagged: %d)   |dgamma| median %.3g max %.3g\n",
            n_elemiss, n_ele, n_elemiss_unfl, q(gerr_ele,0.5),
            isempty(gerr_ele) ? NaN : maximum(gerr_ele))
    @printf("  flagged pencils (dense fallback): %d/%d\n", count(flagged), P)
    return (n_ionmiss + n_elemiss, n_ionmiss_unfl + n_elemiss_unfl)
end

function main()
    # Julia threads drive the pencil-parallel loops (geev ref, finalize, dense fallback); keep
    # BLAS single-threaded per pencil, else 32 Julia threads × an OpenBLAS pool oversubscribe
    # 64 cores and busy-wait in blas_thread_server (observed: 95 threads at 99% CPU, GPUs idle).
    LinearAlgebra.BLAS.set_num_threads(1)
    @printf("BLAS threads=%d  Julia threads=%d\n", LinearAlgebra.BLAS.get_num_threads(),
            Threads.nthreads()); flush(stdout)
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
    @printf("window: gamma in [0, %.3g], |freq| <= %.3g   M=%d Q=%d trust=%d max_rounds=%d\n",
            CFG.re_hi, CFG.im_max, CFG.M, CFG.Q, CFG.trust, CFG.max_rounds); flush(stdout)

    refs = Vector{Vector{ComplexF64}}(undef, P)
    t_ref = @elapsed (Threads.@threads for p in 1:P; refs[p] = run_ref(As[p], Bs[p]); end)
    @printf("geev ref (threaded, all %d): %.1f s total\n", P, t_ref); flush(stdout)

    # Window audit: consumed modes outside [0,re_hi]x[+-im_max] are invisible BY DESIGN.
    consumed = ComplexF64[]
    for p in 1:P
        append!(consumed, topk(refs[p], NMODES))
        bl = branch_leaders(refs[p])
        bl.ion === nothing || push!(consumed, bl.ion)
        bl.ele === nothing || push!(consumed, bl.ele)
    end
    unstable = filter(v -> real(v) > EPS1, consumed)
    nout = count(v -> !(0.0 <= real(v) <= CFG.re_hi && abs(imag(v)) <= CFG.im_max), unstable)
    @printf("window audit: consumed unstable modes: gamma in [%.3g, %.3g], |freq| max %.3g -> %d/%d OUTSIDE window\n",
            minimum(real, unstable), maximum(real, unstable), maximum(abs ∘ imag, unstable),
            nout, length(unstable)); flush(stdout)

    # warm (JIT + cuBLAS handles) on a small slice
    TJLF.certified_si_eigvals_batch(As[1:min(P,4)], Bs[1:min(P,4)];
                                    cfg = CFG, use_gpu = true, dense_fallback = :none)
    println("warm-up done"); flush(stdout)

    # raw (no fallback): exposes silent misses on unflagged pencils
    t_raw = @elapsed raw = sharded(As, Bs, NGPU) do a, b
        TJLF.certified_si_eigvals_batch(a, b; cfg = CFG, use_gpu = true, dense_fallback = :none)
    end
    vals_raw = reduce(vcat, (r[1] for r in raw))
    flags    = reduce(vcat, (r[2] for r in raw))
    nshifts  = reduce(vcat, (r[3] for r in raw))
    reasons  = reduce(vcat, (r[4] for r in raw))
    @printf("\nshifts/pencil: median %d  p90 %d  max %d\n",
            round(Int, quantile(Float64.(nshifts), 0.5)),
            round(Int, quantile(Float64.(nshifts), 0.9)), maximum(nshifts))
    @printf("flag reasons: ok=%d uncovered=%d uncert=%d both=%d\n",
            count(==(:ok), reasons), count(==(:uncovered), reasons),
            count(==(:uncert), reasons), count(==(:both), reasons))
    _, miss_unfl_raw = score(vals_raw, refs, flags, "CERTIFIED-SI raw (no fallback)", t_raw, P)

    # + per-pencil dense CPU fallback (the production configuration)
    t_fb = @elapsed fb = sharded(As, Bs, NGPU) do a, b
        TJLF.certified_si_eigvals_batch(a, b; cfg = CFG, use_gpu = true, dense_fallback = :cpu)
    end
    vals_fb  = reduce(vcat, (r[1] for r in fb))
    flags_fb = reduce(vcat, (r[2] for r in fb))
    miss_fb, miss_unfl_fb = score(vals_fb, refs, flags_fb, "CERTIFIED-SI + dense fallback",
                                  t_fb, P)

    miss_si = -1
    if RUN_SI
        t_si = @elapsed si = sharded(As, Bs, NGPU) do a, b
            TJLF.si_eigvals_batch(a, b; cfg = TJLF.SIConfig(), use_gpu = true)
        end
        vals_si = reduce(vcat, si)
        miss_si, _ = score(vals_si, refs, falses(P), "FIXED-SHIFT SI (13 shifts)", t_si, P)
    end

    @printf("\n================ verdict ================\n")
    @printf("  leader misses: certified+fallback=%d (silent/unflagged: raw=%d, fb=%d)",
            miss_fb, miss_unfl_raw, miss_unfl_fb)
    RUN_SI && @printf("   fixed-SI=%d", miss_si)
    @printf("\n  wall: certified+fallback=%.1fs (%.1f ms/pencil, NGPU=%d)   geev ref=%.1fs\n",
            t_fb, 1e3*t_fb/P, NGPU, t_ref)
    @printf("  speedup vs threaded geev: %.1fx\n", t_ref / t_fb)
    println(miss_unfl_fb == 0 ? "CSI_SILENT_MISSES_ZERO" : "CSI_HAS_SILENT_MISSES")
end

main()
