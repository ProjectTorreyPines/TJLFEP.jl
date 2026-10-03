using Pkg; Pkg.activate(normpath(@__DIR__, "..", ".."))
import CUDA
using TJLFEP, TJLF, Printf
using TJLFEP: preprocess_gacode_inputs, kwscale_scan

# End-to-end: does inner=:batched_si (batched shift-invert eigensolver) reproduce the exact grid
# answer (inner=:threads, dense geev) for the same grid resolution? Compares marginal sfmin / ky /
# width per radius. :batched_si uses the GPU batched solver; :threads is the bit-exact CPU golden.
CASE = normpath(@__DIR__, "..", "..", "examples", "DIIID_202017C42_500ms_v3.1")
NB   = parse(Int, get(ENV, "NB", "16"))
GAC  = joinpath(CASE, "input.gacode")
TGL  = joinpath(CASE, "input_scan20_nb$(NB).TGLFEP")
USE_GPU = get(ENV, "USE_GPU", "1") == "1"

# Moderate grid so the dense CPU golden stays tractable at n>=720 (both paths use the same grid,
# so agreement proves eigenvalue equivalence regardless of resolution).
# NOTE: k_max=2 has NO intermediate k-round, so the hybrid never exercises the batched solver —
# use KMAX>=3 (ideally 4) to actually test it (the k_max=2 "E2E_OK" trap, see README_batched_GPU).
grid = (nfactor = parse(Int, get(ENV,"NFACTOR","6")),
        nefwid  = parse(Int, get(ENV,"NEFWID","4")),
        nkyhat  = parse(Int, get(ENV,"NKYHAT","2")),
        k_max   = parse(Int, get(ENV,"KMAX","2")))

# SI_METHOD=si (fixed-shift subspace iteration) or contour (Beyn integral, residual-certified,
# per-pencil dense fallback). SI_DENSE_FALLBACK toggles the per-RADIUS hybrid guard.
SI_METHOD = Symbol(get(ENV, "SI_METHOD", "si"))
SI_DENSE_FALLBACK = get(ENV, "SI_DENSE_FALLBACK", "1") == "1"

base_ep, prof, _ = preprocess_gacode_inputs(GAC, TGL)
scan_n = Int(base_ep.SCAN_N)
radii = let r = get(ENV, "RADII", "")
    isempty(r) ? collect(1:scan_n) : parse.(Int, split(r, ","))
end
@printf("e2e compare nb=%d grid=%s use_gpu=%s si_method=%s dense_fb=%s radii=%s\n",
        NB, grid, USE_GPU, SI_METHOD, SI_DENSE_FALLBACK, join(radii, ","))

function prep(i)
    ep = deepcopy(base_ep); ep.IR = base_ep.IR_EXP[i]
    ep.WIDTH_IN_FLAG = false; ep.MODE_IN = 2; ep.KY_MODEL = 3; ep.PROCESS_IN = 5
    ep.FACTOR_IN = Float64(base_ep.FACTOR[i]); ep
end
function run1(ep, inner, use_gpu)
    t = @elapsed begin
        _g, epo, = kwscale_scan(ep, prof, false; use_gpu=use_gpu, inner=inner,
                                si_method=SI_METHOD, si_dense_fallback=SI_DENSE_FALLBACK, grid...)
    end
    (; sfmin=Float64(epo.FACTOR_IN), ky=Float64(epo.KYMARK), w=Float64(epo.WIDTH_IN), wall=t)
end

# threads golden can also run on the GPU (serial Xgeev) for a same-device comparison; default CPU.
GOLDEN_GPU = get(ENV, "GOLDEN_GPU", "0") == "1"

@printf("%4s | %-30s | %-30s | %s\n", "IR", "threads (golden)", "batched_si (hybrid)", "dsfmin")
nbad = 0; twall_t = 0.0; twall_b = 0.0; dsfs = Float64[]
for i in radii
    ir = base_ep.IR_EXP[i]
    rt = run1(prep(i), :threads, GOLDEN_GPU)
    rb = run1(prep(i), :batched_si, USE_GPU)
    dsf = abs(rb.sfmin - rt.sfmin) / max(abs(rt.sfmin), 1e-9)
    dky = abs(rb.ky - rt.ky); dw = abs(rb.w - rt.w)
    ok = dsf <= 0.02 && dky <= 1e-6 && dw <= 1e-6
    ok || (global nbad += 1)
    global twall_t += rt.wall; global twall_b += rb.wall; push!(dsfs, dsf)
    @printf("%4d | sf=%9.4g ky=%.3f w=%.3f %5.0fs | sf=%9.4g ky=%.3f w=%.3f %5.0fs | %.2e %s\n",
            ir, rt.sfmin, rt.ky, rt.w, rt.wall, rb.sfmin, rb.ky, rb.w, rb.wall, dsf, ok ? "" : "<-- MISMATCH")
    flush(stdout)
end
sd = sort(dsfs)
@printf("\n================ summary over %d radii ================\n", length(radii))
@printf("  wall total: threads=%.0fs  batched_si=%.0fs  speedup=%.2fx\n",
        twall_t, twall_b, twall_t / max(twall_b, 1e-9))
@printf("  dsfmin: median=%.2e  p90=%.2e  max=%.2e   mismatches(>2%%)=%d/%d\n",
        sd[cld(end,2)], sd[max(1,ceil(Int,0.9*end))], maximum(sd), nbad, length(radii))
println(nbad == 0 ? "E2E_OK" : "E2E_MISMATCH n=$nbad")
