# Per-solve cost: one Float64 TGLF eigensolve (TJLFEP_ky path) vs one
# forward-mode-AD solve (gamma_dgamma_dfactor) at the same operating point, at
# both nb=6 and nb=32. The Dual/Float64 ratio is what converts the AD
# eigensolve-count advantage (~4× fewer) into a wall-time speed-up.
#
#   module load julia/1.11.7 && export JULIA_DEPOT_PATH=$PSCRATCH/.julia
#   julia --project=. dev/ad/micro_benchmark_solve.jl

using TJLFEP
using TJLF
using Printf

const CASE = joinpath(@__DIR__, "..", "..", "examples", "DIIID_202017C42_500ms_v3.1")
const IR, KYHAT_IN, WIDTH_IN, FACTOR = 38, 0.25, 1.0, 2.5

function bench_nb(opts0, prof, nb)
    opts = deepcopy(opts0)
    opts.N_BASIS = nb
    opts.IR = IR; opts.KYHAT_IN = KYHAT_IN; opts.WIDTH_IN = WIDTH_IN; opts.FACTOR_IN = FACTOR

    # Float64: one TGLF eigensolve via the production single-ky path.
    f64() = TJLFEP.TJLFEP_ky(deepcopy(opts), prof, "", 0)
    # Dual: value + dγ/dfactor in one forward pass.
    dual() = gamma_dgamma_dfactor(opts, prof)

    f64(); dual()                       # warm up (compile)
    t_f64 = minimum(@elapsed(f64()) for _ in 1:3)
    t_dual = minimum(@elapsed(dual()) for _ in 1:3)
    @printf("nb=%2d  matrix=%4d  t_f64=%7.3f s  t_dual=%7.3f s  ratio=%.2f×\n",
            nb, 3 * 15 * nb, t_f64, t_dual, t_dual / t_f64)
    return t_f64, t_dual
end

function main()
    opts0, prof, _ = preprocess_gacode_inputs(joinpath(CASE, "input.gacode"),
                                              joinpath(CASE, "input.TGLFEP"))
    println("Per-solve Float64 vs Dual (forward-mode AD):")
    bench_nb(opts0, prof, 6)
    f64_32, dual_32 = bench_nb(opts0, prof, 32)

    # Projected full-grid wall-time speed-up at nb=32 using the measured ratio
    # and the observed eigensolve counts (trad 1024 vs AD 257).
    ratio = dual_32 / f64_32
    proj = (257 * ratio) / (1024 * 1.0)
    @printf("\nProjection (nb=32, counts 1024 vs 257): AD/trad wall ≈ (257·%.2f)/1024 = %.2f → %.2f× %s\n",
            ratio, proj, 1 / proj, proj < 1 ? "faster" : "slower")
end

main()
