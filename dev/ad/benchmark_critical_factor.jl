# Head-to-head: AD-Newton critical-factor vs the traditional Fortran-style
# brute-force (kyhat,width,factor) grid scan, on one DIII-D radius.
#
# Traditional  : `mainsub` → `kwscale_scan` (PROCESS_IN=5). Fixed cost of
#                k_max(4) × nkyhat(4)·nefwid(8)·nfactor(8) = 1024 TGLF eigensolves,
#                producing sfmin = inputsEP.FACTOR_IN on return.
# AD-Newton    : `critical_factor_grid` — for each of the 32 (kyhat,width) grid
#                points, a safeguarded Newton on γ_lead(factor)=thresh using the
#                exact AD dγ/dfactor, then min over the grid.
#
# Reports sfmin agreement, eigensolve counts (the parallelism-independent compute
# metric), and wall time (both threaded, BLAS=1).
#
#   module load julia/1.11.7 && export JULIA_DEPOT_PATH=$PSCRATCH/.julia
#   julia -t 8 --project=. dev/ad/benchmark_critical_factor.jl

using TJLFEP
using Printf

const CASE = joinpath(@__DIR__, "..", "..", "examples", "DIIID_202017C42_500ms_v3.1")
const IR = 38

function main()
    @printf("threads = %d\n", Threads.nthreads())
    opts0, prof, _ = preprocess_gacode_inputs(joinpath(CASE, "input.gacode"),
                                              joinpath(CASE, "input.TGLFEP"))
    opts0.IR = IR
    # Reduce basis for tractable wall-time timing (eigensolve-COUNT speed-up is
    # basis-independent). nb=6 is a supported config (input_scan20_nb6.TGLFEP).
    opts0.N_BASIS = 6
    @printf("N_BASIS = %d (reduced for timing; count speed-up is basis-independent)\n", opts0.N_BASIS)

    # ── Warm up both paths (compile) on a throwaway copy so timings are clean ──
    println("warming up (compilation)...")
    let ep = deepcopy(opts0)
        TJLFEP.mainsub(ep, deepcopy(prof), false; inner = :threads)
    end
    let ep = deepcopy(opts0)
        critical_factor_grid(ep, prof; threaded = true, ae_band = true)
    end

    # ── Traditional brute-force scan ──
    println("\n[traditional] kwscale_scan (1024 eigensolves)...")
    ep_t = deepcopy(opts0)
    t_trad = @elapsed begin
        res_t, _ = TJLFEP.mainsub(ep_t, deepcopy(prof), false; inner = :threads)
    end
    sfmin_trad = res_t[2].FACTOR_IN
    kymark_trad = res_t[2].KYMARK
    width_trad = res_t[2].WIDTH_IN
    n_trad = 4 * 4 * 8 * 8   # k_max × nkyhat × nefwid × nfactor

    # ── AD-Newton grid ──
    println("[AD-Newton] critical_factor_grid...")
    ep_a = deepcopy(opts0)
    local res_a
    t_ad = @elapsed begin
        res_a = critical_factor_grid(ep_a, prof; threaded = true, ae_band = true)
    end

    # ── Report ──
    println("\n================ RESULTS (DIII-D IR=$IR) ================")
    @printf("traditional sfmin = %.6e   (kymark=%.3f width=%.3f)\n", sfmin_trad, kymark_trad, width_trad)
    @printf("AD-Newton  sfmin = %.6e   (kyhat=%.3f width=%.3f)\n", res_a.sfmin, res_a.kyhat, res_a.width)
    @printf("sfmin ratio AD/trad = %.4f\n", res_a.sfmin / sfmin_trad)

    println("\n---- cost ----")
    @printf("eigensolves: traditional = %d   AD-Newton = %d   → %.1f× fewer\n",
            n_trad, res_a.total_evals, n_trad / res_a.total_evals)
    @printf("wall time:   traditional = %.1f s   AD-Newton = %.1f s   → %.2f× %s\n",
            t_trad, t_ad, t_trad / t_ad, t_trad > t_ad ? "faster" : "slower")
    @printf("avg eigensolve: traditional = %.2f s/solve   AD-Newton = %.2f s/solve (Dual)\n",
            t_trad / n_trad, t_ad / max(res_a.total_evals, 1))

    println("\n---- AD per-grid-point marginal factors ----")
    for o in sort(res_a.results; by = x -> x.factor)
        @printf("  kyhat=%.3f width=%.3f  f★=%11.4e  evals=%d  %s\n",
                o.kyhat, o.width, o.factor, o.evals, o.converged ? "" : "(no crossing)")
    end
end

main()
