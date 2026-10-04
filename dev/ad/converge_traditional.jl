# Convergence study of the TRADITIONAL kwscale_scan marginal sfmin.
#
# Question: is the production sfmin (nfactor=8, nefwid=8, nkyhat=4, k_max=4)
# grid-converged, or a coarse-grid artifact of the 8-point factor grid + the
# f_guess linear-extrapolation onset heuristic?  The fine factor scans show the
# clean AE band-entry onset at the marked (kyhat,width) is ~0.008-0.013, well
# below the production sfmin=0.0195 — so this sweep checks where the traditional
# definition itself lands as its grids refine.
#
# We call kwscale_scan directly (mainsub just forwards to it), varying the grid
# kwargs.  sfmin is returned as inputsEP.FACTOR_IN.
#
#   module load julia/1.11.7 && export JULIA_DEPOT_PATH=$PSCRATCH/.julia
#   julia -t 8 --project=. dev/ad/converge_traditional.jl

using TJLFEP
using Printf

const CASE = joinpath(@__DIR__, "..", "..", "examples", "DIIID_202017C42_500ms_v3.1")
const IR = 38

function run_one(opts0, prof, nf, ne, nk, km)
    ep = deepcopy(opts0)
    res = TJLFEP.kwscale_scan(ep, deepcopy(prof), false; inner = :threads,
                              nfactor = nf, nefwid = ne, nkyhat = nk, k_max = km)
    epout = res[2]
    nsolve = km * nk * ne * nf
    return (sfmin = epout.FACTOR_IN, kymark = epout.KYMARK, width = epout.WIDTH_IN, nsolve = nsolve)
end

function main()
    @printf("threads = %d\n", Threads.nthreads())
    opts0, prof, _ = preprocess_gacode_inputs(joinpath(CASE, "input.gacode"),
                                              joinpath(CASE, "input.TGLFEP"))
    opts0.IR = IR
    opts0.N_BASIS = 6
    @printf("N_BASIS = %d  IR = %d\n", opts0.N_BASIS, IR)

    print("warming up (compile)... "); flush(stdout)
    run_one(opts0, prof, 8, 8, 4, 4)
    println("done")

    # (nfactor, nefwid, nkyhat, k_max)
    configs = [
        (8,  8, 4, 4),   # production baseline
        (16, 8, 4, 4),   # refine factor axis
        (32, 8, 4, 4),
        (8,  8, 4, 6),   # more zoom rounds
        (16, 8, 4, 6),
        (32, 8, 4, 8),   # fine factor + deep zoom
        (16, 16, 8, 6),  # also refine (width, kyhat)
    ]

    println("\n  nfactor nefwid nkyhat k_max     sfmin        kymark   width     nsolve   wall(s)")
    for (nf, ne, nk, km) in configs
        t = @elapsed r = run_one(opts0, prof, nf, ne, nk, km)
        @printf("  %5d  %5d  %5d  %4d   %.6e   %.4f   %.4f   %7d   %6.1f\n",
                nf, ne, nk, km, r.sfmin, r.kymark, r.width, r.nsolve, t)
        flush(stdout)
    end
end

main()
