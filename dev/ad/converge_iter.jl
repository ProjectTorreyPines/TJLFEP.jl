# Convergence study of the TRADITIONAL kwscale_scan marginal sfmin on the ITER
# case, which has ROTATIONAL_SUPPRESSION_FLAG=1  ->  GAMMA_THRESH = 0.15·|γ_E/ŝ|
# (a FINITE growth-rate threshold).  Hypothesis: unlike the DIII-D case
# (GAMMA_THRESH=1e-7, where sfmin collapses toward 0 under refinement), the ITER
# marginal is a smooth root γ_keep(factor)=γ*>0 and so should be GRID-CONVERGED.
#
# We rebuild Options/profile/expro exactly as runTHD does, then for a few scan
# radii call kwscale_scan directly with increasing factor resolution / zoom.
#
#   module load julia/1.11.7 && export JULIA_DEPOT_PATH=$PSCRATCH/.julia
#   julia -t 8 --project=. dev/ad/converge_iter.jl

using TJLFEP
using Printf

const DIR = joinpath(@__DIR__, "..", "..", "examples", "ITER")

function build_inputs()
    tglfep = joinpath(DIR, "input.TGLFEP")
    mtglf  = joinpath(DIR, "input.MTGLF")
    expro_ = joinpath(DIR, "input.EXPRO")
    prof = TJLFEP.readMTGLF(mtglf)
    profile = prof[1]
    ir_exp = prof[2]
    Options = TJLFEP.readTGLFEP(tglfep, ir_exp)
    ni, Ti, dlnnidr, dlntidr, cs, rmin_ex, gammaE, gammap, omegaGAM =
        TJLFEP.read_expro_for_alpha(expro_, profile, Options.IS_EP; gacode_file = nothing)
    expro = (ni=ni, Ti=Ti, dlnnidr=dlnnidr, dlntidr=dlntidr, cs=cs, rmin_ex=rmin_ex,
             gammaE=gammaE, gammap=gammap, omegaGAM=omegaGAM)
    TJLFEP._apply_runthd_expro_setup!(Options, profile, expro)
    return Options, profile
end

function run_one(opts, prof, i, nf, ne, nk, km)
    ep = deepcopy(opts)
    ep.IR = ep.IR_EXP[i]
    ep.FACTOR_IN = ep.FACTOR[i]
    res = TJLFEP.kwscale_scan(ep, deepcopy(prof), false; inner = :threads,
                              nfactor = nf, nefwid = ne, nkyhat = nk, k_max = km)
    epout = res[2]
    return (sfmin = epout.FACTOR_IN, kymark = epout.KYMARK, width = epout.WIDTH_IN,
            ir = ep.IR, nsolve = km * nk * ne * nf)
end

function main()
    @printf("threads = %d\n", Threads.nthreads())
    opts, prof = build_inputs()
    @printf("N_BASIS=%d  SCAN_N=%d  IRS=%d  ROT_SUPP=%d  IR_EXP=%s\n",
            opts.N_BASIS, opts.SCAN_N, opts.IRS, opts.ROTATIONAL_SUPPRESSION_FLAG, string(opts.IR_EXP))

    print("warming up (compile)... "); flush(stdout)
    run_one(opts, prof, 1, 8, 8, 4, 4)
    println("done")

    configs = [(8, 8, 4, 4), (16, 8, 4, 6), (32, 8, 4, 8)]
    for i in 1:opts.SCAN_N
        @printf("\n--- scan radius i=%d (IR=%d) ---\n", i, opts.IR_EXP[i])
        println("  nfactor k_max     sfmin        kymark   width    nsolve  wall(s)")
        for (nf, ne, nk, km) in configs
            t = @elapsed r = run_one(opts, prof, i, nf, ne, nk, km)
            @printf("  %5d  %4d   %.6e   %.4f   %.4f   %6d   %5.1f\n",
                    nf, km, r.sfmin, r.kymark, r.width, r.nsolve, t)
            flush(stdout)
        end
    end
end

main()
