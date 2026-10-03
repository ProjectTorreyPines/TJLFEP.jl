# ITER keep-filter scan: does a SECONDARY filter (QL-ratio thresh=0.001, th-pinch)
# gate the marginal onset, pushing it from the frequency band-entry (~2.2) up to
# the brute converged sfmin (~20) at the marked point (kyhat=0.151, w=1.494)?
#
#   module load julia/1.11.7 && export JULIA_DEPOT_PATH=$PSCRATCH/.julia
#   julia -t 8 --project=. dev/ad/check_keep_filters_iter.jl

using TJLFEP
using Printf

const DIR = joinpath(@__DIR__, "..", "..", "examples", "ITER")

function build_inputs()
    prof = TJLFEP.readMTGLF(joinpath(DIR, "input.MTGLF"))
    profile = prof[1]; ir_exp = prof[2]
    Options = TJLFEP.readTGLFEP(joinpath(DIR, "input.TGLFEP"), ir_exp)
    ni, Ti, dlnnidr, dlntidr, cs, rmin_ex, gammaE, gammap, omegaGAM =
        TJLFEP.read_expro_for_alpha(joinpath(DIR, "input.EXPRO"), profile, Options.IS_EP; gacode_file = nothing)
    expro = (ni=ni, Ti=Ti, dlnnidr=dlnnidr, dlntidr=dlntidr, cs=cs, rmin_ex=rmin_ex,
             gammaE=gammaE, gammap=gammap, omegaGAM=omegaGAM)
    TJLFEP._apply_runthd_expro_setup!(Options, profile, expro)
    return Options, profile
end

function dump(opts, prof, ir, ky, w, factors)
    fu = -abs(prof.omegaGAM[ir])
    base = deepcopy(opts); base.IR = ir
    gstar = TJLFEP._gamma_thresh_for(base, prof)
    @printf("\n===== ITER IR=%d kyhat=%.3f width=%.3f  (FREQ_AE_UPPER=%.4e, γ*=%.4e) =====\n", ir, ky, w, fu, gstar)
    @printf("REJECT: TEAR=%s Ipinch=%s Epinch=%s THpinch=%s EPpinch=%s | QL_RATIO_THRESH=%g THETA_SQ_THRESH=%g\n",
            opts.REJECT_TEARING_FLAG, opts.REJECT_I_PINCH_FLAG, opts.REJECT_E_PINCH_FLAG,
            opts.REJECT_TH_PINCH_FLAG, opts.REJECT_EP_PINCH_FLAG, opts.QL_RATIO_THRESH, opts.THETA_SQ_THRESH)
    for f in factors
        ep = deepcopy(base); ep.KYHAT_IN = ky; ep.WIDTH_IN = w; ep.FACTOR_IN = f
        ep.FREQ_AE_UPPER = fu; ep.GAMMA_THRESH = gstar; ep.GAMMA_THRESH_MAX = gstar
        g_out, f_out, _ = TJLFEP.TJLFEP_ky(ep, prof, "", 0)
        @printf("  factor=%.4e\n", f)
        for n in 1:ep.NMODES
            inband = f_out[n] < fu
            @printf("    n=%d  γ=%+.4e  freq=%+.4e  inAE=%s  LKEEP=%s | TEAR=%s Pi=%s Pe=%s Pth=%s PEP=%s QLR=%s TH2=%s\n",
                    n, g_out[n], f_out[n], inband ? "Y" : "n", ep.LKEEP[n] ? "Y" : "n",
                    ep.LTEARING[n] ? "Y" : "n", ep.L_I_PINCH[n] ? "Y" : "n", ep.L_E_PINCH[n] ? "Y" : "n",
                    ep.L_TH_PINCH[n] ? "Y" : "n", ep.L_EP_PINCH[n] ? "Y" : "n", ep.L_QL_RATIO[n] ? "Y" : "n",
                    ep.L_THETA_SQ[n] ? "Y" : "n")
        end
        flush(stdout)
    end
end

function main()
    @printf("threads = %d\n", Threads.nthreads())
    opts, prof = build_inputs()
    @printf("N_BASIS=%d\n", opts.N_BASIS)
    dump(opts, prof, opts.IR_EXP[3], 0.151, 1.494, [1.0, 2.0, 3.0, 5.0, 8.0, 12.0, 16.0, 20.0, 25.0, 30.0])
end

main()
