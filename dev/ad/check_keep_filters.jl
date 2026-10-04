# Decisive diagnostic for the AD-vs-traditional sfmin gap.
#
# The ae_band (frequency-only) AD path finds sfmin=0.0114 at (kyhat=0.25,
# width=1.571); the traditional kwscale_scan finds sfmin=0.0195 at (kymark=0.297,
# width=1.571). Two possible causes:
#   (1) the SECONDARY keep filters (tearing / ion-,electron-,thermal-,EP-pinch /
#       QL-ratio / θ²) reject the low-onset mode the ae_band path keeps — these
#       need the wavefunction + QL weights (IFLUX=true), OR
#   (2) only the (kyhat,width) grid resolution differs (AD uses 0.25; traditional
#       refines to 0.297).
#
# This runs the *full* Float64 TJLFEP_ky (same routine the traditional scan uses)
# at the relevant operating points and dumps, per mode: γ, freq, in-AE-band?, and
# every keep/reject flag, plus the active REJECT_*_FLAG settings. If a mode is
# in-band & growing but LKEEP=false, a secondary filter binds → cause (1).
#
#   module load julia/1.11.7 && export JULIA_DEPOT_PATH=$PSCRATCH/.julia
#   julia -t 8 --project=. dev/ad/check_keep_filters.jl

using TJLFEP
using Printf

const CASE = joinpath(@__DIR__, "..", "..", "examples", "DIIID_202017C42_500ms_v3.1")
const IR = 38

function dump_point(opts0, prof, ky, w, factors)
    fu = -abs(prof.omegaGAM[IR])
    @printf("\n===== kyhat=%.3f width=%.3f  (FREQ_AE_UPPER=%.4e) =====\n", ky, w, fu)
    for f in factors
        ep = deepcopy(opts0)
        ep.IR = IR; ep.KYHAT_IN = ky; ep.WIDTH_IN = w; ep.FACTOR_IN = f
        ep.FREQ_AE_UPPER = fu
        ep.GAMMA_THRESH = 1.0e-7
        ep.GAMMA_THRESH_MAX = 1.0e-7
        g_out, f_out, _ = TJLFEP.TJLFEP_ky(ep, prof, "", 0)
        NM = ep.NMODES
        @printf("  factor=%.4e\n", f)
        for n in 1:NM
            inband = f_out[n] < fu
            @printf("    n=%d  γ=%+.4e  freq=%+.4e  inAEband=%s  LKEEP=%s | TEAR=%s Pi=%s Pe=%s Pth=%s PEP=%s QLR=%s TH2=%s\n",
                    n, g_out[n], f_out[n], inband ? "Y" : "n",
                    ep.LKEEP[n] ? "Y" : "n",
                    ep.LTEARING[n] ? "Y" : "n", ep.L_I_PINCH[n] ? "Y" : "n",
                    ep.L_E_PINCH[n] ? "Y" : "n", ep.L_TH_PINCH[n] ? "Y" : "n",
                    ep.L_EP_PINCH[n] ? "Y" : "n", ep.L_QL_RATIO[n] ? "Y" : "n",
                    ep.L_THETA_SQ[n] ? "Y" : "n")
        end
    end
end

function main()
    @printf("threads = %d\n", Threads.nthreads())
    opts0, prof, _ = preprocess_gacode_inputs(joinpath(CASE, "input.gacode"),
                                              joinpath(CASE, "input.TGLFEP"))
    opts0.IR = IR
    opts0.N_BASIS = 6
    @printf("N_BASIS = %d\n", opts0.N_BASIS)
    @printf("REJECT flags: TEAR=%s I_PINCH=%s E_PINCH=%s TH_PINCH=%s EP_PINCH=%s | QL_RATIO_THRESH=%g THETA_SQ_THRESH=%g\n",
            opts0.REJECT_TEARING_FLAG, opts0.REJECT_I_PINCH_FLAG, opts0.REJECT_E_PINCH_FLAG,
            opts0.REJECT_TH_PINCH_FLAG, opts0.REJECT_EP_PINCH_FLAG,
            opts0.QL_RATIO_THRESH, opts0.THETA_SQ_THRESH)

    # Fine low-γ sweep through the band-entry to see where LKEEP first turns on
    # and which filter (if any) gates it below the traditional onset.
    dump_point(opts0, prof, 0.25, 1.571, [0.012, 0.014, 0.016, 0.018, 0.020, 0.023, 0.026, 0.030, 0.040])
    dump_point(opts0, prof, 0.297, 1.571, [0.0070, 0.0085, 0.010, 0.012, 0.014, 0.016, 0.0180, 0.0195, 0.022, 0.026])
end

main()
