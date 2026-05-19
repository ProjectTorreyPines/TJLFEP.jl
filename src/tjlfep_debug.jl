# Gated debug logging for Julia parity checks (env TJLFEP_DEBUG=1).

const TJLFEP_DEBUG = get(ENV, "TJLFEP_DEBUG", "0") == "1"

"""Print tagged debug line and flush (for Slurm logs)."""
function dbgmsg(args...)
    TJLFEP_DEBUG || return nothing
    print("[TJLFEP_DBG] ")
    println(args...)
    flush(stdout)
    flush(stderr)
    return nothing
end

"""After TJLF_map: mirror Fortran TGLFEP_tglf_map debug fields."""
function debug_dump_tglf_map(inputsEP, inputsPR, inputTJLF)
    TJLFEP_DEBUG || return nothing
    ir = inputsEP.IR
    is = inputsEP.IS_EP + 1
    r_over_a = inputsPR.RMIN[ir] / inputsPR.RMIN[end]
    dbgmsg("tglf_map ir=", ir, " IS_EP=", inputsEP.IS_EP, " is=", is,
        " FACTOR_IN=", inputsEP.FACTOR_IN, " RMIN_LOC=", inputTJLF.RMIN_LOC,
        " KY=", inputTJLF.KY, " KYHAT_IN=", inputsEP.KYHAT_IN,
        " NBASIS_MAX=", inputTJLF.NBASIS_MAX, " NMODES=", inputTJLF.NMODES,
        " AS_e=", inputTJLF.AS[1], " AS_EP=", inputTJLF.AS[is],
        " ZS_e=", inputTJLF.ZS[1], " ZS_EP=", inputTJLF.ZS[is],
        " gamma_thresh=", inputsEP.GAMMA_THRESH,
        " r_over_a=", r_over_a)
    return nothing
end

"""After TJLF.run in TJLFEP_ky: mirror Fortran post-tglf_run debug."""
function debug_dump_ky_postrun(inputsEP, inputTJLF, gamma_out, freq_out, n::Int)
    TJLFEP_DEBUG || return nothing
    n <= length(gamma_out) || return nothing
    dbgmsg("ky_postrun n=", n, " gamma=", gamma_out[n], " freq=", freq_out[n],
        " lkeep=", freq_out[n] < inputsEP.FREQ_AE_UPPER && gamma_out[n] > inputsEP.GAMMA_THRESH)
    return nothing
end

"""First kwscale combo only."""
function debug_dump_kw_combo(inputsEP, i::Int)
    TJLFEP_DEBUG || return nothing
    i == 1 || return nothing
    dbgmsg("kwscale i=", i, " factor_in=", inputsEP.FACTOR_IN,
        " kyhat_in=", inputsEP.KYHAT_IN, " width_in=", inputsEP.WIDTH_IN)
    return nothing
end
