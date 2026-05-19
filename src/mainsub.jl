function mainsub(inputsEP::Options, inputsPR::profile, printout::Bool = true; use_gpu::Bool = false)
    x = inputsEP.PROCESS_IN
    if (x == 1)
        msg = "No"
        return msg
    elseif (x == 2)
        msg = "No"
        return msg
    elseif (x == 3)
        msg = "No"
        return msg
    elseif (x == 4)
        msg = "No"
        return msg
    elseif (x == 5)
        inputsEP.WIDTH_IN_FLAG = false
        inputsEP.MODE_IN = 2
        inputsEP.KY_MODEL = 3
        dbgmsg("mainsub ir=", inputsEP.IR, " suffix=", inputsEP.SUFFIX,
            " SCAN_N=", inputsEP.SCAN_N, " N_BASIS=", inputsEP.N_BASIS)

        growthrate, inputsEP, inputsPR, scalefactor_buffer, wavebuffer_all = kwscale_scan(inputsEP, inputsPR, printout; use_gpu=use_gpu)
        return (growthrate, inputsEP, inputsPR), (scalefactor_buffer, wavebuffer_all)
    elseif (x == 6)
        msg = "No"
        return msg 
    end
end