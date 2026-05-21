# Match Fortran out.scalefactor layout (format 20 / 60 in TGLFEP_kwscale_scan.f90).
function _scalefactor_section_header(kyhat::T, ky_write::T, efwid::T, fmark::T) where {T<:Real}
    return @sprintf(
        "--------------- ky*rho_EP= %7.4f  (ky*rho_s=%7.4f)   width=%7.4f   scalefactor=%10.4f -------------------",
        kyhat, ky_write, efwid, fmark,
    )
end

function _scalefactor_factor_line(fac::T, g::AbstractVector{T}, f::AbstractVector{T}, keep_label::AbstractVector{String}) where {T<:Real}
    parts = [@sprintf("%11.4f", fac)]
    for n in eachindex(keep_label)
        push!(parts, @sprintf(" %14.5E", g[n]), @sprintf(" %14.5E", f[n]), " ", keep_label[n])
    end
    return join(parts)
end

function kwscale_scan(inputsEP::Options{T}, inputsPR::profile{T}, printout::Bool = true; use_gpu::Bool = false) where {T<:Real}
    # These are for testing purposes:
    #baseDirectory = "/Users/benagnew/TJLF.jl/outputs/tglfep_tests/input.MTGLF"
    #inputsPR = readMTGLF(baseDirectory)

    nfactor = 8
    nefwid = 8
    nkyhat = 4
    nkwf = nfactor*nefwid*nkyhat
    k_max = 4
    l_write_out = true

    kyhat_min = 0.0
    kyhat_max = 1.0

    # ikymin = 0
    # ikymax =0
    # iefwmin =0
    # iefwmax =0

    growthrate = zeros(T, nkyhat, nefwid, nfactor, inputsEP.NMODES)
    growthrate_out = zeros(T, nkyhat, nefwid, nfactor, inputsEP.NMODES)

    frequency = zeros(T, nkyhat, nefwid, nfactor, inputsEP.NMODES)
    frequency_out = zeros(T, nkyhat, nefwid, nfactor, inputsEP.NMODES)
    
    lkeep_i = fill(true, (nkyhat, nefwid, nfactor, inputsEP.NMODES))
    ltearing_i = fill(false, (nkyhat, nefwid, nfactor, inputsEP.NMODES))
    l_th_pinch_i = fill(false, (nkyhat, nefwid, nfactor, inputsEP.NMODES))
    l_i_pinch_i = fill(false, (nkyhat, nefwid, nfactor, inputsEP.NMODES))
    l_e_pinch_i = fill(false, (nkyhat, nefwid, nfactor, inputsEP.NMODES))
    l_EP_pinch_i = fill(false, (nkyhat, nefwid, nfactor, inputsEP.NMODES))
    #l_max_outer_panel_i = fill(false, (nkyhat, nefwid, nfactor, inputsEP.NMODES))
    l_QL_ratio_i = fill(false, (nkyhat, nefwid, nfactor, inputsEP.NMODES))
    l_theta_sq_i = fill(false, (nkyhat, nefwid, nfactor, inputsEP.NMODES))
    l_theta_sq_i_out = fill(false, (nkyhat, nefwid, nfactor, inputsEP.NMODES))
    #this variable was not in the code, but used later per the commit
    l_QL_ratio_i_out = fill(false, (nkyhat, nefwid, nfactor, inputsEP.NMODES))
    l_QL_ratio_i_ = fill(false, (nkyhat, nefwid, nfactor, inputsEP.NMODES))

    f0 = zero(T)
    f1 = inputsEP.FACTOR_IN # 1 for first round
    w0 = inputsEP.WIDTH_MIN
    w1 = inputsEP.WIDTH_MAX
    kyhat_min = zero(T)
    kyhat_max = one(T)
    kyhat0 = kyhat_min
    kyhat1 = kyhat_max
    factor = fill(T(NaN), nfactor)
    efwid = fill(T(NaN), nefwid)
    kyhat = fill(T(NaN), nkyhat)

    ikyhat_write = floor(Int, nkyhat/2) # 2
    # println(ikyhat_write)
    iefwid_write = floor(Int, nefwid/2) # 5
    ifactor_write = nfactor # 10
    f_guess_mark = T(1.0E20)
    lkeep_ref = fill(false, (nkyhat, nefwid))
    # for scope purposes:
    inputTJLF = TJLFEP.InputTJLF{T}(inputsPR.NS, 12, true)
    imark_min = 0
    f_guess = fill(T(NaN), (nkyhat, nefwid))
    ikyhat_mark::Int64 = 0
    iefwid_mark::Int64 = 0
    fmark = T(1.0E20)
    scalefactor_buffer = l_write_out && printout ? String[] : nothing
    wavebuffer_all = []
    for k = 1:k_max
        # k doesn't use Threads!
        fill!(factor, T(NaN))
        fill!(efwid, T(NaN))
        fill!(kyhat, T(NaN))
        for i = 1:nfactor
            factor[i] = ((f1-f0)/nfactor)*i+f0
            # k = 1: FACTOR_IN = 1.0
            # [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
        end
        for i = 1:nefwid
            efwid[i] = ((w1-w0)/(nefwid-1))*(i-1)+w0
        end
        
        for i = 1:nkyhat
            if kyhat0 == 0.0
                kyhat[i]= ((kyhat1-kyhat0)/nkyhat)*i + kyhat0
            else 
                kyhat[i]= ((kyhat1-kyhat0)/(nkyhat-1))*(i-1)+kyhat0
            end
            #kyhat[i] = ((kyhat1-kyhat0)/nkyhat)*i+kyhat0
        end
        # eigen_cache is a sequential-seeding optimization incompatible with parallel execution;
        # each thread gets its own local cache initialised to nothing.
        wavebuffer_all = []
        Threads.@threads for i = 1:nkwf
            local_inputsEP = deepcopy(inputsEP)  # thread-local copy; avoids races on FACTOR_IN/KYHAT_IN/WIDTH_IN/LKEEP/etc.
            local_eigen_cache = nothing           # thread-local; no cross-iteration seeding in parallel

            l_wavefunction_out = 0

            # The following 3 statments define each combination of ikyhat, iefwid, and ifactor.
            ikyhat = Int(floor((i-1)/(nefwid*nfactor))+1) # 2
            iefwid = Int(floor(1.0*mod(i-1, nefwid*nfactor)/nfactor)+1) # 5
            ifactor = mod(i-1, nfactor)+1 # 10

            local_inputsEP.FACTOR_IN = factor[ifactor] # Just used in mapping
            local_inputsEP.KYHAT_IN = kyhat[ikyhat] # Set equal to ky_in
            local_inputsEP.WIDTH_IN = efwid[iefwid]
            debug_dump_kw_combo(local_inputsEP, i)

            str_sf = string(Char(mod(floor(Int, local_inputsEP.FACTOR_IN/100.0), 10) + UInt32('0'))) *
                     string(Char(mod(floor(Int, local_inputsEP.FACTOR_IN/10.0), 10) + UInt32('0')))  *
                     string(Char(mod(floor(Int, local_inputsEP.FACTOR_IN), 10) + UInt32('0'))) *
                     "." *
                     string(Char(mod(floor(Int, 10*local_inputsEP.FACTOR_IN), 10) + UInt32('0'))) *
                     string(Char(mod(floor(Int, 100*local_inputsEP.FACTOR_IN), 10) + UInt32('0'))) *
                     string(Char(mod(floor(Int, 1000*local_inputsEP.FACTOR_IN), 10) + UInt32('0')))

            str_wf_file = "out.wavefunction"*coalesce(local_inputsEP.SUFFIX, "")*"_sf"*str_sf

            if ((local_inputsEP.WRITE_WAVEFUNCTION == 1) &&
                (ikyhat == ikyhat_write) &&
                (iefwid == iefwid_write) &&
                (ifactor == ifactor_write) &&
                (k == k_max) && printout)
                l_wavefunction_out = 1
            end

            gamma_out, freq_out, inputTJLF, local_eigen_cache, wavefunction_buffer = TJLFEP_ky(local_inputsEP, inputsPR, str_wf_file, l_wavefunction_out; eigen_cache=local_eigen_cache, use_gpu=use_gpu)
            if wavefunction_buffer !== nothing
                push!(wavebuffer_all, (str_wf_file, wavefunction_buffer))
            end

            for n = 1:local_inputsEP.NMODES
                growthrate[ikyhat,iefwid,ifactor,n] = gamma_out[n]
                frequency[ikyhat,iefwid,ifactor,n]  = freq_out[n]
                lkeep_i[ikyhat,iefwid,ifactor,n] = local_inputsEP.LKEEP[n]
                ltearing_i[ikyhat,iefwid,ifactor,n] = local_inputsEP.LTEARING[n]
                l_th_pinch_i[ikyhat,iefwid,ifactor,n] = local_inputsEP.L_TH_PINCH[n]
                l_i_pinch_i[ikyhat,iefwid,ifactor,n] = local_inputsEP.L_I_PINCH[n]
                l_e_pinch_i[ikyhat,iefwid,ifactor,n] = local_inputsEP.L_E_PINCH[n]
                l_EP_pinch_i[ikyhat,iefwid,ifactor,n] = local_inputsEP.L_EP_PINCH[n]

                l_theta_sq_i[ikyhat,iefwid,ifactor,n] = local_inputsEP.L_THETA_SQ[n]
                #l_max_outer_panel_i[ikyhat,iefwid,ifactor,n] = local_inputsEP.L_MAX_OUTER_PANEL[n]
                l_QL_ratio_i[ikyhat,iefwid,ifactor,n] = local_inputsEP.L_QL_RATIO[n]
            end
            
        end # end Threads.@threads loop

        # This loop creates a 5x10 matrix full of 11. It then runs
        # through all dimensions of lkeep_i, which is a reference matrix
        # telling you where each 
        imark = fill(nfactor+1, (nkyhat, nefwid))
        #if (id == 0)
        #    println(imark, " imark at id: ", id)
        #end

        # The following loop sets the matrix imark finds the first mode for a specific ikyhat, iefwid, ifactor combination
        # that is set to be kept in lkeep_i. This mode is thus representing that specific ifactor value for the run
        # and imark is set to ifactor for that ikyhat and iefwid. The same is also done for the factor values as once the first non-zeros
        # growthrate is found along the factor spectrum, we return back to the next width run. Essentially, we are trying to find
        # all of the first non-zero growthrates at each kyhat, width combo; these are then allotted to the imark matrix which
        # describes whether a kyhat, width combo even has a non-zero growthrate, and if it does, where in the scan did it find that.

        for ikyhat = 1:nkyhat
            for iefwid = 1:nefwid
                for ifactor = 1:nfactor
                    for n = 1:inputsEP.NMODES # For a specified mode, 
                        if (lkeep_i[ikyhat, iefwid, ifactor, n])
                            imark[ikyhat, iefwid] = ifactor
                            # these could range from 1 to 10.
                            break
                        end
                    end # n
                    if (imark[ikyhat,iefwid] <= nfactor) break end
                end # ifactor
            end # iefwid
        end #ikyhat
        
        imark_min = nfactor + 1
        for ikyhat = 1:nkyhat
            for iefwid = 1:nefwid
                imark_min = min(imark[ikyhat, iefwid], imark_min)
            end
        end

        
        fmark = T(1.0E20)
        gmark = zero(T)
        f_guess_mark = T(1.0E20)
        gamma_mark_i_1 = fill(T(NaN), (nkyhat, nefwid))
        gamma_mark_i_2 = fill(T(NaN), (nkyhat, nefwid))
        f_mark_i = fill(T(NaN), (nkyhat, nefwid))
        
        # if there are any unstable modes:
        if (imark_min <= nfactor)

            # default the marked values as the last ones
            ikyhat_mark = nkyhat
            iefwid_mark = nefwid
            # Reset write indices to midpoint defaults each k (matches Fortran)
            ikyhat_write = floor(Int, nkyhat/2)
            iefwid_write = floor(Int, nefwid/2)

            # loop all combos of kyhat, width in imark:
            for ikyhat = 1:nkyhat
                for iefwid = 1:nefwid
                    lkeep_ref[ikyhat, iefwid] = false  # reset each k iteration (matches Fortran)
                    imark_ref = nfactor-1 # ? What is the point ?
                    #println("pass 1: ", imark_ref)
                    
                    # if imark is not the final one of the factor spectrum, set the reference to imark
                    # and set it as kept; otherwise, if it's the last one, set the reference to imark-1
                    # (9) and set as kept, and if it's a default value of 11, set reference to 10.
                    if (imark[ikyhat, iefwid] < nfactor) # < 10
                        imark_ref = imark[ikyhat, iefwid]
                        #println("route 1: ", imark_ref)
                        lkeep_ref[ikyhat, iefwid] = true
                    else # 10 or 11
                        #imark_ref = imark[ikyhat, iefwid] - 1
                        #println("route 2: ", imark_ref)
                        if (imark[ikyhat, iefwid] == nfactor) # == 10
                            lkeep_ref[ikyhat,iefwid] = true
                        end
                    end
                    # Now we set the f_g1 and f_g2 values, which allows us to guess the
                    # value of f_guess based on this equation. This equation is not clear yet to me.
                    f_g1 = factor[imark_ref]
                    f_g2 = factor[imark_ref+1]

                    # This is the portion that is NOT consistent with the Fortran and is
                    # causing problems due to default values:
                    gamma_g1 = zero(T) #-2.0
                    gamma_g2 = T(100.0) #-1.0
                    
                    # For all modes, set gamma_g1 to the maximum growthrate of kept modes at this point of ikyhat, iefwid, and imark_ref
                    # This also sets gamma_g2 to the maximum neighboring factor value of growthrate (imark_ref+1).
                    for n = 1:inputsEP.NMODES
                        if (lkeep_i[ikyhat, iefwid, imark_ref, n])
                            gamma_g1 = max((growthrate[ikyhat, iefwid, imark_ref, n]), gamma_g1)
                        end
                        # if (imark_ref+1 < nfactor+1)
                        if (lkeep_i[ikyhat, iefwid, imark_ref+1, n])
                            gamma_g2 = max((growthrate[ikyhat, iefwid, imark_ref+1, n]), gamma_g2)
                        end
                        # end
                    end
                    f_guess[ikyhat, iefwid] = f_g1 + (inputsEP.GAMMA_THRESH-gamma_g1)*(f_g2-f_g1)/(gamma_g2-gamma_g1)
                    gamma_mark_i_1[ikyhat, iefwid] = gamma_g1
                    gamma_mark_i_2[ikyhat, iefwid] = gamma_g2
                    if imark[ikyhat, iefwid] < nfactor 
                        f_mark_i[ikyhat, iefwid] = f_g1
                    else
                        f_mark_i[ikyhat,iefwid] = f_g2
                    end
                end # iefwid
            end # ikyhat

            # Now we scan across imark again, looking for modes that are:
            # kept; have a max growthrate along the modes that is less than 95% of the neighboring maximum; have a factor value that is 
            # less than or equal to the current maximum... (starts at 1e20)
            # then if these are satisfied, find modes that:
            # have a growthrate that is larger than the current maximum (starts at 0) or
            # have a factor that is smaller than the max
            ikymin =nkyhat
            ikymax = 1
            iefwmin = nefwid
            iefwmax = 1
            for ikyhat = 1:nkyhat
                for iefwid = 1:nefwid
                    if (lkeep_ref[ikyhat, iefwid] && f_mark_i[ikyhat, iefwid] <= fmark)
                        if f_mark_i[ikyhat, iefwid] < fmark
                            # If all of these are satisfied, set new marks away from default:
                            ikyhat_write = ikyhat
                            iefwid_write = iefwid
                            fmark = f_mark_i[ikyhat,iefwid]
                            ikymin = nkyhat
                            ikymax = 1
                            iefwmin = nefwid
                            iefwmax = 1
                            f_guess_mark = T(1.0E20)
                            ikyhat_mark = ikyhat
                            iefwid_mark = iefwid
                        end
                        ikymin = min(ikyhat,ikymin)
                        ikymax = max(ikyhat,ikymax)
                        iefwmin = min(iefwid,iefwmin)
                        iefwmax = max(iefwid,iefwmax)
                        f_guess_mark = min(f_guess[ikyhat,iefwid],f_guess_mark)

                            # then set new maximums (fmark and gmark)
                            # and set f_guess_mark which corresponds to a maximized growthrate and
                            # minimized factor.
 
                            

                            # This loops over all combos of kyhat and width to find a single guess for "f"
                        

                            #=if (id == 0)
                                println("Statements set! for id & ir and k: ", id, " ", inputsEP.IR, " ", k)
                                println("kywrite, wdwrite: ", ikyhat_mark, " ", iefwid_mark)
                                println("After: ", fmark, ", ", gmark, ", ", f_guess_mark)
                            end
                            sleep(0.25)
                            if (id == 1)
                                println("Statements set! for id & ir and k: ", id, " ", inputsEP.IR, " ", k)
                                println("kywrite, wdwrite: ", ikyhat_mark, " ", iefwid_mark)
                                println("After: ", fmark, ", ", gmark, ", ", f_guess_mark)
                            end=#
                        
                    else
                        #println("Statements NOT set! for id & ir: ", id, " ", inputsEP.IR)
                    end
                end # iefwid
            end # ikyhat
        end # if ending (mode found)

        # Next is the writing of the out.scalefactor files.

        g = fill(T(NaN), inputsEP.NMODES)
        f = fill(T(NaN), inputsEP.NMODES)
        keep_label = fill("", inputsEP.NMODES)
        # Set FREQ_AE_UPPER and GAMMA_THRESH on the main inputsEP so that display labels
        # (BT/F/?) and kHz header lines are computed correctly.  The threaded scan used
        # deepcopy'd local structs and never wrote these back to the shared struct.
        inputsEP.FREQ_AE_UPPER = -abs(inputsPR.omegaGAM[inputsEP.IR])
        if inputsEP.ROTATIONAL_SUPPRESSION_FLAG == 1
            r_over_a = inputsPR.RMIN[inputsEP.IR] / inputsPR.RMIN[end]
            inputsEP.GAMMA_THRESH_MAX = abs(inputsPR.gammap[inputsEP.IR]) * 2.0 * (min(1.0 - r_over_a, r_over_a) / inputsPR.RMAJ[inputsEP.IR])
            inputsEP.GAMMA_THRESH = 0.15 * abs(inputsPR.gammaE[inputsEP.IR] / inputsPR.SHEAR[inputsEP.IR])
            inputsEP.GAMMA_THRESH = min(inputsEP.GAMMA_THRESH, inputsEP.GAMMA_THRESH_MAX)
        else
            inputsEP.GAMMA_THRESH = 1.0e-7
            inputsEP.GAMMA_THRESH_MAX = 1.0e-7
        end
        if (l_write_out && printout)
            filename = "out.scalefactor"*coalesce(inputsEP.SUFFIX, "")
            # Header written only once (first k), matching Fortran's file-create-once behavior
            if isempty(scalefactor_buffer)
                push!(scalefactor_buffer, "factor,(gamma(n),freq(n),flag,n=1,nmodes_in)")
                push!(scalefactor_buffer, "flag key:  'K'  mode is kept")
                if (inputsEP.REJECT_TEARING_FLAG == 1) push!(scalefactor_buffer, "           'T' rejected for tearing parity") end
                if (inputsEP.REJECT_I_PINCH_FLAG == 1) push!(scalefactor_buffer, "           'Pi' rejected for ion pinch") end
                if (inputsEP.REJECT_E_PINCH_FLAG == 1) push!(scalefactor_buffer, "           'Pe' rejected for electron pinch") end
                if (inputsEP.REJECT_TH_PINCH_FLAG == 1) push!(scalefactor_buffer, "           'Pth' rejected for thermal pinch") end
                if (inputsEP.REJECT_EP_PINCH_FLAG == 1) push!(scalefactor_buffer, "           'PEP' rejected for EP pinch") end
                push!(scalefactor_buffer, "           'QLR' rejected for QL ratio EP/|chi_i| < $(inputsEP.QL_RATIO_THRESH)")
                push!(scalefactor_buffer, "           'TH2' rejected for <theta^2>  > $(inputsEP.THETA_SQ_THRESH)")
                push!(scalefactor_buffer, "           'F'   rejected for non-AE frequency > $(inputsEP.F_REAL[inputsEP.IR]*inputsEP.FREQ_AE_UPPER) kHz")
                push!(scalefactor_buffer, "           'BT'  mode growth rate is below threshold gamma_thresh = $(inputsEP.F_REAL[inputsEP.IR] * inputsEP.GAMMA_THRESH) kHz")
                push!(scalefactor_buffer, "omega_TAE = $(inputsEP.F_REAL[inputsEP.IR] * inputsPR.OMEGA_TAE[inputsEP.IR]) ;  omega_GAM = $(-inputsEP.F_REAL[inputsEP.IR] * inputsPR.omegaGAM[inputsEP.IR])")
                if inputsEP.REAL_FREQ != 0
                    push!(scalefactor_buffer, "Frequencies in real units, plasma frame [kHz] ; (c_s/a)/(2*pi) = $(inputsEP.F_REAL[inputsEP.IR]) kHz")
                end
            end

            ky_write = kyhat[ikyhat_write]*inputTJLF.ZS[inputsEP.IS_EP+1]/sqrt(inputTJLF.MASS[inputsEP.IS_EP+1]*inputTJLF.TAUS[inputsEP.IS_EP+1])

            push!(scalefactor_buffer, _scalefactor_section_header(kyhat[ikyhat_write], ky_write, efwid[iefwid_write], fmark))

            for ifactor = 1:nfactor
                ikyhat = ikyhat_write
                iefwid = iefwid_write
                for n = 1:inputsEP.NMODES
                    g[n] = inputsEP.F_REAL[inputsEP.IR]*growthrate[ikyhat, iefwid, ifactor, n]
                    f[n] = inputsEP.F_REAL[inputsEP.IR]*frequency[ikyhat, iefwid, ifactor, n]
                    if (lkeep_i[ikyhat,iefwid,ifactor,n])
                        keep_label[n] = " K  "
                    elseif growthrate[ikyhat, iefwid, ifactor, n] < inputsEP.GAMMA_THRESH 
                        keep_label[n] = " BT  "
                    elseif ((ltearing_i[ikyhat,iefwid,ifactor,n]) && (inputsEP.REJECT_TEARING_FLAG == 1))
                        keep_label[n] = " T  "
                    elseif ((l_i_pinch_i[ikyhat,iefwid,ifactor,n]) && (inputsEP.REJECT_I_PINCH_FLAG == 1))
                        keep_label[n] = " Pi "
                    elseif ((l_e_pinch_i[ikyhat,iefwid,ifactor,n]) && (inputsEP.REJECT_E_PINCH_FLAG == 1))
                        keep_label[n] = " Pe "
                    elseif ((l_th_pinch_i[ikyhat,iefwid,ifactor,n]) && (inputsEP.REJECT_TH_PINCH_FLAG == 1))
                        keep_label[n] = " Pth"
                    elseif ((l_EP_pinch_i[ikyhat,iefwid,ifactor,n]) && (inputsEP.REJECT_EP_PINCH_FLAG == 1))
                        keep_label[n] = " PEP"
                    elseif (l_QL_ratio_i[ikyhat,iefwid,ifactor,n])
                        keep_label[n] = " QLR"
                    elseif l_theta_sq_i[ikyhat, iefwid, ifactor, n]
                        keep_label[n] = " TH2"
                    elseif (frequency[ikyhat,iefwid,ifactor,n] > inputsEP.FREQ_AE_UPPER)
                        keep_label[n] = " F  "
                    else
                        keep_label[n] = " ?  "
                    end
                end
                push!(scalefactor_buffer, _scalefactor_factor_line(factor[ifactor], g, f, keep_label))
            end # ifactor
        end 

        # After printing the info just calculated, if fmark is already significantly small, adjust 
        # ranges within the range, otherwise, move to a new order of 10 in the factor.
        # if (fmark < 1.0e10)
        #     # If in the first scan round, set the new max factor to the minimized factor, otherwise, do
        #     # some more specific things to hone in each round (including in kyhat and width)
        if fmark < 1.0E10  # accepted mode with all constraints
            f1 = 2.0 * fmark
            f0 = zero(T)
            delw = (w1 - w0) / (nefwid - 1)
        
            w1 = efwid[iefwmax] + delw
            w0 = efwid[iefwmin] - delw
            if w1 > inputsEP.WIDTH_MAX
                w1 =inputsEP.WIDTH_MAX
            end
            if w0 < inputsEP.WIDTH_MIN
                w0 = inputsEP.WIDTH_MIN
            end
        
            delky = (kyhat1 - kyhat0) / nkyhat
            kyhat1 = kyhat[ikymax] + delky
            kyhat0 = kyhat[ikymin] - delky
            if kyhat1 > kyhat_max
                kyhat1 = kyhat_max
            end
            if kyhat0 < kyhat_min
                kyhat0 = kyhat_min
            end
        
        else
            f0 = f1
            f1 = 10.0*f1
            #println("fmark < 1e10 : Pass ", k)
        end

    end # end of k
    #println(imark_min)
    
    # After the fourth round, if no unstable modes have been found, default the scalefactor which will be used In
    # the calculation in the driver to 10k, the width to its default (I think around 100) and same for kyhat.
    # Otherwise, set the factor to the best guess at ikyhat_mark, iefwid_mark and their corresponding width and kyhat
    # to those marks.
    if (imark_min > nfactor)
        # If, over the scan of k, there's no unstable modes, set each to lowest.
        inputsEP.FACTOR_IN = 10000
        inputsEP.WIDTH_IN = efwid[1]
        inputsEP.KYMARK = kyhat[1]
    
    else
        # If, over the scan of k, there's an unstable mode, set each to each marked point.
        inputsEP.FACTOR_IN = f_guess_mark
        inputsEP.WIDTH_IN = efwid[iefwid_mark]
        inputsEP.KYMARK = kyhat[ikyhat_mark]
    end

    # Return to the driver these values and the final growthrate values of the last scan.
    return growthrate, inputsEP, inputsPR, scalefactor_buffer, wavebuffer_all

    # This function will be done for however many radii you are testing. These values do not interact in the driver. 
end
