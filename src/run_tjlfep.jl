"""
runTHD(tglfepfilepath::String, mtglffilepath::String, exprofilepath::String)

inputs: tglfepfilepath, mtglffilepath, exprofilepath

inputs are used in the threads version of the TJLFEP code for a single run
"""
import TJLF


function runTHD(inputTJLFEP::Vector{InputTJLF}, Options::Options, profile::profile, printout::Bool=false)


    homedir = pwd()


    if (Options.IS_EP == 1)
        Ni = profile.DENS_1
        Ti = profile.TEMP_1
        DLNNDR = profile.DLNNDR_1
        DLNTDR = profile.DLNTDR_1

    elseif (Options.IS_EP == 2)
        Ni = profile.DENS_2
        Ti = profile.TEMP_2
        DLNNDR = profile.DLNNDR_2
        DLNTDR = profile.DLNTDR_2
    elseif (Options.IS_EP == 3)
        Ni = profile.DENS_3
        Ti = profile.TEMP_3
        DLNNDR = profile.DLNNDR_3
        DLNTDR = profile.DLNTDR_3
    elseif (Options.IS_EP == 4)
        Ni = profile.DENS_4
        Ti = profile.TEMP_4
        DLNNDR = profile.DLNNDR_4
        DLNTDR = profile.DLNTDR_4
    else
        println("InputTJLFoptions.IS_EP not within range. Check input.TGLFEP input")
        return 1
    end
    cs= profile.CS
    rmin_ex = profile.RMIN
    gammaE =profile.gammaE
    gammap =profile.gammap
    omegaGAM = profile.omegaGAM
    ir_exp::Vector{Int64} = []
    #= if ir_exp is needed can be coded here=#
    dpdr_EP = fill(NaN, profile.NR)
    if (Options.INPUT_PROFILE_METHOD == 2)
        for i in eachindex(dpdr_EP)
            dpdr_EP[i] = ni[i]*Ti[i]*(DLNNDR[i]+DLNTDR[i])# This has some small changes from old main
        end 
        dpdr_EP_abs = abs.(dpdr_EP)
        dpdr_EP_max = maximum(dpdr_EP_abs)
        dpdr_EP_max_loc = argmax(dpdr_EP_abs)
        n_at_max = ni[dpdr_EP_max_loc]
        if (Options.PROCESS_IN != 5)
            for ir = 1:Options.SCAN_N
                Options.FACTOR = Options.FACTOR*dpdr_EP_max/dpdr_EP_abs[ir_exp[ir]]   
            end
        end
        Options.FACTOR_MAX_PROFILE .= Options.FACTOR
    end

    Options.F_REAL .= 1.0
    if (Options.REAL_FREQ == 1) 
        Options_REAL .= (cs[:]/(rmin_ex[profile.NR]))/(2*pi*1.0e3)
    end

    if (ky_model == 0)
        Options.NTOROIDAL = 4
    else
        Options.NTOROIDAL = 3
    end
        #in actor now
    # if (process_in == 4 || process_in == 5)
    #     inputTJLFEP.NN = nn
    #end
    #
    if (!Options.FACTOR_IN_PROFILE)
        Options.FACTOR = fill(Options.FACTOR_IN, input_tjlf.Options.SCAN_N)
    end
    Options.FACTOR_MAX_PROFILE = Options.FACTOR


    if (Options.INPUT_PROFILE_METHOD == 2)
    # Allotting Ir_exp not from inputTJLFEP.profileGLF.
    Options.IR_EXP = fill(0, Options.SCAN_N)
        for i = 1:Options.SCAN_N
            if (Options.SCAN_N != 1)
                jr_exp = profile.IRS + floor((i-1)*(profile.NR-profile.IRS)/(Options.SCAN_N-1))
            else
                jr_exp = profile.IRS
            end
            Options.IR_EXP[i] = jr_exp
        end
    end


    # Things that need to be set for the function of the driver:
    # str_r & SUFFIX
    # FACTOR_IN
    # color & id (sort of)

    # deepcopy is required so as to avoid overwriting of data:
    n_ir = Options.SCAN_N
    Ts = fill(Options, n_ir)
    Ts[1] = deepcopy(Options)
    for i in 2:n_ir
        Ts[i] = deepcopy(Ts[i-1])
    end
    arrTGLFEP = Ts
    arrPROF = fill(profile, n_ir)
    arrgrowth = fill(fill(NaN,(5, 10, 10, Options.NMODES)), n_ir)

    for i in 1:n_ir
        println(inputTJLFEP[i])
        #try
        arrTGLFEP[i].IR = arrTGLFEP[i].IR_EXP[i]
        ir = arrTGLFEP[i].IR
        str_r = lpad(string(ir), 3, '0')
        arrTGLFEP[i].SUFFIX = "_r"*str_r
        arrTGLFEP[i].FACTOR_IN = arrTGLFEP[i].FACTOR[i]
        input1 = arrTGLFEP[i]
        input2 = arrPROF[i]
        arrgrowth[i], arrTGLFEP[i], arrPROF[i] = TJLFEP.mainsub(input1, input2, inputTJLFEP[i], printout)
        #catch
        #end
    end

    # Right now, width_in is being set directly from first entry from the TGLFEP input
    # rather than some range of widths as TGLFEP does. This can be easily implemented but 
    # most inputs won't need that for now.

    #println(Options.WIDTH_IN, " for id & color ", id_world, " ", color)

    kymark_out::Vector{Float64} = fill(NaN, Options.SCAN_N)
    width::Vector{Float64} = fill(NaN, Options.SCAN_N)
    #local buf_width::Vector{Float64} = fill(NaN, 1)
    #local buf_kymark::Vector{Float64} = fill(NaN, 1)

    if (Options.WIDTH_IN_FLAG)
        # Non-MPI:
        # There are only "3" processes in Threads -- 
        for i = 1:n_ir
            width[i] = arrTGLFEP[i].WIDTH_IN
            kymark_out[i] = arrTGLFEP[i].KYMARK
            
        end

    end

    # Print out basic information about the run (that is common to all radii):
    Options = arrTGLFEP[1]
    if (false)printout
        io2 = open("out.TGLFEP", "w")
        println(io2, "process_in = ", Options.PROCESS_IN)

        if (Options.PROCESS_IN <= 1) println(io2, "mode_in = ", Options.MODE_IN) end
        if ((Options.PROCESS_IN == 4) || (Options.PROCESS_IN == 5)) println(io2, "threshold_flag = ", Options.THRESHOLD_FLAG) end

        println(io2, "ky_mode = ", Options.KY_MODEL)
        println(io2, "--------------------------------------------------------------")
        println(io2, "scan_n = ", Options.SCAN_N)
        println(io2, "irs = ", Options.IRS)
        println(io2, "n_basis = ", Options.N_BASIS)
        println(io2, "scan_method = ", Options.SCAN_METHOD)

        if (Options.WIDTH_IN_FLAG)
            println(io2, "ir,  width")
            for i = 1:Options.SCAN_N
                println(io2, Options.IRS+i-1, " ", width[i])
            end
        else
            println(io2, "ir,  width,  kymark")
            for i = 1:Options.SCAN_N
                println(io2, Options.IRS+i-1, " ", width[i], " ", kymark_out[i])
            end
        end

        println(io2, "--------------------------------------------------------------")
        println(io2, "factor_in_profile = ", Options.FACTOR_IN_PROFILE)
        if (Options.FACTOR_IN_PROFILE)
            for i = 1:Options.SCAN_N
                println(io2, Options.FACTOR[i])
            end
        else
            println(io2, Options.FACTOR[1])
        end

        println(io2, "width_in_flag = ", Options.WIDTH_IN_FLAG)
        if (!Options.WIDTH_IN_FLAG) println(io2, "width_min = ", Options.WIDTH_MIN, " width_max = ", Options.WIDTH_MAX) end
        close(io2)
    end

    # Lastly is the printing of the density threshold:
    #local buf_factor::Vector{Float64} = fill(NaN, 1)

    # Now continue on to radii-dependent part:
    if ((Options.PROCESS_IN == 4) || (Options.PROCESS_IN == 5))
        # Threads Process won't happen until later:


        # Original MPI:
        #=if (id_world < Options.SCAN_N && id_world != 0)
            if (Options.THRESHOLD_FLAG == 0)
                #buf_factor = [Options.FACTOR_IN]
                MPI.Send([Options.FACTOR_IN], 0, id_world, MPI.COMM_WORLD)
                #Options.FACTOR_IN = buf_factor[1]
            end
            #MPI.Send(Options) only process-In 4 here.
        end # id == 0
        MPI.Barrier(TJLFEP_COMM_IN)=#
        
        if (printout)
        io3 = open("out.TGLFEP", "a")
        println(io3, "**************************************************************")
        println(io3, "************** The critical EP density gradient **************")
        println(io3, "**************************************************************")
        end

        SFmin = fill(0.0, Options.SCAN_N)
        SFmin_out = fill(0.0, profile.NR)
        dndr_crit = fill(NaN, Options.SCAN_N)
        dndr_crit_out = fill(NaN, profile.NR)
        dpdr_crit = fill(NaN, Options.SCAN_N)
        dpdr_crit_out = fill(NaN, profile.NR)

        if (Options.THRESHOLD_FLAG == 0)
            # Threads Correct:
            # Set the returned guesses for "f_guess" as SFmin[i]:
            for i = 1:n_ir
                SFmin[i] = arrTGLFEP[i].FACTOR_IN
            end
            # MPI Original:
            #SFmin[1] = Options.FACTOR_IN
            #println(io3, Options.FACTOR_IN, " factor_in before")
            #println("Before MPI.Recv! for factor_in.")
            #=for i = 1:Options.SCAN_N-1
                buf_factor = [NaN]
                MPI.Recv!(buf_factor, i, i, MPI.COMM_WORLD)
                SFmin[i+1] = buf_factor[1]
                #println(io3, SFmin[i_1], " factor_in after and ", Options.FACTOR_IN, " buf_factor after")
                #println(io3, SFmin) # before buf_factor
                #SFmin[i+1] = buf_factor[1]
                #println(io3, SFmin) # after buf_factor, before comp.out
            
            end=#
            if (printout)
                println("After MPI.Recv! for factor_in")
                println(io3, "--------------------------------------------------------------")
                println(io3, "SFmin")
            end
        
        # Next is TGLFEP_complete_output(SFmin, SFmin_out, ir_min, ir_max, l_accept_profile)
        # This function's goal is to determine whether 

            SFmin, SFmin_out, ir_min, ir_max, l_accept_profile = tjlfep_complete_output(SFmin, Options, profile)
        
            if (printout)
                println(io3, SFmin, " SFmin after buf and coutput") # after comp.out
            end
            #println(io3, Options.FACTOR_MAX_PROFILE)

            # We've received the altered profile (interpolated and accepted or not).
            # If the minimum radius is not the first one...
            if ((ir_min-Options.IRS+1) > 1)
                # Originally, this had no concern for accessing out-of-range values;
                # it does now:

                # Set any scans before this point's factor values to the factor_max_profile values respectively (?? Why so physically)
                # If you're running a normal amount of scans, this will pretty much never be done, right?
                if (ir_min-Options.IRS > Options.SCAN_N)
                    SFmin[1:Options.SCAN_N] = Options.FACTOR_MAX_PROFILE[1:Options.SCAN_N]
                else # original alone:
                    SFmin[1:ir_min-Options.IRS] = Options.FACTOR_MAX_PROFILE[1:ir_min-Options.IRS]
                end

                # If the starting radius is greater than 1, set the values before it in the interpolated profile to
                # the first value of the max_profile
                if (Options.IRS > 1) SFmin_out[1:Options.IRS-1] .= Options.FACTOR_MAX_PROFILE[1] end

                # This one doesn't look right...
                # This says that in the interpolated profile, you should set any values from the initial to the first point (minus 1)
                # to the same factor_max_profile. But this factor_max_profile is not an interpolation...
                # Is this being done as a default value? FACTOR_MAX_PROFILE if just scan_n of the same value for the case I've been testing
                # hence it's return of 1.0 (*) when rejected. The problem is that this doesn't make much sense for SFmin_out...

                # The problem is that this ignores defaults again. If factor_max_profile is accessed outside of scan_n, they should be set to 0...
                if (ir_min-Options.IRS > Options.SCAN_N)
                    # The +1 on SFmin_out exists because the original is keeping a spacing of 1 between the two...
                    SFmin_out[Options.IRS:Options.SCAN_N+1] = Options.FACTOR_MAX_PROFILE[1:Options.SCAN_N]

                    # This is a test input:
                    # SFmin_out[Options.SCAN_N+1:ir_min-Options.IRS] = 0.0
                else
                    SFmin_out[Options.IRS:ir_min-1] = Options.FACTOR_MAX_PROFILE[1:ir_min-Options.IRS]
                end
            end

            # Perform a similar maneuver for above the maximum:
            if ((ir_max-Options.IRS+1) < Options.SCAN_N)
                SFmin[ir_max-Options.IRS+2:Options.SCAN_N] = Options.FACTOR_MAX_PROFILE[ir_max-Options.IRS+2:Options.SCAN_N]
                if (Options.IRS+Options.SCAN_N-1 < profile.NR) SFmin_out[Options.IRS+Options.SCAN_N:profile.NR] .= Options.FACTOR_MAX_PROFILE[Options.SCAN_N] end
                SFmin_out[ir_max+1:Options.IRS+Options.SCAN_N-1] = Options.FACTOR_MAX_PROFILE[ir_max-Options.IRS+2:Options.SCAN_N]
            end

            if (printout)
                println(io3, SFmin, " SFmin after Max assign")
            end

            if (printout)
                for i = 1:inputTGLFEP.SCAN_N
                    if (l_accept_profile[i])
                        println(io3, SFmin[i])
                    else
                        println(io3, SFmin_out[i+inputTGLFEP.IRS-1], "   (*)")
                    end
                    println(io3, "--------------------------------------------------------------")
                    println(io3, l_accept_profile)
                    println(io3, ir_min, " : ", ir_max)
                end
            end

            # Calculate the density critical gradient at each of the scanned radii.
            if (inputTGLFEP.INPUT_PROFILE_METHOD == 2)
                dndr_crit .= 10000.0
                for i = 1:inputTGLFEP.SCAN_N
                    # If SFmin[i] is not the default non-rejected value, multiply the scalefactor by the density and the density gradient at that point 
                    # for the energetic ion. 
                    # If SFmin[i] is the default or >= 9k, check if it is one of the factor_max_profile ones, and if so, calculate it with that.
                    # otherwise, leave it at 10k.
                    if (SFmin[i] < 9000.0)
                        dndr_crit[i] = SFmin[i]*ni[Int(inputTGLFEP.IR_EXP[i])]*dlnnidr[Int(inputTGLFEP.IR_EXP[i])]
                    elseif ((i < ir_min-inputTGLFEP.IRS+1) || (i > ir_max-inputTGLFEP.IRS+1))
                        dndr_crit[i] = inputTGLFEP.FACTOR_MAX_PROFILE[i]*ni[Int(inputTGLFEP.IR_EXP[i])]*dlnnidr[Int(inputTGLFEP.IR_EXP[i])]
                    end
                end
                # Interpolate and accept are reject needed values of this profile:
                dndr_crit, dndr_crit_out, ir_dum_1, ir_dum_2, l_accept_profile = tjlfep_complete_output(dndr_crit, inputTGLFEP, profile)
                
                if (printout)
                    io4 = open("alpha_dndr_crit.input", "w")
                    println(io4, "Density critical gradient (10^19/m^4)")
                    println(io4, dndr_crit_out)
                    close(io4)
                end
            end

            if (inputTGLFEP.INPUT_PROFILE_METHOD == 2)
                dpdr_crit .= 10000.0
                dpdr_EP[:] .= ni[:].*Ti[:].*(dlnnidr[:].+dlntidr[:]).*0.16022
                for i = 1:inputTGLFEP.SCAN_N
                    if (SFmin[i] < 9000.0)
                        if ((inputTGLFEP.PROCESS_IN == 4) || (inputTGLFEP.PROCESS_IN == 5))
                            case = inputTGLFEP.SCAN_METHOD
                            if (case == 1)
                                dpdr_scale = SFmin[i]
                            elseif (case == 2)
                                dpdr_scale = ((SFmin[i]*dlnnidr[inputTGLFEP.IR_EXP[i]]+dlntidr[inputTGLFEP.IR_EXP[i]]) /
                                (dlnnidr[inputTGLFEP.IR_EXP[i]]+dlntidr[inputTGLFEP.IR_EXP[i]]))
                            end
                            dpdr_crit[i] = dpdr_scale*dpdr_EP[inputTGLFEP.IR_EXP[i]]
                        end # 4 || 5
                    end # < 9000
                end # over scan_n
                dpdr_crit, dpdr_crit_out, ir_dum_1, ir_dum_2, l_accept_profile = tjlfep_complete_output(dpdr_crit, inputTGLFEP, profile)
                
                if (printout)
                    io5 = open("alpha_dpdr_crit.input", "w")
                    println(io5, "Pressure critical gradient (10 kPa/m)")
                    println(io5, dpdr_crit_out)
                    close(io5)
                end
            end # end prof. method 2

            if (printout)
                println(io3, "--------------------------------------------------------------")
                println(io3, "The EP density threshold n_EP/n_e (%) for gamma_AE = 0")
                for i = 1:inputTGLFEP.SCAN_N
                    println(io3, SFmin[i]*profile.AS[inputTGLFEP.IRS+i-1, profile.IS]*100.0) #percent
                end
            end
            
            if (printout)
                println(io3, "--------------------------------------------------------------")
                println(io3, "The EP beta crit (%) = beta_e*(n_EP_th/n_e)*(T_EP/T_e)")
                for i = 1:inputTGLFEP.SCAN_N
                    if (profile.GEOMETRY_FLAG == 0)
                        println(io3, SFmin[i]*profile.BETAE[inputTGLFEP.IRS+i-1]*100.0*profile.AS[inputTGLFEP.IRS+i-1, profile.IS]*profile.TAUS[inputTGLFEP.IRS+i-1, profile.IS]) #percent
                    else
                        println(io3, SFmin[i]*profile.BETAE[inputTGLFEP.IRS+i-1]*100.0*profile.AS[inputTGLFEP.IRS+i-1, profile.IS]*profile.TAUS[inputTGLFEP.IRS+i-1, profile.IS]*profile.KAPPA[inputTGLFEP.IRS+i-1]^2) #percent
                    end
                end
            end

            # there is a process_in == 4 addition I won't be doing quite yet.
        else # ThreshFlag != 0
            # Skipping for now as I want to test just threshold flag == 0 first
        end # ThreshFlag
    end # process 4 || 5

    if (printout)
        close(io3)
    end

    #MPI.Finalize()
    # The rest of the fortran is deallocations and a process_in == 6 route

    # Checking out some info:
    # @profview TJLFEP.mainsub(inputTGLFEP, profile, TJLFEP_COMM_IN)
    # using BenchmarkTools
    # @btime TJLFEP.mainsub(inputTGLFEP, profile, TJLFEP_COMM_IN)
    # The following list is for unpopulated fields and things I want to look into today:
    # inputTJLFEP.USE_AVE_ION_GRID
    # inputTJLFEP.WIDTH_SPECTRUM
    # inputTJLFEP.FIND_EIGEN
    # inputTJLFEPRLNP_CUTOFF
    # inputTJLFEP.BETA_LOC
    # inputTJLFEP.DAMP_PSI
    # inputTJLFEP.DAMP_SIG
    # inputTJLFEP.WDIA_TRAPPED

    #ADD whatever terms you want to output into a dictionary for analysis

    TJLFEPoutput = Dict{String, Any}(
        "width" => width,
        "kymark_out" => kymark_out,
        "SFmin" => SFmin,
        "dpdr_crit_out" => dpdr_crit_out,
        "dndr_crit_out" => dndr_crit_out
    )
    graph = true
    if (graph)
        x_value1 = range(1, stop=profile.NR, length = dict_out["dpdr_crit_out"])
        x_value2 = range(1, stop=profile.NR, length = dict_out["dndr_crit_out"])

        # First plot for the "Fortran" dataset
        plot(x_value1, dict_out["dpdr_crit_out"], xlabel="Radius", ylabel="dp/dr", 
        title="Critical Pressure Gradient", legend=:topright, linewidth=2, linestyle=:solid)

        # Second plot for the "Julia" dataset
        plot(x_value2, dict_out["dndr_crit_out"], xlabel="Radius", ylabel="dn/dr", 
        title="Critical Density Gradient", legend=:topright, linewidth=2, linestyle=:solid)
    end
    return TJLFEPoutput
end

#=
"""
runMPIs(tglfepfilepaths::Vector{String}, mtglffilepaths::Vector{String}, exprofilepaths::Vector{String})

Options: tglfepfilepaths, mtglffilepaths, exprofilepaths

These Options are vectors of individual filepaths for the MPI version of the TJLFEP code. It will run the code for each input.
"""
function runMPIs(tglfepfilepaths::Vector{String}, mtglffilepaths::Vector{String}, exprofilepaths::Vector{String})
for i = 1:length(tglfepfilepaths)
    runMPI(tglfepfilepaths[i], mtglffilepaths[i], exprofilepaths[i])
end
end
=#

# runTHDs(tglfepfilepaths::Vector{String}, mtglffilepaths::Vector{String}, exprofilepaths::Vector{String})

# Options: tglfepfilepaths, mtglffilepaths, exprofilepaths

# These inputs are vectors of individual filepaths for the threads version of the TJLFEP code. It will run the code for each input.
# It is probably advised to only use this function on runs that use the same # of scans, 
# """
##Have not updated For Fuse
# function runTHDs(tglfepfilepaths::Vector{String}, mtglffilepaths::Vector{String}, exprofilepaths::Vector{String})
# nruns = length(tglfepfilepaths)

# # Set return vectors:
# widths = fill([], nruns)
# kymark_outs = fill([], nruns)
# SFmins = fill([], nruns)
# dpdr_crit_outs = fill([], nruns)
# dndr_crit_outs = fill([], nruns)

# for i = 1:nruns
#     width, kymark_out, SFmin, dpdr_crit_out, dndr_crit_out = runTHD(tglfepfilepaths[i], mtglffilepaths[i], exprofilepaths[i])

#     widths[i] = width
#     kymark_outs[i] = kymark_out
#     SFmins[i] = SFmin
#     dpdr_crit_outs[i] = dpdr_crit_out
#     dndr_crit_outs[i] = dndr_crit_out
# end
# return widths, kymark_outs, SFmins, dpdr_crit_outs, dndr_crit_outs
# end



"""
checkInput(inputTJLF::InputTJLF)

description:
check that the InputTJLF struct is properly populated
"""
function checkInput(inputTJLFEP::InputTJLFEP)
field_names = fieldnames(InputTJLFEP)
for field_name in field_names
    field_value = getfield(inputTJLFEP, field_name)
    if typeof(field_value)<:Missing
        @assert !ismissing(field_value) "Did not properly populate inputTJLFEP for $field_name = $field_value"
    end
    if typeof(field_value)<:Real
        @assert !isnan(field_value) "Did not properly populate inputTJLFEP for $field_name = $field_value"
    end
    if typeof(field_value)<:Vector && field_name!=:KY_SPECTRUM && field_name!=:EIGEN_SPECTRUM
        for val in field_value
            @assert !isnan(val) "Did not properly populate inputTJLFEP for array $field_name = $val"
        end
    end
end
if !inputTJLFEP.FIND_EIGEN
    @assert !inputTJLFEP.FIND_WIDTH "If FIND_EIGEN false, FIND_WIDTH should also be false"
end
end

function checkInput(inputTJLFEPVector::Vector{InputTJLFEP})
for inputTJLFEP in inputTJLFEPVector
    field_names = fieldnames(inputTJLFEP)
    for field_name in field_names
        field_value = getfield(inputTJLFEP, field_name)
        if typeof(field_value)<:Missing
            @assert !ismissing(field_value) "Did not properly populate inputTJLFEP for $field_name = $field_value"
        end
        if typeof(field_value)<:Real
            @assert !isnan(field_value) "Did not properly populate inputTJLFEP for $field_name = $field_value"
        end
        if typeof(field_value)<:Vector && field_name!=:KY_SPECTRUM && field_name!=:EIGEN_SPECTRUM
            for val in field_value
                @assert !isnan(val) "Did not properly populate inputTJLFEP for array $field_name = $val"
            end
        end
    end
    if !inputTJLFEP.FIND_EIGEN
        @assert !inputTJLFEP.FIND_WIDTH "If FIND_EIGEN false, FIND_WIDTH should also be false"
    end
end
end
