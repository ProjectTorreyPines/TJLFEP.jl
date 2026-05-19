using Distributed

"""Unpack `mainsub` return for PROCESS_IN=5: `((growth, ep, mt), buffers)`."""
function _unpack_mainsub!(ret)
    (growth, ep, mt), _buffers = ret
    return growth, ep, mt
end

"""Run mainsub for one scan index (radius). Used by threads and distributed loops."""
function _runTHD_radius!(
    i::Int,
    arrTGLFEP,
    arrMTGLF,
    arrgrowth,
    printout::Bool,
    use_gpu::Bool;
    stdout_lock::Union{ReentrantLock,Nothing} = nothing,
)
    arrTGLFEP[i].IR = arrTGLFEP[i].IR_EXP[i]
    ir = arrTGLFEP[i].IR
    arrTGLFEP[i].SUFFIX = "_r" * lpad(string(ir), 3, '0')
    arrTGLFEP[i].FACTOR_IN = arrTGLFEP[i].FACTOR[i]
    if stdout_lock !== nothing
        lock(stdout_lock) do
            println("=============================================================")
            println("pre mainsub")
            println("i is ", i, " ir is ", ir)
            println("=============================================================")
        end
    end
    growth, ep, mt = _unpack_mainsub!(
        TJLFEP.mainsub(arrTGLFEP[i], arrMTGLF[i], printout; use_gpu=use_gpu))
    arrgrowth[i] = growth
    arrTGLFEP[i] = ep
    arrMTGLF[i] = mt
    return nothing
end

function _resolve_runTHD_parallel(parallel::Symbol)
    parallel == :auto && return nworkers() > 1 ? :distributed : :threads
    return parallel
end

function runTHD(tglfepfilepath::String, mtglffilepath::String, exprofilepath::String;
                printout::Bool=false, use_gpu::Bool=false, parallel::Symbol=:auto)

    # Auto-detect device via TJLF.pick_device(:auto); shadows the use_gpu parameter.
    # Thread safety: Threads.@threads runs each iteration in a separate Julia task.
    # CUDA.jl v5 assigns per-task streams, so concurrent GPU calls are stream-isolated.
    # use_gpu = TJLF.pick_device(:auto) === :gpu
    # processor = use_gpu ? "GPU" : "CPU"
    # println("TJLFEP runTHD: using $processor")

    # Default values for EXPRO:
    ni = TJLFEP.exproConst.ni
    Ti = TJLFEP.exproConst.Ti
    dlnnidr = TJLFEP.exproConst.dlnnidr
    dlntidr = TJLFEP.exproConst.dlntidr
    cs = TJLFEP.exproConst.cs
    rmin_ex = TJLFEP.exproConst.rmin_ex
    omegaGAM = TJLFEP.exproConst.omegaGAM
    gammaE = TJLFEP.exproConst.gammaE
    gammap = TJLFEP.exproConst.gammap
    # These should be set from the working directory, but these test cases are good for now:

    homedir = pwd()

    iEPexist::Bool = false
    iMPexist::Bool = false
    iEXPexist::Bool = false



    iEPexist = isfile(tglfepfilepath)
    iMPexist = isfile(mtglffilepath)
    iEXPexist = isfile(exprofilepath)


    @assert iEPexist != false "Requested TGLFEP input file path does not exist"
    @assert iMPexist != false "Requested MTGLF input file path does not exist"
    @assert iEXPexist != false "Requested EXPRO input file path does not exist"

    inputEPfile = tglfepfilepath
    inputMPfile = mtglffilepath
    inputEXPfile = exprofilepath

    #inputEPfile = "/Users/benagnew/TJLF.jl/outputs/tglfep_tests/input.TGLFEP"
    #inputMPfile = "/Users/benagnew/TJLF.jl/outputs/tglfep_tests/input.MTGLF"

    # Set up profile struct:
    prof = TJLFEP.readMTGLF(inputMPfile)
    profile = prof[1]
    ir_exp = prof[2]

    # Set up TGLFEP struct:
    Options = TJLFEP.readTGLFEP(inputEPfile, ir_exp)

    # IS_EP in input.TGLFEP is the GACODE ion index (EP species); EXPRO file has electron as species 1.
    gacode_dump = get(ENV, "GACODE_DUMP", nothing)
    ni, Ti, dlnnidr, dlntidr, cs, rmin_ex, gammaE, gammap, omegaGAM = TJLFEP.read_expro_for_alpha(
        inputEXPfile, profile, Options.IS_EP; gacode_file=gacode_dump,
    )

    profile.gammaE = gammaE
    profile.gammap = gammap
    profile.omegaGAM = omegaGAM

    println("options ir_exp is ", Options.IR_EXP)

    # If IR_EXP was not saved in the MTGLF file (e.g. legacy test cases), fall back to linear spacing
    if isempty(Options.IR_EXP)
        Options.IR_EXP = fill(0, Options.SCAN_N)
        for i = 1:Options.SCAN_N
            if (Options.SCAN_N != 1)
                jr_exp = profile.IRS + floor((i-1)*(profile.NR-profile.IRS)/(Options.SCAN_N-1))
            else
                jr_exp = profile.IRS
            end
            Options.IR_EXP[i] = jr_exp
        end
        println("IR_EXP not found in file, using linear spacing: ", Options.IR_EXP)
    end

    dpdr_EP = fill(NaN, profile.NR)
    if (Options.INPUT_PROFILE_METHOD == 2)
        for i in eachindex(dpdr_EP)
            dpdr_EP[i] = ni[i]*Ti[i]*(dlnnidr[i]+dlntidr[i])# This has some small changes from old main
        end
        #println(Options.FACTOR)
        dpdr_EP_abs = abs.(dpdr_EP)
        dpdr_EP_max = maximum(dpdr_EP_abs)
        dpdr_EP_max_loc = argmax(dpdr_EP_abs)
        n_at_max = ni[dpdr_EP_max_loc]
        if (Options.PROCESS_IN != 5)
            for ir = 1:Options.SCAN_N
                # Options.FACTOR = Options.FACTOR*dpdr_EP_max/dpdr_EP_abs[Options.IR_EXP[ir]] 
                # matches fortran
                Options.FACTOR[ir] = Options.FACTOR[ir]*dpdr_EP_max/dpdr_EP_abs[Options.IR_EXP[ir]] 
            end
        end
        Options.FACTOR_MAX_PROFILE .= Options.FACTOR
    end

    Options.F_REAL .= 1.0
    if (Options.REAL_FREQ == 1) 
        Options.F_REAL .= (cs[:]/(rmin_ex[profile.NR]))/(2*pi*1.0e3)
    end

    # deepcopy is required so as to avoid overwriting of data:
    n_ir = Options.SCAN_N
    Ts = fill(Options, n_ir)
    Ts[1] = deepcopy(Options)
    for i in 2:n_ir
        Ts[i] = deepcopy(Ts[i-1])
    end
    arrTGLFEP = Ts
    arrMTGLF = Vector{typeof(profile)}(undef, n_ir)
    arrMTGLF[1] = deepcopy(profile)
    for i in 2:n_ir
        arrMTGLF[i] = deepcopy(arrMTGLF[i-1])
    end
    arrgrowth = fill(fill(NaN,(5, 10, 10, Options.NMODES)), n_ir)

    par = _resolve_runTHD_parallel(parallel)
    if par === :threads
        stdout_lock = ReentrantLock()
        Threads.@threads for i in 1:n_ir
            _runTHD_radius!(i, arrTGLFEP, arrMTGLF, arrgrowth, printout, use_gpu; stdout_lock=stdout_lock)
        end
    elseif par === :distributed
        pmap_outputs = pmap(i -> begin
            ep = deepcopy(arrTGLFEP[i])
            mt = deepcopy(arrMTGLF[i])
            ep.IR = ep.IR_EXP[i]
            ir = ep.IR
            ep.SUFFIX = "_r" * lpad(string(ir), 3, '0')
            ep.FACTOR_IN = ep.FACTOR[i]
            println("worker $(myid()) on $(gethostname()): i=$i ir=$ir start")
            flush(stdout)
            ret = TJLFEP.mainsub(ep, mt, printout; use_gpu=use_gpu)
            println("worker $(myid()) on $(gethostname()): i=$i ir=$ir done")
            flush(stdout)
            return ret
        end, 1:n_ir)
        results = [_unpack_mainsub!(p) for p in pmap_outputs]
        all_buffers = [p[2] for p in pmap_outputs]
        for (i, (growth, ep, mt)) in enumerate(results)
            arrgrowth[i] = growth
            arrTGLFEP[i] = ep
            arrMTGLF[i] = mt
        end
        if printout
            for i in 1:n_ir
                sf_buf, wf_buf_all = all_buffers[i]
                suffix_i = coalesce(arrTGLFEP[i].SUFFIX, "")
                if sf_buf !== nothing && !isempty(sf_buf)
                    open("out.scalefactor" * suffix_i, "w") do io
                        for line in sf_buf
                            println(io, line)
                        end
                    end
                end
                if wf_buf_all !== nothing && !isempty(wf_buf_all)
                    for (str_wf_file, wfbuf) in wf_buf_all
                        if wfbuf !== nothing && !isempty(wfbuf)
                            open(str_wf_file, "w") do io
                                for line in wfbuf
                                    println(io, line)
                                end
                            end
                        end
                    end
                end
            end
        end
    else
        error("parallel must be :auto, :threads, or :distributed (got $parallel)")
    end

    # IS was set on each arrMTGLF[i] deepcopy by TJLF_map; copy it back to the
    # original profile which is used below for indexing (e.g. profile.AS[..., profile.IS])
    profile.IS = arrMTGLF[1].IS

    Options = arrTGLFEP[1]
    
    kymark_out::Vector{Float64} = fill(NaN, Options.SCAN_N)
    width::Vector{Float64} = fill(NaN, Options.SCAN_N)

    if (!Options.WIDTH_IN_FLAG)
        # Non-MPI:
        # There are only "3" processes in Threads -- 
        for i = 1:n_ir
            width[i] = arrTGLFEP[i].WIDTH_IN
            kymark_out[i] = arrTGLFEP[i].KYMARK
            
        end
    end

    # Options = arrTGLFEP[1]
    outTGLFEP_buffer = String[]
    if (printout)
        push!(outTGLFEP_buffer, "process_in = $(Options.PROCESS_IN)")
        if (Options.PROCESS_IN <= 1) push!(outTGLFEP_buffer, "mode_in = $(Options.MODE_IN)") end
        if ((Options.PROCESS_IN == 4) || (Options.PROCESS_IN == 5)) push!(outTGLFEP_buffer, "threshold_flag = $(Options.THRESHOLD_FLAG)") end
        push!(outTGLFEP_buffer, "ky_mode = $(Options.KY_MODEL)")
        push!(outTGLFEP_buffer, "--------------------------------------------------------------")
        push!(outTGLFEP_buffer, "scan_n = $(Options.SCAN_N)")
        push!(outTGLFEP_buffer, "irs = $(Options.IRS)")
        push!(outTGLFEP_buffer, "n_basis = $(Options.N_BASIS)")
        push!(outTGLFEP_buffer, "scan_method = $(Options.SCAN_METHOD)")
        if (Options.WIDTH_IN_FLAG)
            push!(outTGLFEP_buffer, "ir,  width")
            for i = 1:Options.SCAN_N
                push!(outTGLFEP_buffer, "$(Options.IRS+i-1) $(width[i])")
            end
        else
            push!(outTGLFEP_buffer, "ir,  width,  kymark")
            for i = 1:Options.SCAN_N
                push!(outTGLFEP_buffer, "$(Options.IRS+i-1) $(width[i]) $(kymark_out[i])")
            end
        end
        push!(outTGLFEP_buffer, "--------------------------------------------------------------")
        push!(outTGLFEP_buffer, "factor_in_profile = $(Options.FACTOR_IN_PROFILE)")
        if (Options.FACTOR_IN_PROFILE)
            for i = 1:Options.SCAN_N
                push!(outTGLFEP_buffer, string(Options.FACTOR[i]))
            end
        else
            push!(outTGLFEP_buffer, string(Options.FACTOR[1]))
        end
        push!(outTGLFEP_buffer, "width_in_flag = $(Options.WIDTH_IN_FLAG)")
        if (!Options.WIDTH_IN_FLAG) push!(outTGLFEP_buffer, "width_min = $(Options.WIDTH_MIN) width_max = $(Options.WIDTH_MAX)") end
    end
    
    # Initialize output arrays
    # kymark_out and width already defined above
    
    # Now continue on to radii-dependent part:
    if ((Options.PROCESS_IN == 4) || (Options.PROCESS_IN == 5))
        if (printout)
            push!(outTGLFEP_buffer, "**************************************************************")
            push!(outTGLFEP_buffer, "************** The critical EP density gradient **************")
            push!(outTGLFEP_buffer, "**************************************************************")
        end

        SFmin = fill(0.0, Options.SCAN_N)
        SFmin_out = fill(0.0, profile.NR)
        dndr_crit = fill(NaN, Options.SCAN_N)
        dndr_crit_out = fill(NaN, profile.NR)
        dpdr_crit = fill(NaN, Options.SCAN_N)
        dpdr_crit_out = fill(NaN, profile.NR)

        if (Options.THRESHOLD_FLAG == 0)
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
                push!(outTGLFEP_buffer, "--------------------------------------------------------------")
                push!(outTGLFEP_buffer, "SFmin")
            end
        
        # Next is TGLFEP_complete_output(SFmin, SFmin_out, ir_min, ir_max, l_accept_profile)
        # This function's goal is to determine whether 

            SFmin, SFmin_out, ir_min, ir_max, l_accept_profile = tjlfep_complete_output(SFmin, Options, profile)
        
            if (printout)
                push!(outTGLFEP_buffer, string(SFmin, " SFmin after buf and coutput")) # after comp.out
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
                push!(outTGLFEP_buffer, string(SFmin, " SFmin after Max assign"))
            end

            if (printout)
                for i = 1:Options.SCAN_N
                    if (l_accept_profile[i])
                        push!(outTGLFEP_buffer, string(SFmin[i]))
                    else
                        push!(outTGLFEP_buffer, string(SFmin_out[i+Options.IRS-1], "   (*)"))
                    end
                    push!(outTGLFEP_buffer, "--------------------------------------------------------------")
                    push!(outTGLFEP_buffer, string(l_accept_profile))
                    push!(outTGLFEP_buffer, string(ir_min, " : ", ir_max))
                end
            end

            # Calculate the density critical gradient at each of the scanned radii.
            if (Options.INPUT_PROFILE_METHOD == 2)
                dndr_crit .= 10000.0
                for i = 1:Options.SCAN_N
                    # If SFmin[i] is not the default non-rejected value, multiply the scalefactor by the density and the density gradient at that point 
                    # for the energetic ion. 
                    # If SFmin[i] is the default or >= 9k, check if it is one of the factor_max_profile ones, and if so, calculate it with that.
                    # otherwise, leave it at 10k.
                    if (SFmin[i] < 9000.0)
                        dndr_crit[i] = SFmin[i]*ni[Int(Options.IR_EXP[i])]*dlnnidr[Int(Options.IR_EXP[i])]
                    elseif ((i < ir_min-Options.IRS+1) || (i > ir_max-Options.IRS+1))
                        dndr_crit[i] = Options.FACTOR_MAX_PROFILE[i]*ni[Int(Options.IR_EXP[i])]*dlnnidr[Int(Options.IR_EXP[i])]
                    end
                end
                # Interpolate and accept are reject needed values of this profile:
                dndr_crit, dndr_crit_out, ir_dum_1, ir_dum_2, l_accept_profile = tjlfep_complete_output(dndr_crit, Options, profile)
                
                if (printout)
                    io4 = open("alpha_dndr_crit.input", "w")
                    println(io4, "Density critical gradient (10^19/m^4)")
                    println(io4, dndr_crit_out)
                    close(io4)
                end
            end

            if (Options.INPUT_PROFILE_METHOD == 2)
                dpdr_crit .= 10000.0
                nr = profile.NR
                for i in 1:nr
                    dpdr_EP[i] = ni[i] * Ti[i] * (dlnnidr[i] + dlntidr[i]) * 0.16022
                end
                for i = 1:Options.SCAN_N
                    if (SFmin[i] < 9000.0)
                        if ((Options.PROCESS_IN == 4) || (Options.PROCESS_IN == 5))
                            case = Options.SCAN_METHOD
                            if (case == 1)
                                dpdr_scale = SFmin[i]
                            elseif (case == 2)
                                dpdr_scale = ((SFmin[i]*dlnnidr[Options.IR_EXP[i]]+dlntidr[Options.IR_EXP[i]]) /
                                (dlnnidr[Options.IR_EXP[i]]+dlntidr[Options.IR_EXP[i]]))
                            end
                            dpdr_crit[i] = dpdr_scale*dpdr_EP[Options.IR_EXP[i]]
                        end # 4 || 5
                    end # < 9000
                end # over scan_n
                dpdr_crit, dpdr_crit_out, ir_dum_1, ir_dum_2, l_accept_profile = tjlfep_complete_output(dpdr_crit, Options, profile)
                
                if (printout)
                    io5 = open("alpha_dpdr_crit.input", "w")
                    println(io5, "Pressure critical gradient (10 kPa/m)")
                    println(io5, dpdr_crit_out)
                    close(io5)
                end
            end # end prof. method 2

            if (printout)
                push!(outTGLFEP_buffer, "--------------------------------------------------------------")
                push!(outTGLFEP_buffer, "The EP density threshold n_EP/n_e (%) for gamma_AE = 0")
                for i = 1:Options.SCAN_N
                    push!(outTGLFEP_buffer, string(SFmin[i]*profile.AS[Options.IRS+i-1, profile.IS]*100.0))
                end
            end
            
            if (printout)
                push!(outTGLFEP_buffer, "--------------------------------------------------------------")
                push!(outTGLFEP_buffer, "The EP beta crit (%) = beta_e*(n_EP_th/n_e)*(T_EP/T_e)")
                for i = 1:Options.SCAN_N
                    if coalesce(profile.GEOMETRY_FLAG, 1) == 0
                        push!(outTGLFEP_buffer, string(SFmin[i]*profile.BETAE[Options.IRS+i-1]*100.0*profile.AS[Options.IRS+i-1, profile.IS]*profile.TAUS[Options.IRS+i-1, profile.IS]))
                    else
                        push!(outTGLFEP_buffer, string(SFmin[i]*profile.BETAE[Options.IRS+i-1]*100.0*profile.AS[Options.IRS+i-1, profile.IS]*profile.TAUS[Options.IRS+i-1, profile.IS]*profile.KAPPA[Options.IRS+i-1]^2))
                    end
                end
            end

            # there is a process_in == 4 addition I won't be doing quite yet.
        else # ThreshFlag != 0
            # Skipping for now as I want to test just threshold flag == 0 first
        end # ThreshFlag
    end # process 4 || 5

    # At the very end, write the buffer to file once
    if (printout)
        open("out.TGLFEP", "w") do io
            for line in outTGLFEP_buffer
                println(io, line)
            end
        end
    end

    return width, kymark_out, SFmin, dpdr_crit_out, dndr_crit_out
end  # End of string-based runTHD
"""
checkInput(inputTJLF::InputTJLF)

description:
check that the InputTJLF struct is properly populated
"""
function checkInput(inputTJLFEP::InputTJLF)
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

function checkInput(inputTJLFEPVector::Vector{InputTJLF})
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