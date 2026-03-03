Base.@kwdef mutable struct InputTGLF
    SIGN_BT::Union{Int,Missing} = missing
    SIGN_IT::Union{Int,Missing} = missing
    NS::Union{Int,Missing} = missing
    ZMAJ_LOC::Union{Float64,Missing} = missing
    DRMINDX_LOC::Union{Float64,Missing} = missing
    DZMAJDX_LOC::Union{Float64,Missing} = missing
    S_DELTA_LOC::Union{Float64,Missing} = missing
    ZETA_LOC::Union{Float64,Missing} = missing
    S_ZETA_LOC::Union{Float64,Missing} = missing

    MASS_1::Union{Float64,Missing} = missing
    ZS_1::Union{Float64,Missing} = missing
    AS_1::Union{Float64,Missing} = missing
    TAUS_1::Union{Float64,Missing} = missing

    MASS_2::Union{Float64,Missing} = missing
    ZS_2::Union{Float64,Missing} = missing
    VPAR_2::Union{Float64,Missing} = missing
    VPAR_SHEAR_2::Union{Float64,Missing} = missing

    MASS_3::Union{Float64,Missing} = missing
    ZS_3::Union{Float64,Missing} = missing
    RLTS_3::Union{Float64,Missing} = missing
    TAUS_3::Union{Float64,Missing} = missing
    VPAR_3::Union{Float64,Missing} = missing
    VPAR_SHEAR_3::Union{Float64,Missing} = missing

    # TGLF-NN uses 3 species
    # This is why parameters for species 1:3 are sorted differently than 4:10
    MASS_4::Union{Float64,Missing} = missing
    AS_4::Union{Float64,Missing} = missing
    ZS_4::Union{Float64,Missing} = missing
    RLNS_4::Union{Float64,Missing} = missing
    RLTS_4::Union{Float64,Missing} = missing
    TAUS_4::Union{Float64,Missing} = missing
    VPAR_4::Union{Float64,Missing} = missing
    VPAR_SHEAR_4::Union{Float64,Missing} = missing

    MASS_5::Union{Float64,Missing} = missing
    AS_5::Union{Float64,Missing} = missing
    ZS_5::Union{Float64,Missing} = missing
    RLNS_5::Union{Float64,Missing} = missing
    RLTS_5::Union{Float64,Missing} = missing
    TAUS_5::Union{Float64,Missing} = missing
    VPAR_5::Union{Float64,Missing} = missing
    VPAR_SHEAR_5::Union{Float64,Missing} = missing

    MASS_6::Union{Float64,Missing} = missing
    AS_6::Union{Float64,Missing} = missing
    ZS_6::Union{Float64,Missing} = missing
    RLNS_6::Union{Float64,Missing} = missing
    RLTS_6::Union{Float64,Missing} = missing
    TAUS_6::Union{Float64,Missing} = missing
    VPAR_6::Union{Float64,Missing} = missing
    VPAR_SHEAR_6::Union{Float64,Missing} = missing

    MASS_7::Union{Float64,Missing} = missing
    AS_7::Union{Float64,Missing} = missing
    ZS_7::Union{Float64,Missing} = missing
    RLNS_7::Union{Float64,Missing} = missing
    RLTS_7::Union{Float64,Missing} = missing
    TAUS_7::Union{Float64,Missing} = missing
    VPAR_7::Union{Float64,Missing} = missing
    VPAR_SHEAR_7::Union{Float64,Missing} = missing

    MASS_8::Union{Float64,Missing} = missing
    AS_8::Union{Float64,Missing} = missing
    ZS_8::Union{Float64,Missing} = missing
    RLNS_8::Union{Float64,Missing} = missing
    RLTS_8::Union{Float64,Missing} = missing
    TAUS_8::Union{Float64,Missing} = missing
    VPAR_8::Union{Float64,Missing} = missing
    VPAR_SHEAR_8::Union{Float64,Missing} = missing

    MASS_9::Union{Float64,Missing} = missing
    AS_9::Union{Float64,Missing} = missing
    ZS_9::Union{Float64,Missing} = missing
    RLNS_9::Union{Float64,Missing} = missing
    RLTS_9::Union{Float64,Missing} = missing
    TAUS_9::Union{Float64,Missing} = missing
    VPAR_9::Union{Float64,Missing} = missing
    VPAR_SHEAR_9::Union{Float64,Missing} = missing

    MASS_10::Union{Float64,Missing} = missing
    AS_10::Union{Float64,Missing} = missing
    ZS_10::Union{Float64,Missing} = missing
    RLNS_10::Union{Float64,Missing} = missing
    RLTS_10::Union{Float64,Missing} = missing
    TAUS_10::Union{Float64,Missing} = missing
    VPAR_10::Union{Float64,Missing} = missing
    VPAR_SHEAR_10::Union{Float64,Missing} = missing

    AS_2::Union{Float64,Missing} = missing
    AS_3::Union{Float64,Missing} = missing
    BETAE::Union{Float64,Missing} = missing
    DEBYE::Union{Float64,Missing} = missing
    DELTA_LOC::Union{Float64,Missing} = missing
    DRMAJDX_LOC::Union{Float64,Missing} = missing
    KAPPA_LOC::Union{Float64,Missing} = missing
    P_PRIME_LOC::Union{Float64,Missing} = missing
    Q_LOC::Union{Float64,Missing} = missing
    Q_PRIME_LOC::Union{Float64,Missing} = missing
    RLNS_1::Union{Float64,Missing} = missing
    RLNS_2::Union{Float64,Missing} = missing
    RLNS_3::Union{Float64,Missing} = missing
    RLTS_1::Union{Float64,Missing} = missing
    RLTS_2::Union{Float64,Missing} = missing
    RMAJ_LOC::Union{Float64,Missing} = missing
    RMIN_LOC::Union{Float64,Missing} = missing
    S_KAPPA_LOC::Union{Float64,Missing} = missing
    TAUS_2::Union{Float64,Missing} = missing
    VEXB_SHEAR::Union{Float64,Missing} = missing
    VPAR_1::Union{Float64,Missing} = missing
    VPAR_SHEAR_1::Union{Float64,Missing} = missing
    XNUE::Union{Float64,Missing} = missing
    ZEFF::Union{Float64,Missing} = missing

    # switches
    UNITS::Union{String,Missing} = missing
    ALPHA_ZF::Union{Float64,Missing} = missing
    USE_MHD_RULE::Union{Bool,Missing} = missing
    NKY::Union{Int,Missing} = missing
    SAT_RULE::Union{Int,Missing} = missing
    KYGRID_MODEL::Union{Int,Missing} = missing
    NMODES::Union{Int,Missing} = missing
    NBASIS_MIN::Union{Int,Missing} = missing
    NBASIS_MAX::Union{Int,Missing} = missing
    XNU_MODEL::Union{Int,Missing} = missing
    USE_AVE_ION_GRID::Union{Bool,Missing} = missing
    ALPHA_QUENCH::Union{Int,Missing} = missing
    ALPHA_MACH::Union{Float64,Missing} = missing
    WDIA_TRAPPED::Union{Float64,Missing} = missing
    USE_BPAR::Union{Bool,Missing} = missing
    USE_BPER::Union{Bool,Missing} = missing

    _Qgb::Union{Float64,Missing} = missing

    # missing
    USE_BISECTION::Bool = true
    USE_INBOARD_DETRAPPED::Bool = false
    NEW_EIKONAL::Bool = true
    FIND_WIDTH::Bool = true
    IFLUX::Bool = true
    ADIABATIC_ELEC::Bool = false

    NWIDTH::Int = 21
    NXGRID::Int = 16
    VPAR_MODEL::Int = 0
    VPAR_SHEAR_MODEL::Int = 1
    IBRANCH::Int = -1

    KY::Float64 = 0.3
    ALPHA_E::Float64 = 1.0
    ALPHA_P::Float64 = 1.0
    XNU_FACTOR::Float64 = 1.0
    DEBYE_FACTOR::Float64 = 1.0
    RLNP_CUTOFF::Float64 = 18.0
    WIDTH::Float64 = 1.65
    WIDTH_MIN::Float64 = 0.3
    BETA_LOC::Float64 = 1.0
    KX0_LOC::Float64 = 1.0
    PARK::Float64 = 1.0
    GHAT::Float64 = 1.0
    GCHAT::Float64 = 1.0
    WD_ZERO::Float64 = 0.1
    LINSKER_FACTOR::Float64 = 0.0
    GRADB_FACTOR::Float64 = 0.0
    FILTER::Float64 = 2.0
    THETA_TRAPPED::Float64 = 0.7
    ETG_FACTOR::Float64 = 1.25
    DAMP_PSI::Float64 = 0.0
    DAMP_SIG::Float64 = 0.0

end


mutable struct InputTJLF{T<:Real}

    UNITS::Union{String,Missing}

    USE_BPER::Union{Bool,Missing}
    USE_BPAR::Union{Bool,Missing}
    USE_MHD_RULE::Union{Bool,Missing}
    USE_BISECTION::Union{Bool,Missing}
    USE_INBOARD_DETRAPPED::Union{Bool,Missing}
    USE_AVE_ION_GRID::Union{Bool,Missing}
    NEW_EIKONAL::Union{Bool,Missing} ## this seems useless, the flag has to be both true and false to do anything
    FIND_WIDTH::Union{Bool,Missing}
    IFLUX::Union{Bool,Missing}
    ADIABATIC_ELEC::Union{Bool,Missing}

    SAT_RULE::Union{Int,Missing}
    NS::Union{Int,Missing}
    NMODES::Union{Int,Missing}
    NWIDTH::Union{Int,Missing}
    NBASIS_MAX::Union{Int,Missing}
    NBASIS_MIN::Union{Int,Missing}
    NXGRID::Union{Int,Missing}
    NKY::Union{Int,Missing}
    KYGRID_MODEL::Union{Int,Missing}
    XNU_MODEL::Union{Int,Missing}
    VPAR_MODEL::Union{Int,Missing}
    IBRANCH::Union{Int,Missing}

    ZS::Union{Vector{T},Missing}
    MASS::Union{Vector{T},Missing}
    RLNS::Union{Vector{T},Missing}
    RLTS::Union{Vector{T},Missing}
    TAUS::Union{Vector{T},Missing}
    AS::Union{Vector{T},Missing}
    VPAR::Union{Vector{T},Missing}
    VPAR_SHEAR::Union{Vector{T},Missing}

    # NOT IN TGLF
    WIDTH_SPECTRUM::Union{Vector{T},Missing}
    KY_SPECTRUM::Union{Vector{T},Missing}
    EIGEN_SPECTRUM::Union{Vector{ComplexF64},Missing}
    FIND_EIGEN::Union{Bool,Missing}
    # NOT IN TGLF

    SIGN_BT::Union{Int,Missing}
    SIGN_IT::Union{Int,Missing}
    KY::Union{T,Missing}

    VEXB_SHEAR::Union{T,Missing}
    BETAE::Union{T,Missing}
    XNUE::Union{T,Missing}
    ZEFF::Union{T,Missing}
    DEBYE::Union{T,Missing}

    ALPHA_MACH::Union{T,Missing}
    ALPHA_E::Union{T,Missing}
    ALPHA_P::Union{T,Missing}
    ALPHA_QUENCH::Union{Int,Missing}
    ALPHA_ZF::Union{T,Missing}
    XNU_FACTOR::Union{T,Missing}
    DEBYE_FACTOR::Union{T,Missing}
    ETG_FACTOR::Union{T,Missing}
    RLNP_CUTOFF::Union{T,Missing}

    WIDTH::Union{T,Missing}
    WIDTH_MIN::Union{T,Missing}

    RMIN_LOC::Union{T,Missing}
    RMAJ_LOC::Union{T,Missing}
    ZMAJ_LOC::Union{T,Missing}
    DRMINDX_LOC::Union{T,Missing}
    DRMAJDX_LOC::Union{T,Missing}
    DZMAJDX_LOC::Union{T,Missing}
    Q_LOC::Union{T,Missing}
    KAPPA_LOC::Union{T,Missing}
    S_KAPPA_LOC::Union{T,Missing}
    DELTA_LOC::Union{T,Missing}
    S_DELTA_LOC::Union{T,Missing}
    ZETA_LOC::Union{T,Missing}
    S_ZETA_LOC::Union{T,Missing}
    P_PRIME_LOC::Union{T,Missing}
    Q_PRIME_LOC::Union{T,Missing}
    BETA_LOC::Union{T,Missing}
    KX0_LOC::Union{T,Missing}
    DAMP_PSI::Union{T,Missing}
    DAMP_SIG::Union{T,Missing}
    WDIA_TRAPPED::Union{T,Missing}
    PARK::Union{T,Missing}
    GHAT::Union{T,Missing}
    GCHAT::Union{T,Missing}
    WD_ZERO::Union{T,Missing}
    LINSKER_FACTOR::Union{T,Missing}
    GRADB_FACTOR::Union{T,Missing}
    FILTER::Union{T,Missing}
    THETA_TRAPPED::Union{T,Missing}
    SMALL::Union{T,Missing}

    USE_TRANSPORT_MODEL::Union{Bool, Missing}

    #For list-format inputs:
    # function InputTJLF{T}(inP::Bool) where {T<:Real}
    #     if inP
    #         new("CGYRO", false, false, false, true, false, true, true, true, true, false, 2, 3, 5, 21, 6, 4, 16, 12, 4, 3, 0, -1, [-1.0, 1.0, 6.0], [0.0002723125672605524, 1.0, 6.0], [0.9691383387573976, 1.078021414201318, 0.0733427635614379], [3.332037619158914, 2.0626412607995435, 2.0626412607995435], [1.0, 1.3661261082028286, 1.3661261082028286], [1.0, 0.8075398023805694, 0.030988644410732645], [0.30611236015079274, 0.30611236015079274, 0.30611236015079274], [1.5491649356389778, 1.5491649356389778, 1.5491649356389778], [1.65, 1.65, 1.65, 1.65, 0.9685467847385054, 0.7035623639735143, 0.6324554384084142, 0.591251466897806, 0.5250292077902518, 1.65, 1.65, 1.65, 1.65, 1.65, 1.65, 1.65, 1.551199191728457, 1.5603157485179238, 1.5643992899672103, 1.5643992899672103, 1.4832394778484315, 0.6324554384084142, 0.5390399200630244, 0.4845607815598427], [0.05994688615887238, 0.11989377231774476, 0.17984065847661712, 0.2397875446354895, 0.2997344307943619, 0.5395219754298515, 0.6594157477475961, 0.7793095200653409, 0.8992032923830857, 1.0190970647008304, 1.138990837018575, 1.25888460933632, 1.1989377231774476, 1.60763400555729, 2.1556474918269477, 2.890468908319081, 3.875777714879744, 5.196960552637156, 6.968510831252536, 9.343950702231606, 12.529135254288084, 16.80009186935815, 22.52694069387401, 30.205969167995068], ComplexF64[0.01903754432811539 - 0.03822460700618263im, 0.066744785730153 - 0.08638900959186772im, 0.12700366682079575 - 0.13466694388477374im, 0.17881630543473276 - 0.16774618757919235im, 0.21216108783162443 - 0.18009019917758487im, 0.33316817373508417 - 0.3198305942745321im, 0.34745089281285046 - 0.3955067169459169im, 0.3399817917167648 - 0.4616241685193869im, 0.3137572066001085 - 0.5077973469905399im, 0.28555930064941276 + 0.45333962118452586im, 0.34112248363191056 + 0.501904780504562im, 0.3914434522233301 + 0.5529917255055463im, 0.36420457375571014 + 0.5235784686127289im, 0.531883979920564 + 0.7034735644682284im, 0.7500443094115281 + 0.9612053548664551im, 1.0612410507865202 + 1.3223708279218578im, 1.2154632046734721 + 1.7048988962925355im, 1.8556234007580341 + 2.2437573929973844im, 2.891191856550966 + 3.1140102214419576im, 4.213253106128079 + 4.483889088628243im, 5.747152586787458 + 6.01284811393025im, 7.352163708936412 + 7.770351222625782im, 9.652984589008142 + 10.27500211396812im, 11.970055852084016 + 14.27560461246188im], true, -1, 1, 0.3, 0.148365431821359, 0.0009809454014984833, 0.2658337070903717, 1.9296593, 0.029821537734289975, 0.0, 1.0, 1.0, 0, -1.0, 1.0, 1.0, 1.25, 18.0, 1.65, 0.3, 0.8896452200962354, 2.8058920740841784, 0.0, 1.0, -0.19752155788650919, 0.0, 3.3106313319155714, 1.6054967596315595, 0.39307251195418547, 0.21740011375976812, 0.7746695322421236, -0.05113765116526302, -0.2388377806334241, -0.0011605489188390146, 35.13388509054382, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.1, 0.0, 0.0, 2.0, 0.7, 1.0e-12)        else
    #     end
    # end

    function InputTJLF{T}(ns::Int, nky::Int, dflt::Bool) where {T<:Real}
        if dflt
            new("GYRO",
            false,false,true,true,false,missing,true,false,true,false,
            0,2,2,21,4,4,32,12,0,2,0,-1,
            fill(NaN,(ns)),fill(NaN,(ns)),fill(NaN,(ns)),fill(NaN,(ns)),fill(NaN,(ns)),fill(NaN,(ns)),fill(NaN,(ns)),fill(NaN,(ns)),
            fill(NaN,(nky)),fill(NaN,(nky)),fill(NaN*im,(nky)),missing,
            1.0,1.0,0.3,0.0,0.0,0.0,1.0,0.0,0.0,1.0,
            1.0,0.0,1.0,1.0,1.0,1.25,18.0,1.65,0.3,0.5,
            3.0,0.0,1.0,0.0,0.0,2.0,1.0,16.0,0.0,0.0,
            0.0,0.0,0.0,16.0,0.0,0.0,0.0,0.0,0.0,1.0,
            1.0,1.0,0.1,0.0,0.0,0.0,0.7,1.0e-13,true)
        else
            new("",
            missing,missing,missing,missing,missing,missing,missing,missing,missing,missing,
            missing,missing,missing,missing,missing,missing,missing,missing,missing,missing,missing,missing,
            fill(NaN,(ns)),fill(NaN,(ns)),fill(NaN,(ns)),fill(NaN,(ns)),fill(NaN,(ns)),fill(NaN,(ns)),fill(NaN,(ns)),fill(NaN,(ns)),
            fill(NaN,(nky)),fill(NaN,(nky)),fill(NaN*im,(nky)),missing,
            0.0,0.0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,
            NaN,0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,
            NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,
            NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,
            NaN,NaN,NaN,NaN,NaN,NaN,NaN,1.0e-13,true)
        end
    end

    # function InputTJLF{T}(ns::Int, nky::Int) where {T<:Real}
    #     InputTJLF{T}(ns, nky, false)
    # end

    # create InputTJLF struct given a InputTGLF struct
    # function InputTJLF{T}(inputTGLF::InputTGLF) where {T<:Real}

    #     nky = get_ky_spectrum_size(inputTGLF.NKY,inputTGLF.KYGRID_MODEL)
    #     inputTJLF = InputTJLF{T}(inputTGLF.NS, nky)

    #     for fieldname in fieldnames(inputTGLF)
    #         if occursin(r"\d", String(fieldname)) || fieldname == :_Qgb # species parameter
    #             continue
    #         end
    #         setfield!(inputTJLF, fieldname, getfield(inputTGLF, fieldname))
    #     end
    #     for i in 1:inputTGLF.NS
    #         inputTJLF.ZS[i] = getfield(inputTGLF, Symbol("ZS_", i))
    #         inputTJLF.AS[i] = getfield(inputTGLF, Symbol("AS_", i))
    #         inputTJLF.MASS[i] = getfield(inputTGLF, Symbol("MASS_", i))
    #         inputTJLF.RLNS[i] = getfield(inputTGLF, Symbol("RLNS_", i))
    #         inputTJLF.RLTS[i] = getfield(inputTGLF, Symbol("RLTS_", i))
    #         inputTJLF.TAUS[i] = getfield(inputTGLF, Symbol("TAUS_", i))
    #         inputTJLF.VPAR[i] = getfield(inputTGLF, Symbol("VPAR_", i))
    #         inputTJLF.VPAR_SHEAR[i] = getfield(inputTGLF, Symbol("VPAR_SHEAR_", i))
    #     end
    #     inputTJLF.WIDTH_SPECTRUM .= inputTJLF.WIDTH

    #     checkInput(inputTJLF)

    #     return inputTJLF
    # end
end

mutable struct Options{T<:Real} # This acts as the interface module of Fortran, essentially. It reads the TGLFEP file
    PROCESS_IN::Union{Int, Missing} # May need to be a MODE_IN::Union{T, Missing} for PROCESS_IN <= 1
    THRESHOLD_FLAG::Union{Int, Missing}
    N_BASIS::Union{Int, Missing}
    SCAN_METHOD::Union{Int, Missing}
    REJECT_I_PINCH_FLAG::Union{Int, Missing} # 5
    REJECT_E_PINCH_FLAG::Union{Int, Missing}
    REJECT_TH_PINCH_FLAG::Union{Int, Missing}
    REJECT_EP_PINCH_FLAG::Union{Int, Missing}
    REJECT_TEARING_FLAG::Union{Int, Missing}
    ROTATIONAL_SUPPRESSION_FLAG::Union{Int, Missing} # 10
    QL_RATIO_THRESH::Union{T, Missing}
    THETA_SQ_THRESH::Union{T, Missing}
    Q_SCALE::Union{T, Missing}
    WRITE_WAVEFUNCTION::Union{Int, Missing} # Maybe not needed for my purposes? Not sure yet
    KY_MODEL::Union{Int, Missing} #15
    SCAN_N::Union{Int} 
    IRS::Union{Int, Missing}
    FACTOR_IN_PROFILE::Union{Bool, Missing}
    FACTOR_IN::Union{T, Missing}
    WIDTH_IN_FLAG::Union{Bool} #20
    WIDTH_IN::Union{T, Missing} 
    WIDTH_MIN::Union{T, Missing}
    WIDTH_MAX::Union{T, Missing}
    INPUT_PROFILE_METHOD::Union{Int, Missing} # I very likely will not use this (?). It seems that INPUT_PROFILE_METHOD is meant for 
    # determining the type of input that will be used: input.profile or input.profiles. EXPRO is the method used to do the latter
    # but I cannot find the runction in Fortran for expro_read because it is imported from somewhere else. 
    IS_EP::Union{Int, Missing} #25
    M_I::Union{T, Missing} 
    Q_I::Union{T, Missing}
    M_EP::Union{T, Missing}
    Q_EP::Union{T, Missing}
    IR::Union{Int, Missing} #30
    NMODES::Union{Int}  
    MODE_IN::Union{Int, Missing} # Likely will not be used for now. There are some default settings dependent on the choice of process-in. = 2 for process_in = 5: the case we are interested in.
    FACTOR_MAX_PROFILE::Union{Vector{T}, Missing}
    FACTOR_MAX::Union{T, Missing}
    KY_IN::Union{T, Missing} #35
    KYMARK::Union{T, Missing} 
    # All of the next 7 vectors require the value nr to be recognized for the array. This is already done in the profile struct.
    # I need to decide what I want to do with this then. These are for process_in = 6 apparently.
    DEP_TRACE_LOC::Union{T, Missing}
    DEP_TRACE::Union{Vector{T}, Missing}
    DEP_TRACE_COMPLETE::Union{Vector{T}, Missing}
    DEP_TRACE_SCALE::Union{Matrix{T}, Missing} #40
    QL_RATIO_LOC::Union{T, Missing} 
    QL_RATIO::Union{Vector{T}, Missing}
    QL_RATIO_COMPLETE::Union{Vector{T}, Missing}
    QL_RATIO_SCAN::Union{Matrix{T}, Missing}
    CHI_GB::Union{Vector{T}, Missing} #45
    IR_EXP::Union{Vector{Int64}, Missing}  
    NBASIS::Union{Int, Missing}
    NTOROIDAL::Union{Int, Missing}
    KYHAT_IN::Union{T, Missing}
    SCAN_FACTOR::Union{T, Missing} #50 Default to 1
    JTSCALE::Union{Int, Missing} 
    JTSCALE_MAX::Union{Int, Missing} # Default to 1
    TSCALE_INTERVAL::Union{T, Missing} # Default to 1.0
    N_ION::Union{Int, Missing} # Equal to NS-1?
    FREQ_CUTOFF::Union{T, Missing} #55 Default to -0.35
    FREQ_AE_UPPER::Union{T, Missing} 
    GAMMA_THRESH::Union{T,Missing}   
    GAMMA_THRESH_MAX::Union{T, Missing} #58
    NN::Union{Int, Missing}
    FACTOR_OUT::Union{Vector{T}, Missing} #60 vector of nn elements
    # May or may not need id_2, np_2, id_3, np_3 ?
    # Ignoring suffix and str_wf_file right now.
    # As with l_print, l-WRITE_WAVEFUNCTION, l_wavefunction_out
    LKEEP::Vector{Bool} # All the L-starting parts have 4 dimensions. See kwscale_scan for that.
    LTEARING::Vector{Bool} # This needs to be defined somewhat differently then.
    L_TH_PINCH::Vector{Bool} 
    L_I_PINCH::Vector{Bool}
    L_E_PINCH::Vector{Bool}#65
    L_EP_PINCH::Vector{Bool} 
    L_QL_RATIO::Vector{Bool}
    L_THETA_SQ::Vector{Bool} 
    FACTOR::Union{Vector{T}, Missing}
    SUFFIX::Union{String, Missing} #70
    F_REAL::Union{Vector{T}, Missing} 
    REAL_FREQ::Int
    WIDTH::Union{Vector{T}, Missing}#73
    # I might want a default constructor for default values. I'm going to make one (no paramters)
    # which follows the 'default' case that TGLFEP used in my first test of TGLFEP (OMFIT).

    # There are a bunch of cases where the dimension of the fills or the existence of a value (WIDTH_IN and WIDTH_MIN/MAX)
    # depends on other input variables:

    # If WIDTH_IN_FLAG is true the fill needs (?) to be SCAN_N long of that same value and the min/max are missing
    # There is a bit of confusion on if the same goes for the FACTOR_IN_PROFILE. 

    # There are few things that depend on the profile reading: nr
    function Options{T}(nscan_in::Int64, widthin::Bool, nn::Int64, nr::Int64, jtscale_max::Int64, nmodes::Int64) where {T<:Real}
        if(widthin)
            new(missing, missing, missing, missing, missing, missing, missing, missing, missing, missing,
            NaN, NaN, NaN, missing, missing, nscan_in, missing, missing, NaN, widthin, NaN,
            NaN, NaN, missing, missing, NaN, NaN, NaN, NaN, missing, nmodes, missing, fill(NaN, nscan_in),
            NaN, NaN, NaN, NaN, fill(NaN, nr), fill(NaN, nr), fill(NaN, (jtscale_max, nr)), NaN, fill(NaN, nr), fill(NaN, nr),
            fill(NaN, (jtscale_max, nr)), fill(NaN, nr), fill(0, nscan_in), missing, missing, NaN, 1, missing, 1, 
            1.0, missing, -0.35, NaN, NaN, NaN, nn, fill(NaN, nn), fill(false, 4), fill(false, 4),
            fill(false, 4), fill(false, 4), fill(false, 4),
            fill(false, 4), fill(false, 4), fill(false, 4),
            fill(NaN, nscan_in), missing, fill(NaN, nr), 0, fill(NaN, nscan_in))
        else
            new(missing, missing, missing, missing, missing, missing, missing, missing, missing, missing,
            NaN, NaN, NaN, missing, missing, nscan_in, missing, missing, NaN, widthin, 0.0,
            NaN, NaN, missing, missing, NaN, NaN, NaN, NaN, missing, nmodes, missing, fill(NaN, nscan_in),
            NaN, NaN, NaN, NaN, fill(NaN, nr), fill(NaN, nr), fill(NaN, (jtscale_max, nr)), NaN, fill(NaN, nr), fill(NaN, nr),
            fill(NaN, (jtscale_max, nr)), fill(NaN, nr), fill(0, nscan_in), missing, missing, NaN, 1, missing, 1, 
            1.0, missing, -0.35, NaN, NaN, NaN, nn, fill(NaN, nn), fill(false, 4), fill(false, 4),
            fill(false, 4), fill(false, 4), fill(false, 4),
            fill(false, 4), fill(false, 4), fill(false, 4),
            fill(NaN, nscan_in), missing, fill(NaN, nr), 0, fill(NaN, nscan_in))
        end
    end
end

mutable struct profile{T<:Real}
    # This struct is an intermediary between the processes in read_inputs
    
    SIGN_BT::Union{T, Missing} #
    SIGN_IT::Union{T, Missing} #

    NR::Union{Int, Missing} #
    NS::Union{Int, Missing} #
    GEOMETRY_FLAG::Union{Int, Missing} #5   Constant 1 for TJLFEPnotanything else.
    ROTATION_FLAG::Union{Int, Missing} #

    ZS::Union{Vector{T}, Missing} #
    MASS::Union{Vector{T}, Missing} #

    # This next section is especially why the profile struct must exist
    AS::Union{Matrix{T}, Missing} #
    TAUS::Union{Matrix{T}, Missing} #10
    RLNS::Union{Matrix{T}, Missing} #
    RLTS::Union{Matrix{T}, Missing} #
    VPAR::Union{Matrix{T}, Missing} #
    VPAR_SHEAR::Union{Matrix{T}, Missing} #

    RMIN::Union{Vector{T}, Missing}#15
    RMAJ::Union{Vector{T}, Missing}
    SHIFT::Union{Vector{T}, Missing}
    Q::Union{Vector{T}, Missing}
    SHEAR::Union{Vector{T}, Missing}
    ALPHA::Union{Vector{T}, Missing}#20
    Q_PRIME::Union{Vector{T}, Missing}
    P_PRIME::Union{Vector{T}, Missing}
    KAPPA::Union{Vector{T}, Missing}
    S_KAPPA::Union{Vector{T}, Missing}
    DELTA::Union{Vector{T}, Missing}#25
    S_DELTA::Union{Vector{T}, Missing}
    ZETA::Union{Vector{T}, Missing}
    S_ZETA::Union{Vector{T}, Missing}
    ZEFF::Union{Vector{T}, Missing}
    BETAE::Union{Vector{T}, Missing}#30
    RHO_STAR::Union{Vector{T}, Missing}
    OMEGA_TAE::Union{Vector{T}, Missing}
    omegaGAM::Union{Vector{T}, Missing}
    gammaE::Union{Vector{T}, Missing}
    gammap::Union{Vector{T}, Missing}#35
    B_UNIT::Union{Vector{T}, Missing}

    IS::Union{Int, Missing}
    IRS::Union{Int, Missing}

    A_QN::Union{T, Missing} # quasineutrality scale factor for non-EP species
    N_ION::Union{Int, Missing} #40

    DENS_1::Union{Vector{T},Missing} 
    TEMP_1::Union{Vector{T},Missing}
    DLNNDR_1::Union{Vector{T},Missing} 
    DLNTDR_1::Union{Vector{T},Missing} #44

    DENS_2::Union{Vector{T},Missing} 
    TEMP_2::Union{Vector{T},Missing} 
    DLNNDR_2::Union{Vector{T},Missing} 
    DLNTDR_2::Union{Vector{T},Missing} #48

    DENS_3::Union{Vector{T},Missing} 
    TEMP_3::Union{Vector{T},Missing}
    DLNNDR_3::Union{Vector{T},Missing}
    DLNTDR_3::Union{Vector{T},Missing}  #52

    DENS_4::Union{Vector{T},Missing} 
    TEMP_4::Union{Vector{T},Missing} 
    DLNNDR_4::Union{Vector{T},Missing} 
    DLNTDR_4::Union{Vector{T},Missing}  #56

    DENS_5::Union{Vector{T},Missing} 
    TEMP_5::Union{Vector{T},Missing} 
    DLNNDR_5::Union{Vector{T},Missing} 
    DLNTDR_5::Union{Vector{T},Missing}  # 60

    DENS_6::Union{Vector{T},Missing} 
    TEMP_6::Union{Vector{T},Missing} 
    DLNNDR_6::Union{Vector{T},Missing} 
    DLNTDR_6::Union{Vector{T},Missing}  #64

    DENS_7::Union{Vector{T},Missing} 
    TEMP_7::Union{Vector{T},Missing} 
    DLNNDR_7::Union{Vector{T},Missing} 
    DLNTDR_7::Union{Vector{T},Missing}  #68

    DENS_8::Union{Vector{T},Missing} 
    TEMP_8::Union{Vector{T},Missing} 
    DLNNDR_8::Union{Vector{T},Missing} 
    DLNTDR_8::Union{Vector{T},Missing}  #72

    DENS_9::Union{Vector{T},Missing} 
    TEMP_9::Union{Vector{T},Missing} 
    DLNNDR_9::Union{Vector{T},Missing} 
    DLNTDR_9::Union{Vector{T},Missing}  #76

    DENS_10::Union{Vector{T},Missing} 
    TEMP_10::Union{Vector{T},Missing} 
    DLNNDR_10::Union{Vector{T},Missing} 
    DLNTDR_10::Union{Vector{T},Missing}  #80

    CS::Union{Vector{T},Missing}  #81

    # function profile()
    #     return profile{Float64}()
    # end
    # As of right now, I don't believe there needs to be parameters, but the vectors
    # are probably the most of concern there. 
    function profile{T}(nr::Int, ns::Int) where (T<:Real)
        # new(NaN, NaN, nr, ns, 1, missing, fill(NaN, ns), fill(NaN, ns),
        new(NaN, NaN, missing, missing, 1, missing, fill(NaN, ns), fill(NaN, ns),
        fill(NaN, (nr, ns)), fill(NaN, (nr, ns)), fill(NaN, (nr, ns)), fill(NaN, (nr, ns)), 
        fill(NaN, (nr, ns)), fill(NaN, (nr, ns)), fill(NaN, nr), fill(NaN, nr), fill(NaN, nr), 
        fill(NaN, nr), fill(NaN, nr), fill(NaN, nr), fill(NaN, nr), fill(NaN, nr), fill(NaN, nr), 
        fill(NaN, nr), fill(NaN, nr), fill(NaN, nr), fill(NaN, nr), fill(NaN, nr), fill(NaN, nr), 
        fill(NaN, nr), fill(NaN, nr), fill(NaN, nr), fill(NaN, nr), fill(1.0E-7, nr), fill(NaN, nr), 
        fill(NaN, nr), missing, missing, NaN, missing, fill(NaN,(ns)), fill(NaN,(ns)), fill(NaN,(ns)),
        fill(NaN,(ns)), fill(NaN,(ns)), fill(NaN,(ns)), fill(NaN,(ns)), fill(NaN,(ns)), fill(NaN,(ns)),
        fill(NaN,(ns)), fill(NaN,(ns)), fill(NaN,(ns)), fill(NaN,(ns)), fill(NaN,(ns)), fill(NaN,(ns)),
        fill(NaN,(ns)), fill(NaN,(ns)), fill(NaN,(ns)), fill(NaN,(ns)), fill(NaN,(ns)), fill(NaN,(ns)),
        fill(NaN,(ns)), fill(NaN,(ns)), fill(NaN,(ns)), fill(NaN,(ns)), fill(NaN,(ns)), fill(NaN,(ns)),
        fill(NaN,(ns)), fill(NaN,(ns)), fill(NaN,(ns)), fill(NaN,(ns)), fill(NaN,(ns)), fill(NaN,(ns)),
        fill(NaN,(ns)), fill(NaN,(ns)), fill(NaN,(ns)), fill(NaN,(ns)), fill(NaN,(ns)), fill(NaN,(ns)),
        fill(NaN,(ns)), fill(NaN,(nr)))
    end
end

mutable struct InputTJLFEP{T<:Real}
    InputTJLF::InputTJLF
    Options::Options
    profile::profile

    # function InputTJLFEP()
    # end
end

