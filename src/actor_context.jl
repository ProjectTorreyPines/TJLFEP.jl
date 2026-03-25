# Packages imported from parent TJLF module via TJLFEP.jl
# (TurbulentTransport, InputTGLF already available)
# Import types defined within this TJLFEP module
using FUSE
import FUSE: ParametersActor, ParametersAllActors, SingleAbstractActor, logging_actor_init, @actor_parameters_struct
import SimulationParameters: Switch, Entry
import GACODE
# Note: TJLFEP types (InputTJLFEP, Options, profile) are defined in this module

#= ========= =#
#  ActorTGLF  #
#= ========= =#
Base.@kwdef mutable struct FUSEparameters__ActorTGLF{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    model::Switch{Symbol} = Switch{Symbol}([:TGLF, :TGLFNN, :TJLF, :TJLFEP], "-", "Implementation of TGLF"; default=:TGLFNN)
    sat_rule::Switch{Symbol} = Switch{Symbol}([:sat0, :sat0quench, :sat1, :sat1geo, :sat2, :sat3], "-", "Saturation rule"; default=:sat1)
    electromagnetic::Entry{Bool} = Entry{Bool}("-", "Electromagnetic or electrostatic"; default=true)
    tglfnn_model::Entry{String} = Entry{String}(
        "-",
        "Use a user specified TGLF-NN model stored in TGLFNN/models";
        check=x -> @assert x in TurbulentTransport.available_models() "ActorTGLF.tglfnn_model must be one of:\n  \"$(join(TurbulentTransport.available_models(),"\"\n  \""))\""
    )
    rho_transport::Entry{AbstractVector{T}} = Entry{AbstractVector{T}}("-", "rho_tor_norm values to compute tglf fluxes on"; default=0.25:0.1:0.85)
    warn_nn_train_bounds::Entry{Bool} = Entry{Bool}("-", "Raise warnings if querying cases that are certainly outside of the training range"; default=false)
    custom_input_files::Entry{Union{Vector{<:InputTGLF},Vector{<:InputTJLF}}} =
        Entry{Union{Vector{<:InputTGLF},Vector{<:InputTJLF}}}("-", "Sets up the input file that will be run with the custom input file as a mask")
    # lump_ions::Entry{Bool} = Entry{Bool}("-", "Lumps the fuel species (D,T) as well as the impurities together"; default=true)
    lump_ions::Entry{Bool} = Entry{Bool}("-", "Lumps the fuel species (D,T) as well as the impurities together"; default=false)
end

mutable struct ActorTGLF{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorTGLF{P}
    input_tglfs::Union{Vector{<:InputTGLF},Vector{<:InputTJLF}}
    flux_solutions::Vector{GACODE.FluxSolution{D}}
    # TJLFEP-specific fields (only used when par.model == :TJLFEP)
    tjlfep_options::Union{Options{Float64}, Nothing}
    tjlfep_profile::Union{profile{Float64}, Nothing}
end

"""
    ActorTGLF(dd::IMAS.dd, act::ParametersAllActors; kw...)

Evaluates the TGLF predicted turbulence
"""
function ActorTGLF(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorTGLF(dd, act.ActorTGLF; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorTGLF(dd::IMAS.dd, par::FUSEparameters__ActorTGLF; kw...)
    logging_actor_init(ActorTGLF)
    par = par(kw...)
    if par.model ∈ [:TGLF, :TGLFNN]
        input_tglfs = Vector{InputTGLF}(undef, length(par.rho_transport))
    #elseif par.model == :TJLF
    elseif par.model ∈ [:TJLF, :TJLFEP]
        input_tglfs = Vector{InputTJLF}(undef, length(par.rho_transport))
    end
    return ActorTGLF(dd, par, input_tglfs, GACODE.FluxSolution{D}[], nothing, nothing)
end

"""
    _step(actor::ActorTGLF)

Runs TGLF actor to evaluate the turbulence flux on a vector of gridpoints
"""
function _step(actor::ActorTGLF)
    par = actor.par
    dd = actor.dd

    input_tglfs, extraEP = InputTGLF(dd, par.rho_transport, par.sat_rule, par.electromagnetic, par.lump_ions)

    for k in eachindex(par.rho_transport)
        input_tglf = input_tglfs[k]
        if par.model ∈ [:TGLF, :TGLFNN]
            actor.input_tglfs[k] = input_tglf
        elseif par.model == :TJLF
            if !isassigned(actor.input_tglfs, k) # this is done to keep memory of the widths
                nky = TJLF.get_ky_spectrum_size(input_tglf.NKY, input_tglf.KYGRID_MODEL)
                actor.input_tglfs[k] = InputTJLF{Float64}(input_tglf.NS, nky)
                actor.input_tglfs[k].WIDTH_SPECTRUM .= 1.65
                actor.input_tglfs[k].FIND_WIDTH = true # first case should find the widths
            end
            update_input_tjlf!(actor.input_tglfs[k], input_tglf)

        elseif par.model == :TJLFEP
            # TJLFEP operates on ALL radii at once, not per radius like TGLF
            # Only initialize once on first iteration
            if k == 1
                nr = extraEP["NR"]
                ns = extraEP["NS"]
                nmodes = 4
                jtscale_max = 1
                nn = 15
                nscan_in = nr  # TJLFEP scans all radial points
                widthin = false
                
                # Initialize TJLFEP structs once (they're global to the actor)
                actor.tjlfep_options = Options{Float64}(nscan_in, widthin, nn, nr, jtscale_max, nmodes)
                actor.tjlfep_profile = profile{Float64}(nr, ns)
                
                # Populate profile struct with full radial data from extraEP
                populate_tjlfep_profile!(actor.tjlfep_profile, extraEP, nr, ns)
            end
            
            if !isassigned(actor.input_tglfs, k)
                nky = 12
                actor.input_tglfs[k] = InputTJLF{Float64}(input_tglf.NS, nky)
                actor.input_tglfs[k].WIDTH_SPECTRUM .= 1.65
                actor.input_tglfs[k].FIND_WIDTH = true
            end
            
            # Update this specific radial point's TJLF input
            update_input_tjlfep!(actor.input_tglfs[k], input_tglf, actor.tjlfep_options)
        end

        # Overwrite TGLF / TJLF parameters with the custom parameters mask
        if !ismissing(par, :custom_input_files)
            for field_name in fieldnames(typeof(actor.input_tglfs[k]))
                if !ismissing(getproperty(par.custom_input_files[k], field_name))
                    setproperty!(actor.input_tglfs[k], field_name, getproperty(par.custom_input_files[k], field_name))
                end
            end
        end
    end

    if par.model == :TGLFNN
        actor.flux_solutions = TurbulentTransport.run_tglfnn(actor.input_tglfs; par.warn_nn_train_bounds, model_filename=model_filename(par))

    elseif par.model == :TGLF
        actor.flux_solutions = TurbulentTransport.run_tglf(actor.input_tglfs)
    

    elseif par.model == :TJLFEP
        # TJLFEP.runTHD expects Vector{InputTJLF}, Options, and profile
        # It runs once for all radii, not per radius
        converted_input_tglfs = map(convert_to_TJLFEP, actor.input_tglfs)
        TJLFEPoutput = TJLFEP.runTHD(converted_input_tglfs, actor.tjlfep_options, actor.tjlfep_profile, false)
        
        # Convert TJLFEP output to flux_solutions format
        # TODO: Extract flux data from TJLFEPoutput and populate actor.flux_solutions
        @warn "TJLFEP flux extraction not yet implemented. TJLFEPoutput structure: $(typeof(TJLFEPoutput))"
        
    elseif par.model == :TJLF
        QL_fluxes_out = TJLF.run_tjlf(actor.input_tglfs)
        actor.flux_solutions =
            [GACODE.FluxSolution{D}(TJLF.Qe(QL_flux_out), TJLF.Qi(QL_flux_out), TJLF.Γe(QL_flux_out), TJLF.Γi(QL_flux_out), TJLF.Πi(QL_flux_out)) for QL_flux_out in QL_fluxes_out]
    end

    return actor
end

"""
    _finalize(actor::ActorTGLF)

Writes results to dd.core_transport
"""
function _finalize(actor::ActorTGLF)
    dd = actor.dd
    par = actor.par
    cp1d = dd.core_profiles.profiles_1d[]
    eqt = dd.equilibrium.time_slice[]

    model = resize!(dd.core_transport.model, :anomalous; wipe=false)
    model.identifier.name = string(par.model) * " " * model_filename(par)
    m1d = resize!(model.profiles_1d)
    m1d.grid_flux.rho_tor_norm = par.rho_transport

    GACODE.flux_gacode_to_imas((:electron_energy_flux, :ion_energy_flux, :electron_particle_flux, :ion_particle_flux, :momentum_flux), actor.flux_solutions, m1d, eqt, cp1d)

    return actor
end

function model_filename(par::FUSEparameters__ActorTGLF)
    if par.model == :TGLFNN
        filename = par.tglfnn_model
    else
        filename = string(par.sat_rule) * "_" * (par.electromagnetic ? "em" : "es")
    end
    return filename
end

"""
    update_input_tjlf!(input_tglf::InputTGLF)

Modifies an InputTJLF from a InputTGLF
"""
function update_input_tjlf!(input_tjlf::InputTJLF, input_tglf::InputTGLF)
    input_tjlf.NWIDTH = 21

    for fieldname in fieldnames(typeof(input_tglf))
        if occursin(r"\d", String(fieldname)) || fieldname == :_Qgb # species parameter
            continue
        end
        setfield!(input_tjlf, fieldname, getfield(input_tglf, fieldname))
    end

    for i in 1:input_tglf.NS
        input_tjlf.ZS[i] = getfield(input_tglf, Symbol("ZS_", i))
        input_tjlf.AS[i] = getfield(input_tglf, Symbol("AS_", i))
        input_tjlf.MASS[i] = getfield(input_tglf, Symbol("MASS_", i))
        input_tjlf.RLNS[i] = getfield(input_tglf, Symbol("RLNS_", i))
        input_tjlf.RLTS[i] = getfield(input_tglf, Symbol("RLTS_", i))
        input_tjlf.TAUS[i] = getfield(input_tglf, Symbol("TAUS_", i))
        input_tjlf.VPAR[i] = getfield(input_tglf, Symbol("VPAR_", i))
        input_tjlf.VPAR_SHEAR[i] = getfield(input_tglf, Symbol("VPAR_SHEAR_", i))
    end

    # Defaults
    input_tjlf.KY = 0.3
    input_tjlf.ALPHA_E = 1.0
    input_tjlf.ALPHA_P = 1.0
    input_tjlf.XNU_FACTOR = 1.0
    input_tjlf.DEBYE_FACTOR = 1.0
    input_tjlf.RLNP_CUTOFF = 18.0
    input_tjlf.WIDTH = 1.65
    input_tjlf.WIDTH_MIN = 0.3
    input_tjlf.BETA_LOC = 0.0
    input_tjlf.KX0_LOC = 1.0
    input_tjlf.PARK = 1.0
    input_tjlf.GHAT = 1.0
    input_tjlf.GCHAT = 1.0
    input_tjlf.WD_ZERO = 0.1
    input_tjlf.LINSKER_FACTOR = 0.0
    input_tjlf.GRADB_FACTOR = 0.0
    input_tjlf.FILTER = 2.0
    input_tjlf.THETA_TRAPPED = 0.7
    input_tjlf.ETG_FACTOR = 1.25
    input_tjlf.DAMP_PSI = 0.0
    input_tjlf.DAMP_SIG = 0.0

    input_tjlf.FIND_EIGEN = true
    input_tjlf.NXGRID = 16

    input_tjlf.ADIABATIC_ELEC = false
    input_tjlf.VPAR_MODEL = 0
    input_tjlf.NEW_EIKONAL = true
    input_tjlf.USE_BISECTION = true
    input_tjlf.USE_INBOARD_DETRAPPED = false
    input_tjlf.IFLUX = true
    input_tjlf.IBRANCH = -1
    input_tjlf.KX0_LOC = 0.0

    # for now settings
    input_tjlf.ALPHA_ZF = -1  # smooth   

    # check converison
    TJLF.checkInput(input_tjlf)

    return input_tjlf
end

function Base.show(io::IO, ::MIME"text/plain", input::Union{InputTGLF,InputTJLF})
    for field_name in fieldnames(typeof(input))
        println(io, " $field_name = $(getfield(input,field_name))")
    end
end



"""
    update_input_tjlfep!(input_tjlf::InputTJLF, input_tglf::InputTGLF, options::Options)

Modifies an InputTJLF from a InputTGLF for TJLFEP usage
"""
function update_input_tjlfep!(input_tjlf::InputTJLF, input_tglf::InputTGLF, options::Options)
    input_tjlf.NWIDTH = 21

    for fieldname in fieldnames(typeof(input_tglf))
        if occursin(r"\d", String(fieldname)) || fieldname == :_Qgb # species parameter
            continue
        end
        setfield!(input_tjlf, fieldname, getfield(input_tglf, fieldname))
    end

    for i in 1:input_tglf.NS
        input_tjlf.ZS[i] = getfield(input_tglf, Symbol("ZS_", i))
        input_tjlf.AS[i] = getfield(input_tglf, Symbol("AS_", i))
        input_tjlf.MASS[i] = getfield(input_tglf, Symbol("MASS_", i))
        input_tjlf.RLNS[i] = getfield(input_tglf, Symbol("RLNS_", i))
        input_tjlf.RLTS[i] = getfield(input_tglf, Symbol("RLTS_", i))
        input_tjlf.TAUS[i] = getfield(input_tglf, Symbol("TAUS_", i))
        input_tjlf.VPAR[i] = getfield(input_tglf, Symbol("VPAR_", i))
        input_tjlf.VPAR_SHEAR[i] = getfield(input_tglf, Symbol("VPAR_SHEAR_", i))
    end

    # Defaults
    input_tjlf.KY = 0.3
    input_tjlf.ALPHA_E = 1.0
    input_tjlf.ALPHA_P = 1.0
    input_tjlf.XNU_FACTOR = 1.0
    input_tjlf.DEBYE_FACTOR = 1.0
    input_tjlf.RLNP_CUTOFF = 18.0
    input_tjlf.WIDTH = 1.65
    input_tjlf.WIDTH_MIN = 0.3
    input_tjlf.BETA_LOC = 0.0
    input_tjlf.KX0_LOC = 1.0
    input_tjlf.PARK = 1.0
    input_tjlf.GHAT = 1.0
    input_tjlf.GCHAT = 1.0
    input_tjlf.WD_ZERO = 0.1
    input_tjlf.LINSKER_FACTOR = 0.0
    input_tjlf.GRADB_FACTOR = 0.0
    input_tjlf.FILTER = 2.0
    input_tjlf.THETA_TRAPPED = 0.7
    input_tjlf.ETG_FACTOR = 1.25
    input_tjlf.DAMP_PSI = 0.0
    input_tjlf.DAMP_SIG = 0.0

    input_tjlf.FIND_EIGEN = true
    input_tjlf.NXGRID = 16

    input_tjlf.ADIABATIC_ELEC = false
    input_tjlf.VPAR_MODEL = 0
    input_tjlf.NEW_EIKONAL = true
    input_tjlf.USE_BISECTION = true
    input_tjlf.USE_INBOARD_DETRAPPED = false
    input_tjlf.IFLUX = true
    input_tjlf.IBRANCH = -1
    input_tjlf.KX0_LOC = 0.0

    # for now settings
    input_tjlf.ALPHA_ZF = -1  # smooth   

    # TJLF-EP specific settings
    input_tjlf.SAT_RULE = 0
    
    # Set Options parameters (these are used by TJLFEP.mainsub)
    options.MODE_IN = 2
    options.PROCESS_IN = 5
    options.THRESHOLD_FLAG = 0
    options.N_BASIS = 2
    options.SCAN_METHOD = 1
    options.REJECT_I_PINCH_FLAG = 0
    options.REJECT_E_PINCH_FLAG = 0
    options.REJECT_TH_PINCH_FLAG = 1
    options.REJECT_EP_PINCH_FLAG = 0
    options.REJECT_TEARING_FLAG = 1
    options.ROTATIONAL_SUPPRESSION_FLAG = 1
    options.QL_RATIO_THRESH = 0.001
    options.THETA_SQ_THRESH = 100.0
    options.Q_SCALE = 1.0
    options.WRITE_WAVEFUNCTION = 1
    options.KY_MODEL = 2
    options.IRS = 2
    options.FACTOR_IN_PROFILE = false
    options.FACTOR_IN = 1.0
    options.WIDTH_IN_FLAG = false
    options.WIDTH_MIN = 1.0
    options.WIDTH_MAX = 2.0
    options.INPUT_PROFILE_METHOD = 2
    options.REAL_FREQ = 1

    return input_tjlf
end

"""
    populate_tjlfep_profile!(prof::profile, extraEP::Dict, nr::Int, ns::Int)

Populates TJLFEP profile struct from extraEP dictionary containing full radial data.
Note: This only sets the multi-radial EP-specific fields. Geometry and other fields
should be populated from the InputTGLF data via readMTGLF or similar methods.
"""
function populate_tjlfep_profile!(prof::profile, extraEP::Dict, input::InputTGLFs, nr::Int, ns::Int)
    
    # Copy full radial arrays from extraEP to profile struct
    println("prof ns is $ns and nr is $nr")
    prof.NS = ns
    prof.NR = nr
    prof.RMIN = extraEP["RMIN"]
    prof.omegaGAM = extraEP["omegaGAM"]
    prof.OMEGA_TAE = extraEP["OMEGA_TAE"]
    prof.RHO_STAR = extraEP["RHO_STAR"]
    prof.gammaE = extraEP["gammaE"]
    prof.gammap = extraEP["gammap"]
    prof.IRS = 2  # Standard starting radius index

    # Geometry and physics fields
    prof.SIGN_BT = extraEP["SIGN_BT"]
    prof.SIGN_IT = extraEP["SIGN_IT"]
    prof.GEOMETRY_FLAG = 1
    prof.ROTATION_FLAG = 0  # VPAR set to 0.0 in TJLF_map when ROTATION_FLAG == 0
    prof.RMAJ    = extraEP["RMAJ"]
    prof.SHIFT   = extraEP["SHIFT"]
    prof.Q       = extraEP["Q"]
    prof.SHEAR   = extraEP["SHEAR"]
    prof.Q_PRIME = extraEP["Q_PRIME"]
    prof.P_PRIME = extraEP["P_PRIME"]
    prof.KAPPA   = extraEP["KAPPA"]
    prof.S_KAPPA = extraEP["S_KAPPA"]
    prof.DELTA   = extraEP["DELTA"]
    prof.S_DELTA = extraEP["S_DELTA"]
    prof.ZETA    = extraEP["ZETA"]
    prof.S_ZETA  = extraEP["S_ZETA"]
    prof.BETAE   = extraEP["BETAE"]
    prof.ZEFF    = extraEP["ZEFF"]
    
    # Extract charge numbers from extraEP, which contains all 4 species (e, D, T, EP)
    # extraEP["ZS"] should have length >= 4

    prof.ZS = extraEP["ZS"]
    println("prof.ZS: ", prof.ZS)
    
    prof.MASS = extraEP["MASS"]
    println("prof.MASS: ", prof.MASS)

    prof.N_ION = extraEP["N_ION"]
    
    # Species data - populate matrices with full radial profiles
    # Profile struct uses [nr, ns] indexing
    # NOTE: AS and TAUS need to be NORMALIZED (ni/ne and Ti/Te)
    # The raw values are in extraEP, so we need to normalize them
    
    ne_full = extraEP["DENS_1"]  # Electron density
    Te_full = extraEP["TEMP_1"]  # Electron temperature
    
    for s in 1:ns
        dens_key = "DENS_$s"
        temp_key = "TEMP_$s"
        dlnndr_key = "DLNNDR_$s"
        dlntdr_key = "DLNTDR_$s"
        
        # All species data (including EP) come from extraEP
        if haskey(extraEP, dens_key)
            for ir in 1:nr
                if s == 1
                    # Electrons: AS=1, TAUS=1 by definition
                    prof.AS[ir, s] = 1.0
                    prof.TAUS[ir, s] = 1.0
                else
                    # Ions and EP: normalize to electrons
                    prof.AS[ir, s] = extraEP[dens_key][ir] / ne_full[ir]
                    prof.TAUS[ir, s] = extraEP[temp_key][ir] / Te_full[ir]
                end
                prof.RLNS[ir, s] = extraEP[dlnndr_key][ir]  # Density gradient
                prof.RLTS[ir, s] = extraEP[dlntdr_key][ir]  # Temperature gradient
            end
        end
    end
    
    return prof
end
function convert_to_TJLFEP(input::TJLF.InputTJLF{Float64})
    return TJLFEP.InputTJLF{Float64}(
        input.UNITS,  # Adjust as needed
        input.USE_BPER,
        input.USE_BPAR,
        input.USE_MHD_RULE,
        input.USE_BISECTION,
        input.USE_INBOARD_DETRAPPED,
        input.USE_AVE_ION_GRID,
        input.NEW_EIKONAL,
        input.FIND_WIDTH,
        input.IFLUX,
        input.ADIABATIC_ELEC,
        input.SAT_RULE,
        input.NS,
        input.NMODES,
        input.NWIDTH,
        input.NBASIS_MAX,
        input.NBASIS_MIN,
        input.NXGRID,
        input.NKY,
        input.KYGRID_MODEL,
        input.XNU_MODEL,
        input.VPAR_MODEL,
        input.IBRANCH,
        input.ZS,
        input.MASS,
        input.RLNS,
        input.RLTS,
        input.TAUS,
        input.AS,
        input.VPAR,
        input.VPAR_SHEAR,
        input.WIDTH_SPECTRUM,
        input.KY_SPECTRUM,
        input.EIGEN_SPECTRUM,
        input.FIND_EIGEN,
        input.SIGN_BT,
        input.SIGN_IT,
        input.KY,
        input.VEXB_SHEAR,
        input.BETAE,
        input.XNUE,
        input.ZEFF,
        input.DEBYE,
        input.ALPHA_MACH,
        input.ALPHA_E,
        input.ALPHA_P,
        input.ALPHA_QUENCH,
        input.ALPHA_ZF,
        input.XNU_FACTOR,
        input.DEBYE_FACTOR,
        input.ETG_FACTOR,
        input.RLNP_CUTOFF,
        input.WIDTH,
        input.WIDTH_MIN,
        input.RMIN_LOC,
        input.RMAJ_LOC,
        input.ZMAJ_LOC,
        input.DRMINDX_LOC,
        input.DRMAJDX_LOC,
        input.DZMAJDX_LOC,
        input.Q_LOC,
        input.KAPPA_LOC,
        input.S_KAPPA_LOC,
        input.DELTA_LOC,
        input.S_DELTA_LOC,
        input.ZETA_LOC,
        input.S_ZETA_LOC,
        input.P_PRIME_LOC,
        input.Q_PRIME_LOC,
        input.BETA_LOC,
        input.KX0_LOC,
        input.DAMP_PSI,
        input.DAMP_SIG,
        input.WDIA_TRAPPED,
        input.PARK,
        input.GHAT,
        input.GCHAT,
        input.WD_ZERO,
        input.LINSKER_FACTOR,
        input.GRADB_FACTOR,
        input.FILTER,
        input.THETA_TRAPPED,
        input.SMALL,
        input.USE_TRANSPORT_MODEL
    )
end