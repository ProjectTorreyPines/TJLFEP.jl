# Packages imported from parent TJLF module via TJLFEP.jl
# (IMAS, GACODE, TurbulentTransport, InputTGLF are already available)
#  Base.@kwdef mutable struct EPconst (rho::AbstractVector{Float64}){T<:Real}
#      #::Union{Vector{T}, Missing}
#     DENS_1::Union{Vector{T}, Missing} = fill(missing, nr)
#     TEMP_1::Union{Vector{T}, Missing} = fill(missing, nr)
#     DLNNDR_1::Union{Vector{T}, Missing} = fill(missing, nr)
#     DLNTDR_1::Union{Vector{T}, Missing} = fill(missing, nr)
#     DENS_2::Union{Vector{T}, Missing} = fill(missing, nr)#5
#     TEMP_2::Union{Vector{T}, Missing} = fill(missing, nr)
#     DLNNDR_2::Union{Vector{T}, Missing} = fill(missing, nr)
#     DLNTDR_2::Union{Vector{T}, Missing} = fill(missing, nr)
#     DENS_3::Union{Vector{T}, Missing} = fill(missing, nr)
#     TEMP_3::Union{Vector{T}, Missing} = fill(missing, nr) #10
#     DLNNDR_3::Union{Vector{T}, Missing} = fill(missing, nr)
#     DLNTDR_3::Union{Vector{T}, Missing} = fill(missing, nr)
#     DENS_4::Union{Vector{T}, Missing} = fill(missing, nr)
#     TEMP_4::Union{Vector{T}, Missing} = fill(missing, nr)
#     DLNNDR_4::Union{Vector{T}, Missing} = fill(missing, nr) #15
#     DLNTDR_4::Union{Vector{T}, Missing} = fill(missing, nr)
#     DENS_5::Union{Vector{T}, Missing} = fill(missing, nr)
#     TEMP_5::Union{Vector{T}, Missing} = fill(missing, nr)
#     DLNNDR_5::Union{Vector{T}, Missing} = fill(missing, nr)
#     DLNTDR_5::Union{Vector{T}, Missing} = fill(missing, nr)#20
#     DENS_6::Union{Vector{T}, Missing} = fill(missing, nr)
#     TEMP_6::Union{Vector{T}, Missing} = fill(missing, nr)
#     DLNNDR_6::Union{Vector{T}, Missing} = fill(missing, nr)
#     gammaE::Union{Vector{T}, Missing} = fill(missing, nr) #25
#     gammap::Union{Vector{T}, Missing} = fill(missing, nr)
#     omegaGAM::Union{Vector{T}, Missing} = fill(missing, nr)
#     OMEGA_TAE::Union{Vector{T}, Missing} = fill(missing, nr)
#     RHO_STAR::Union{Vector{T}, Missing} = fill(missing, nr) #30
#     CS::Union{Vector{T}, Missing} = fill(missing, nr)
#     RMIN::Union{Vector{T}, Missing} = fill(missing, nr) 
#     NR::Union{Int, Missing} = missing #32
#     NS::Union{Int, Missing} = missing
#     end
 

# mutable struct EPconst{T<:Real}
#     DENS_1::Union{Vector{T}, Missing}
#     TEMP_1::Union{Vector{T}, Missing}
#     DLNNDR_1::Union{Vector{T}, Missing}
#     DLNTDR_1::Union{Vector{T}, Missing}
#     DENS_2::Union{Vector{T}, Missing}
#     TEMP_2::Union{Vector{T}, Missing}
#     DLNNDR_2::Union{Vector{T}, Missing}
#     DLNTDR_2::Union{Vector{T}, Missing}
#     DENS_3::Union{Vector{T}, Missing}
#     TEMP_3::Union{Vector{T}, Missing}
#     DLNNDR_3::Union{Vector{T}, Missing}
#     DLNTDR_3::Union{Vector{T}, Missing}
#     DENS_4::Union{Vector{T}, Missing}
#     TEMP_4::Union{Vector{T}, Missing}
#     DLNNDR_4::Union{Vector{T}, Missing}
#     DLNTDR_4::Union{Vector{T}, Missing}
#     DENS_5::Union{Vector{T}, Missing}
#     TEMP_5::Union{Vector{T}, Missing}
#     DLNNDR_5::Union{Vector{T}, Missing}
#     DLNTDR_5::Union{Vector{T}, Missing}
#     DENS_6::Union{Vector{T}, Missing}
#     TEMP_6::Union{Vector{T}, Missing}
#     DLNNDR_6::Union{Vector{T}, Missing}
#     DLNTDR_6::Union{Vector{T}, Missing}
#     gammaE::Union{Vector{T}, Missing}
#     gammap::Union{Vector{T}, Missing}
#     omegaGAM::Union{Vector{T}, Missing}
#     OMEGA_TAE::Union{Vector{T}, Missing}
#     RHO_STAR::Union{Vector{T}, Missing}
#     CS::Union{Vector{T}, Missing}
#     RMIN::Union{Vector{T}, Missing}
#     NR::Union{Int, Missing}
#     NS::Union{Int, Missing}

#     function EPconst{T}(nr::Int) where {T<:Real}
#         new("",
#         fill(NaN, nr),  # NaN will be automatically treated as type T
#         fill(NaN, nr),
#         fill(NaN, nr),
#         fill(NaN, nr),
#         fill(NaN, nr),
#         fill(NaN, nr),
#         fill(NaN, nr),
#         fill(NaN, nr),
#         fill(NaN, nr),
#         fill(NaN, nr),
#         fill(NaN, nr),
#         fill(NaN, nr),
#         fill(NaN, nr),
#         fill(NaN, nr),
#         fill(NaN, nr),
#         fill(NaN, nr),
#         fill(NaN, nr),
#         fill(NaN, nr),
#         fill(NaN, nr),
#         fill(NaN, nr),
#         fill(NaN, nr),
#         fill(NaN, nr),
#         fill(NaN, nr),
#         fill(NaN, nr),
#         fill(NaN, nr),
#         fill(NaN, nr),
#         fill(NaN, nr),
#         fill(NaN, nr),
#         nr,missing)
            
#     end
# end

# NOTE: InputTGLF is imported from TJLF via TurbulentTransport
# The struct definition below was causing type conflicts, so it's been removed.
# The functions below extend TJLF.InputTGLF with additional constructors.
mutable struct InputTGLFs
    tglfs::Vector  # Allow any InputTGLF type (TJLF.InputTGLF or compatible)
end

function Base.setproperty!(inputTGLFs::InputTGLFs, field::Symbol, value::AbstractVector{<:Any})
    @assert length(value) == length(getfield(inputTGLFs, :tglfs))
    for (k, inputTGLF) in enumerate(getfield(inputTGLFs, :tglfs))
        setproperty!(inputTGLF, field, value[k])
    end
    return value
end

function Base.setproperty!(inputTGLFs::InputTGLFs, field::Symbol, value::Any)
    for (k, inputTGLF) in enumerate(getfield(inputTGLFs, :tglfs))
        setproperty!(inputTGLF, field, value)
    end
    return value
end

function Base.getproperty(inputTGLFs::InputTGLFs, field::Symbol)
    if field === :tglfs
        return getfield(inputTGLFs, :tglfs)
    end
    data = typeof(getfield(getfield(inputTGLFs, :tglfs)[1], field))[]
    for inputTGLF in getfield(inputTGLFs, :tglfs)
        push!(data, getproperty(inputTGLF, field))
    end
    return data
end

function Base.getindex(inputTGLFs::InputTGLFs, index::Int)
    return getfield(inputTGLFs, :tglfs)[index]
end

function diff(A::InputTGLF, B::InputTGLF)
    differences = Symbol[]
    for field in fieldnames(typeof(A))
        if getfield(A, field) === getfield(B, field)
            # pass
        else
            push!(differences, field)
        end
    end
    return differences
end

Base.@kwdef mutable struct InputCGYRO
    EQUILIBRIUM_MODEL::Union{Int,Missing} = missing
    RMIN::Union{Float64,Missing} = missing
    RMAJ::Union{Float64,Missing} = missing
    SHIFT::Union{Float64,Missing} = missing
    KAPPA::Union{Float64,Missing} = missing
    S_KAPPA::Union{Float64,Missing} = missing
    DELTA::Union{Float64,Missing} = missing
    S_DELTA::Union{Float64,Missing} = missing
    ZETA::Union{Float64,Missing} = missing
    S_ZETA::Union{Float64,Missing} = missing
    ZMAG::Union{Float64,Missing} = missing
    DZMAG::Union{Float64,Missing} = missing
    Q::Union{Float64,Missing} = missing
    S::Union{Float64,Missing} = missing
    BTCCW::Union{Float64,Missing} = missing
    IPCCW::Union{Float64,Missing} = missing
    UDSYMMETRY_FLAG::Union{Float64,Missing} = missing

    SHAPE_COS0::Union{Float64,Missing} = missing
    SHAPE_S_COS0::Union{Float64,Missing} = missing
    SHAPE_COS1::Union{Float64,Missing} = missing
    SHAPE_S_COS1::Union{Float64,Missing} = missing
    SHAPE_COS2::Union{Float64,Missing} = missing
    SHAPE_S_COS2::Union{Float64,Missing} = missing
    SHAPE_COS3::Union{Float64,Missing} = missing
    SHAPE_S_COS3::Union{Float64,Missing} = missing
    SHAPE_SIN3::Union{Float64,Missing} = missing
    SHAPE_S_SIN3::Union{Float64,Missing} = missing

    PROFILE_MODEL::Union{Int,Missing} = missing
    QUASINEUTRAL_FLAG::Union{Int,Missing} = missing
    NONLINEAR_FLAG::Union{Int,Missing} = missing
    ZF_TEST_MODE::Union{Int,Missing} = missing
    SILENT_FLAG::Union{Int,Missing} = missing
    AMP::Union{Float64,Missing} = missing
    AMP0::Union{Float64,Missing} = missing

    N_FIELD::Union{Int,Missing} = missing
    BETAE_UNIT::Union{Float64,Missing} = missing
    BETAE_UNIT_SCALE::Union{Float64,Missing} = missing
    BETA_STAR_SCALE::Union{Float64,Missing} = missing
    LAMBDA_DEBYE::Union{Float64,Missing} = missing
    LAMBDA_DEBYE_SCALE::Union{Float64,Missing} = missing

    N_RADIAL::Union{Int,Missing} = missing
    BOX_SIZE::Union{Int,Missing} = missing
    N_TOROIDAL::Union{Int,Missing} = missing
    KY::Union{Float64,Missing} = missing
    N_THETA::Union{Int,Missing} = missing
    N_XI::Union{Int,Missing} = missing
    N_ENERGY::Union{Int,Missing} = missing
    E_MAX::Union{Float64,Missing} = missing

    UP_RADIAL::Union{Float64,Missing} = missing
    UP_THETA::Union{Float64,Missing} = missing
    UP_ALPHA::Union{Float64,Missing} = missing
    NUP_RADIAL::Union{Int,Missing} = missing
    NUP_THETA::Union{Int,Missing} = missing
    NUP_ALPHA::Union{Int,Missing} = missing
    VELOCITY_ORDER::Union{Int,Missing} = missing
    UPWIND_SINGLE_FLAG::Union{Int,Missing} = missing

    DELTA_T_METHOD::Union{Int,Missing} = missing
    DELTA_T::Union{Float64,Missing} = missing
    ERROR_TOL::Union{Float64,Missing} = missing
    MAX_TIME::Union{Float64,Missing} = missing
    FREQ_TOL::Union{Float64,Missing} = missing
    PRINT_STEP::Union{Float64,Missing} = missing
    RESTART_STEP::Union{Float64,Missing} = missing

    N_SPECIES::Union{Int,Missing} = missing
    Z_1::Union{Float64,Missing} = missing
    MASS_1::Union{Float64,Missing} = missing
    DENS_1::Union{Float64,Missing} = missing
    TEMP_1::Union{Float64,Missing} = missing
    DLNNDR_1::Union{Float64,Missing} = missing
    DLNTDR_1::Union{Float64,Missing} = missing

    Z_2::Union{Float64,Missing} = missing
    MASS_2::Union{Float64,Missing} = missing
    DENS_2::Union{Float64,Missing} = missing
    TEMP_2::Union{Float64,Missing} = missing
    DLNNDR_2::Union{Float64,Missing} = missing
    DLNTDR_2::Union{Float64,Missing} = missing

    Z_3::Union{Float64,Missing} = missing
    MASS_3::Union{Float64,Missing} = missing
    DENS_3::Union{Float64,Missing} = missing
    TEMP_3::Union{Float64,Missing} = missing
    DLNNDR_3::Union{Float64,Missing} = missing
    DLNTDR_3::Union{Float64,Missing} = missing

    Z_4::Union{Float64,Missing} = missing
    MASS_4::Union{Float64,Missing} = missing
    DENS_4::Union{Float64,Missing} = missing
    TEMP_4::Union{Float64,Missing} = missing
    DLNNDR_4::Union{Float64,Missing} = missing
    DLNTDR_4::Union{Float64,Missing} = missing

    Z_5::Union{Float64,Missing} = missing
    MASS_5::Union{Float64,Missing} = missing
    DENS_5::Union{Float64,Missing} = missing
    TEMP_5::Union{Float64,Missing} = missing
    DLNNDR_5::Union{Float64,Missing} = missing
    DLNTDR_5::Union{Float64,Missing} = missing

    Z_6::Union{Float64,Missing} = missing
    MASS_6::Union{Float64,Missing} = missing
    DENS_6::Union{Float64,Missing} = missing
    TEMP_6::Union{Float64,Missing} = missing
    DLNNDR_6::Union{Float64,Missing} = missing
    DLNTDR_6::Union{Float64,Missing} = missing

    Z_7::Union{Float64,Missing} = missing
    MASS_7::Union{Float64,Missing} = missing
    DENS_7::Union{Float64,Missing} = missing
    TEMP_7::Union{Float64,Missing} = missing
    DLNNDR_7::Union{Float64,Missing} = missing
    DLNTDR_7::Union{Float64,Missing} = missing

    Z_8::Union{Float64,Missing} = missing
    MASS_8::Union{Float64,Missing} = missing
    DENS_8::Union{Float64,Missing} = missing
    TEMP_8::Union{Float64,Missing} = missing
    DLNNDR_8::Union{Float64,Missing} = missing
    DLNTDR_8::Union{Float64,Missing} = missing

    Z_9::Union{Float64,Missing} = missing
    MASS_9::Union{Float64,Missing} = missing
    DENS_9::Union{Float64,Missing} = missing
    TEMP_9::Union{Float64,Missing} = missing
    DLNNDR_9::Union{Float64,Missing} = missing
    DLNTDR_9::Union{Float64,Missing} = missing

    Z_10::Union{Float64,Missing} = missing
    MASS_10::Union{Float64,Missing} = missing
    DENS_10::Union{Float64,Missing} = missing
    TEMP_10::Union{Float64,Missing} = missing
    DLNNDR_10::Union{Float64,Missing} = missing
    DLNTDR_10::Union{Float64,Missing} = missing

    Z_11::Union{Float64,Missing} = missing
    MASS_11::Union{Float64,Missing} = missing
    DENS_11::Union{Float64,Missing} = missing
    TEMP_11::Union{Float64,Missing} = missing
    DLNNDR_11::Union{Float64,Missing} = missing
    DLNTDR_11::Union{Float64,Missing} = missing

    NU_EE::Union{Float64,Missing} = missing
    COLLISION_MODEL::Union{Int64,Missing} = missing
    COLLISION_FIELD_MODEL::Union{Int64,Missing} = missing
    COLLISION_MOM_RESTORE::Union{Int64,Missing} = missing
    COLLISION_ENE_RESTORE::Union{Int64,Missing} = missing
    COLLISION_ENE_DIFFUSION::Union{Int64,Missing} = missing
    COLLISION_KPERP::Union{Int64,Missing} = missing
    COLLISION_PRECISION_MODE::Union{Int64,Missing} = missing
    GPU_BIGMEM_FLAG::Union{Int64,Missing} = missing

    ROTATION_MODEL::Union{Int64,Missing} = missing
    GAMMA_E::Union{Float64,Missing} = missing
    GAMMA_P::Union{Float64,Missing} = missing
    MACH::Union{Float64,Missing} = missing
    GAMMA_E_SCALE::Union{Float64,Missing} = missing
    GAMMA_P_SCALE::Union{Float64,Missing} = missing
    MACH_SCALE::Union{Float64,Missing} = missing

    N_GLOBAL::Union{Int64,Missing} = missing
    NU_GLOBAL::Union{Float64,Missing} = missing
    FIELD_PRINT_FLAG::Union{Int64,Missing} = missing
    MOMENT_PRINT_FLAG::Union{Int64,Missing} = missing
    GFLUX_PRINT_FLAG::Union{Int64,Missing} = missing
    H_PRINT_FLAG::Union{Int64,Missing} = missing

end

Base.@kwdef mutable struct InputQLGYRO
    N_PARALLEL::Union{Int,Missing} = missing
    N_RUNS::Union{Int,Missing} = missing
    GAMMA_E::Union{Float64,Missing} = missing
    CODE::Union{Int,Missing} = missing
    NKY::Union{Int,Missing} = missing
    KYGRID_MODEL::Union{Int,Missing} = missing
    KY::Union{Float64,Missing} = missing
    SAT_RULE::Union{Int,Missing} = missing
end

function scan(input_tglf::InputTGLF; kw...)
    # Convert keyword arguments to a dictionary for easier handling
    kw_dict = Dict(kw)

    # Base case: if no keywords are left, return the current input_tglf
    if isempty(kw_dict)
        return [deepcopy(input_tglf)]
    end

    # Extract one keyword and its values
    key, values = pop!(kw_dict)
    intglfs = InputTGLF[]

    # Iterate over the values for the current keyword
    for v in values
        tmp = deepcopy(input_tglf)
        setproperty!(tmp, key, v)

        # Recursively call scan for the rest of the keywords
        intglfs = vcat(intglfs, scan(tmp; kw_dict...))
    end

    return intglfs
end

function Base.show(io::IO, ::MIME"text/plain", input_tglf::InputTGLF)
    for fname in sort(collect(fieldnames(typeof(input_tglf))))
        value = getfield(input_tglf, fname)
        if value !== missing && (!isdigit(string(fname)[end]) || (isdigit(string(fname)[end]) && parse(Int, split(string(fname), "_")[end]) <= input_tglf.NS))
            println(io, "$fname = $(value)")
        end
    end
end

"""
    InputTGLFEP(dd::IMAS.dd, rho::AbstractVector{Float64}, sat::Symbol=:sat0, electromagnetic::Bool=false, lump_ions::Bool=true)

Evaluate TGLF input parameters at given radii
"""
function InputTGLFEP(dd::IMAS.dd, rho::AbstractVector{Float64}, sat::Symbol=:sat0, electromagnetic::Bool=false, lump_ions::Bool=false)
    eqt = dd.equilibrium.time_slice[]
    cp1d = dd.core_profiles.profiles_1d[]
    gridpoint_cp = [argmin(abs.(cp1d.grid.rho_tor_norm .- ρ)) for ρ in rho]
    return InputTGLF_EP(eqt, cp1d, gridpoint_cp, sat, electromagnetic, lump_ions)
end

"""
    InputTGLFEP(dd::IMAS.dd, gridpoint_cp::AbstractVector{Int}, sat::Symbol=:sat0, electromagnetic::Bool=false, lump_ions::Bool=true)

Evaluate TGLF input parameters at given core profiles grid indexes
"""
function InputTGLFEP(dd::IMAS.dd, gridpoint_cp::AbstractVector{Int}, sat::Symbol=:sat0, electromagnetic::Bool=false, lump_ions::Bool=false)
    eqt = dd.equilibrium.time_slice[]
    cp1d = dd.core_profiles.profiles_1d[]
    return InputTGLF_EP(eqt, cp1d, gridpoint_cp, sat, electromagnetic, lump_ions)
end

function InputTGLF_EP(
    eqt::IMAS.equilibrium__time_slice,
    cp1d::IMAS.core_profiles__profiles_1d,
    gridpoint_cp::AbstractVector{Int},
    sat::Symbol,
    electromagnetic::Bool,
    lump_ions::Bool)

    e = IMAS.cgs.e # statcoul
    k = IMAS.cgs.k # erg/eV
    mp = IMAS.cgs.mp # g
    me = IMAS.cgs.me # g
    md = IMAS.cgs.md # g
    m_to_cm = IMAS.cgs.m_to_cm
    m³_to_cm³ = IMAS.cgs.m³_to_cm³
    T_to_Gauss = IMAS.cgs.T_to_Gauss

    eq1d = eqt.profiles_1d

    if lump_ions
        ions = IMAS.lump_ions_as_bulk_and_impurity(cp1d)
    else
        ions = cp1d.ion
    end

    Rmaj = IMAS.interp1d(eq1d.rho_tor_norm, m_to_cm * 0.5 * (eq1d.r_outboard .+ eq1d.r_inboard)).(cp1d.grid.rho_tor_norm)

    rmin = GACODE.r_min_core_profiles(eqt.profiles_1d, cp1d.grid.rho_tor_norm)

    q_profile = IMAS.interp1d(eq1d.rho_tor_norm, eq1d.q).(cp1d.grid.rho_tor_norm)

    if !ismissing(eq1d, :elongation)
        kappa = IMAS.interp1d(eq1d.rho_tor_norm, eq1d.elongation).(cp1d.grid.rho_tor_norm)
    else
        kappa = zero(cp1d.grid.rho_tor_norm)
    end

    if !ismissing(eq1d, :triangularity_lower) && !ismissing(eq1d, :triangularity_upper)
        delta = IMAS.interp1d(eq1d.rho_tor_norm, 0.5 * (eq1d.triangularity_lower + eq1d.triangularity_upper)).(cp1d.grid.rho_tor_norm)
    else
        delta = zero(cp1d.grid.rho_tor_norm)
    end

    if !ismissing(eq1d, :squareness_lower_inner) && !ismissing(eq1d, :squareness_lower_outer) && !ismissing(eq1d, :squareness_upper_inner) &&
       !ismissing(eq1d, :squareness_upper_outer)
        tmp = 0.25 * (eq1d.squareness_lower_inner .+ eq1d.squareness_lower_outer .+ eq1d.squareness_upper_inner .+ eq1d.squareness_upper_outer)
        zeta = IMAS.interp1d(eq1d.rho_tor_norm, tmp).(cp1d.grid.rho_tor_norm)
    else
        zeta = zero(cp1d.grid.rho_tor_norm)
    end

    a = rmin[end]
    q = q_profile[gridpoint_cp]

    Bt = eqt.global_quantities.vacuum_toroidal_field.b0
    bunit = IMAS.interp1d(eq1d.rho_tor_norm, GACODE.bunit(eqt) .* T_to_Gauss).(cp1d.grid.rho_tor_norm)[gridpoint_cp]
    input_tglf = InputTGLFs([TurbulentTransport.InputTGLF() for k in eachindex(gridpoint_cp)])  # Use local InputTGLFs struct

    signb = sign(Bt)
    signq = sign.(q)
    input_tglf.SIGN_BT = signb
    input_tglf.SIGN_IT = signb .* signq

    input_tglf.NS = length(ions) + 1 + 3 # add 1 to include electrons, 3 to include EP
    #creating instancefor extraEP
    # nr = Int(length(gridpoint_cp))
    # println(gridpoint_cp)
    # println(length(gridpoint_cp))
    # println(nr)
    # println(typeof(nr))
    
    # extraEP = EPconst(nr)


    c_s = GACODE.c_s(cp1d)[gridpoint_cp]
    w0 = -1 * cp1d.rotation_frequency_tor_sonic
    w0p = IMAS.gradient(rmin, w0)
    gamma_p = -Rmaj[gridpoint_cp] .* w0p[gridpoint_cp]
    gamma_e = -rmin[gridpoint_cp] ./ q .* w0p[gridpoint_cp]
    mach = Rmaj[gridpoint_cp] .* w0[gridpoint_cp] ./ c_s
    input_tglf.VPAR_1 = -input_tglf.SIGN_IT .* mach
    input_tglf.VPAR_SHEAR_1 = -input_tglf.SIGN_IT .* (a ./ c_s) .* gamma_p
    input_tglf.VEXB_SHEAR = -gamma_e .* (a ./ c_s)

    # ===== EP Data: Calculate for ALL radii (TJLFEP needs full profiles) =====
    # Get FULL radial grid (not just selected points)
    nr_full = length(cp1d.grid.rho_tor_norm)
    ns_full = length(ions) + 1 + 3  # +1 for electrons, +3 for ep
    
    # Preallocate EP arrays for all radii and species (6 max species for compatibility)
    ep_species_data = Dict{String, Vector{Float64}}()
    
    # Full radial profiles for thermal electrons
    Te_full = cp1d.electrons.temperature
    dlntedr_full = -IMAS.calc_z(rmin, Te_full, :backward)
    
    ne_full = cp1d.electrons.density_thermal ./ m³_to_cm³
    dlnnedr_full = -IMAS.calc_z(rmin, ne_full, :backward)
    
    # Store full radial data for TJLFEP
    ep_species_data["DENS_1"] = ne_full
    ep_species_data["TEMP_1"] = Te_full
    ep_species_data["DLNNDR_1"] = dlnnedr_full
    ep_species_data["DLNTDR_1"] = dlntedr_full
    ep_species_data["ZS_1"] = fill(-1.0, nr_full)

    # Extract values at selected grid points for TGLF
    Te = Te_full[gridpoint_cp]
    dlntedr = dlntedr_full[gridpoint_cp]
    ne = ne_full[gridpoint_cp]
    dlnnedr = dlnnedr_full[gridpoint_cp]

    Ze = -1.0

    setproperty!(input_tglf, Symbol("ZS_1"), Ze)
    setproperty!(input_tglf, Symbol("MASS_1"), me / md)
    setproperty!(input_tglf, Symbol("TAUS_1"), Te ./ Te)
    setproperty!(input_tglf, Symbol("AS_1"), ne ./ ne)
    setproperty!(input_tglf, Symbol("VPAR_1"), input_tglf.VPAR_1)
    setproperty!(input_tglf, Symbol("VPAR_SHEAR_1"), input_tglf.VPAR_SHEAR_1)
    setproperty!(input_tglf, Symbol("RLNS_1"), a * dlnnedr)
    setproperty!(input_tglf, Symbol("RLTS_1"), a * dlntedr)

    # Full radial profiles for fast electrons
    ne_full = cp1d.electrons.density_fast ./ m³_to_cm³
    fast_zero = ne_full .== 0
    dlnnedr_full = ifelse.(fast_zero, 0.0, -IMAS.calc_z(rmin, max.(ne_full, eps(Float64)), :backward))

    Te_full = (cp1d.electrons.pressure_fast_parallel ./ 3 + 2 * cp1d.electrons.pressure_fast_perpendicular ./ 3) ./ max.(ne_full, eps(Float64))
    Te_full[fast_zero] .= 0.0
    dlntedr_full = ifelse.(fast_zero, 0.0, -IMAS.calc_z(rmin, max.(Te_full, eps(Float64)), :backward))
    dlntedr_full[fast_zero] .= 0.0
    
    # Store full radial data for TJLFEP
    ep_species_data["DENS_4"] = ne_full
    ep_species_data["TEMP_4"] = Te_full
    ep_species_data["DLNNDR_4"] = dlnnedr_full
    ep_species_data["DLNTDR_4"] = dlntedr_full
    ep_species_data["ZS_4"] = fill(-1.0, nr_full)  # fast electrons

    # Extract values at selected grid points for TGLF
    Te_fast = Te_full[gridpoint_cp]
    dlntedr_fast = dlntedr_full[gridpoint_cp]
    nef = ne_full[gridpoint_cp]
    dlnnedr_fast = dlnnedr_full[gridpoint_cp]

    println("ne vs nef:", ne, nef)

    Ze = -1.0
    setproperty!(input_tglf, Symbol("ZS_4"), Ze)
    setproperty!(input_tglf, Symbol("MASS_4"), me / md)
    setproperty!(input_tglf, Symbol("TAUS_4"), Te_fast ./ Te)
    setproperty!(input_tglf, Symbol("AS_4"), nef ./ ne)
    setproperty!(input_tglf, Symbol("VPAR_4"), input_tglf.VPAR_1)
    setproperty!(input_tglf, Symbol("VPAR_SHEAR_4"), input_tglf.VPAR_SHEAR_1)
    setproperty!(input_tglf, Symbol("RLNS_4"), a * dlnnedr_fast)
    setproperty!(input_tglf, Symbol("RLTS_4"), a * dlntedr_fast)

    println("Te 4 is ", size(Te_fast), Te_fast)
    println("dlntedr 4 is ", size(dlntedr_fast), dlntedr_fast)
    println("dlnnedr 4 is ", size(dlnnedr_fast), dlnnedr_fast)

    # Ion species (2 through NS) - full radial profiles
    for iion in eachindex(ions)
        species = iion + 1

        z_n = ions[iion].element[1].z_n
        println("============================================================")


        # Full radial profiles for this ion species
        Ti_full = ions[iion].temperature
        dlntidr_full = -IMAS.calc_z(rmin, Ti_full, :backward)
        
        ni_full = ions[iion].density_thermal ./ m³_to_cm³
        dlnnidr_full = -IMAS.calc_z(rmin, ni_full, :backward)
        
        # Store full radial data for TJLFEP
        ep_species_data["DENS_$species"] = ni_full
        ep_species_data["TEMP_$species"] = Ti_full
        ep_species_data["DLNNDR_$species"] = dlnnidr_full
        ep_species_data["DLNTDR_$species"] = dlntidr_full
        ep_species_data["ZS_$species"] = IMAS.avgZ(Float64(z_n), Ti_full)

        # println("ni_full $species:", ni_full)
        # println("Ti_full $species:", Ti_full)
        # println("dlntidr_full $species:", dlntidr_full)
        # println("dlnnidr_full $species:", dlnnidr_full)

        # Extract values at selected grid points for TGLF
        Ti = Ti_full[gridpoint_cp]
        dlntidr = dlntidr_full[gridpoint_cp]
        ni = ni_full[gridpoint_cp]
        dlnnidr = dlnnidr_full[gridpoint_cp]

        println("ni for $iion is ", size(ni), ni)

        Zi = IMAS.avgZ(Float64(ions[iion].element[1].z_n), Ti)

        # println("Ti for $iion is ", size(Ti), Ti)
        println("Zi for $iion is ", size(Zi), Zi)

        setproperty!(input_tglf, Symbol("ZS_$species"), Zi)
        setproperty!(input_tglf, Symbol("MASS_$species"), (ions[iion].element[1].a)[1] * mp / md)
        setproperty!(input_tglf, Symbol("TAUS_$species"), Ti ./ Te)
        setproperty!(input_tglf, Symbol("AS_$species"), ni ./ ne)
        setproperty!(input_tglf, Symbol("VPAR_$species"), input_tglf.VPAR_1)
        setproperty!(input_tglf, Symbol("VPAR_SHEAR_$species"), input_tglf.VPAR_SHEAR_1)
        setproperty!(input_tglf, Symbol("RLNS_$species"), a * dlnnidr)
        setproperty!(input_tglf, Symbol("RLTS_$species"), a * dlntidr)


        fast = iion + 1 + 3

        # Full radial profiles for this ion species
        ni_full = ions[iion].density_fast ./ m³_to_cm³
        fast_zero = ni_full .== 0
        dlnnidr_full = ifelse.(fast_zero, 0.0, -IMAS.calc_z(rmin, max.(ni_full, eps(Float64)), :backward))

        Ti_full = (ions[iion].pressure_fast_parallel ./ 3 + 2 * ions[iion].pressure_fast_perpendicular ./ 3) ./ max.(ni_full, eps(Float64))
        Ti_full[fast_zero] .= 0.0
        dlntidr_full = ifelse.(fast_zero, 0.0, -IMAS.calc_z(rmin, max.(Ti_full, eps(Float64)), :backward))
        dlntidr_full[fast_zero] .= 0.0
        # print("Ti for fast $iion is ", size(Ti_full), Ti_full)
        
        # Store full radial data for TJLFEP
        ep_species_data["DENS_$fast"] = ni_full
        ep_species_data["TEMP_$fast"] = Ti_full
        ep_species_data["DLNNDR_$fast"] = dlnnidr_full
        ep_species_data["DLNTDR_$fast"] = dlntidr_full
        # Same pattern as thermal ZS (line above), but floor at 1 where fast density=0
        # to prevent ZS=0 → bb=taus*mass*(ky/0)²=Inf → singular matrix
        ep_species_data["ZS_$fast"] = max.(IMAS.avgZ(Float64(z_n), Ti_full), 1.0)

        # Extract values at selected grid points for TGLF
        Ti = Ti_full[gridpoint_cp]
        dlntidr = dlntidr_full[gridpoint_cp]
        ni = ni_full[gridpoint_cp]
        dlnnidr = dlnnidr_full[gridpoint_cp]

        println("ni for fast $fast is ", size(ni), ni)

        Zi = IMAS.avgZ(Float64(ions[iion].element[1].z_n), Ti)
        # Floor at 1 where fast density=0 → Ti=0 → avgZ=0 → ZS=0 → singular matrix
        Zi = max.(Zi, 1.0)

        setproperty!(input_tglf, Symbol("ZS_$fast"), Zi)
        setproperty!(input_tglf, Symbol("MASS_$fast"), (ions[iion].element[1].a)[1] * mp / md)
        setproperty!(input_tglf, Symbol("TAUS_$fast"), Ti ./ Te)
        setproperty!(input_tglf, Symbol("AS_$fast"), ni ./ ne)
        setproperty!(input_tglf, Symbol("VPAR_$fast"), input_tglf.VPAR_1)
        setproperty!(input_tglf, Symbol("VPAR_SHEAR_$fast"), input_tglf.VPAR_SHEAR_1)
        setproperty!(input_tglf, Symbol("RLNS_$fast"), a * dlnnidr)
        setproperty!(input_tglf, Symbol("RLTS_$fast"), a * dlntidr)

    end

    println("============================================================")

    println("Te 6 is ", ep_species_data["TEMP_6"][gridpoint_cp])
    println("dlntedr 6 is ", ep_species_data["DLNTDR_6"][gridpoint_cp])
    println("dlnnedr 6 is ", ep_species_data["DLNNDR_6"][gridpoint_cp])

    input_tglf.BETAE = 8.0 * pi .* ne .* k .* Te ./ bunit .^ 2
    loglam = 24.0 .- log.(sqrt.(ne) ./ Te)
    input_tglf.XNUE = a ./ c_s * sqrt.(ions[1].element[1].a) .* e^4 .* pi .* ne .* loglam ./ (sqrt.(me) .* (k .* Te) .^ 1.5)
    input_tglf.ZEFF = cp1d.zeff[gridpoint_cp]
    rho_s = GACODE.rho_s(cp1d, eqt)[gridpoint_cp]
    input_tglf.DEBYE = 7.43e2 * sqrt.(Te ./ ne) ./ rho_s
    input_tglf.RMIN_LOC = rmin[gridpoint_cp] ./ a
    input_tglf.RMAJ_LOC = Rmaj[gridpoint_cp] ./ a
    input_tglf.ZMAJ_LOC = 0
    input_tglf.DRMINDX_LOC = 1.0

    drmaj = IMAS.gradient(rmin, Rmaj)

    input_tglf.DRMAJDX_LOC = drmaj[gridpoint_cp]
    input_tglf.DZMAJDX_LOC = 0.0

    input_tglf.Q_LOC = abs.(q)

    #filling in EPconst


    # extraEP.E_OMEGA_TAE = sqrt.(2.0 ./ input_tglf.BETAE) ./ 2.0 ./ q ./ Rmaj
    # extraEP.E_omegaGAM = (1.0 .÷ Rmaj) .* sqrt.(1.0 .+ input_tglf.TAUS_2) .÷ (1.0 .+ 1.0 .÷ (2.0 .* q))
    # extra_EP.E_RHO_STAR = rho_s
    # extraEP.E_CS = c_s
    # extraEP.E_RMIN = rmin
    # extraEP.E_gammaE = gamma_e
    # extraEP.E_gammap = gamma_p
    # extraEP.E_NS = input_tglf.NS


    omegaGAM =  (1.0 .÷ input_tglf.RMAJ_LOC) .* sqrt.(1.0 .+ input_tglf.TAUS_2) .÷ (1.0 .+ 1.0 .÷ (2.0 .* q))
    OMEGA_TAE = sqrt.(2.0 ./ input_tglf.BETAE) ./ 2.0 ./ q ./ input_tglf.RMAJ_LOC
    RHO_STAR = rho_s
    CS = c_s
    RMIN = rmin
    NS = input_tglf.NS[1]
    

    # ===== Calculate full radial profiles for EP (before extracting selected points) =====
    c_s_full = GACODE.c_s(cp1d)
    rho_s_full = GACODE.rho_s(cp1d, eqt)
    
    # Use thermal electron profiles (not fast electrons, which overwrote ne_full/Te_full above)
    ne_thermal_full = ep_species_data["DENS_1"]
    Te_thermal_full = ep_species_data["TEMP_1"]

    # OMEGA_TAE and omegaGAM need full radial betae
    betae_full = 8.0 * pi .* ne_thermal_full .* k .* Te_thermal_full ./ (IMAS.interp1d(eq1d.rho_tor_norm, GACODE.bunit(eqt) .* T_to_Gauss).(cp1d.grid.rho_tor_norm)) .^ 2
    
    # Full radial rotation profiles for gamma_e and gamma_p
    w0_full = -cp1d.rotation_frequency_tor_sonic
    w0p_full = IMAS.gradient(rmin, w0_full)
    Rmaj_full = IMAS.interp1d(eq1d.rho_tor_norm, m_to_cm * 0.5 * (eq1d.r_outboard .+ eq1d.r_inboard)).(cp1d.grid.rho_tor_norm)
    q_full = q_profile
    gamma_e_full = -rmin ./ q_full .* w0p_full
    gamma_p_full = -Rmaj_full .* w0p_full
    
    # Calculate OMEGA_TAE and omegaGAM for full radial grid
    OMEGA_TAE_full = sqrt.(2.0 ./ betae_full) ./ (2.0 .* q_full .* Rmaj_full ./ a)
    
    # For omegaGAM, we need TAUS_2 at all radii
    if length(ions) >= 1
        Ti_ion1_full = ions[1].temperature
        omegaGAM_full = (a ./ Rmaj_full) .* sqrt.(1.0 .+ Ti_ion1_full ./ Te_thermal_full) ./ (1.0 .+ 1.0 ./ (2.0 .* q_full))
    else
        omegaGAM_full = zeros(nr_full)
    end

    input_tglf.KAPPA_LOC = kappa[gridpoint_cp]

    skappa = rmin .* IMAS.gradient(rmin, kappa) ./ kappa
    sdelta = rmin .* IMAS.gradient(rmin, delta)
    szeta = rmin .* IMAS.gradient(rmin, zeta)

    input_tglf.S_KAPPA_LOC = skappa[gridpoint_cp]
    input_tglf.DELTA_LOC = delta[gridpoint_cp]
    input_tglf.S_DELTA_LOC = sdelta[gridpoint_cp]
    input_tglf.ZETA_LOC = zeta[gridpoint_cp]
    input_tglf.S_ZETA_LOC = szeta[gridpoint_cp]

    press = cp1d.pressure_thermal
    Pa_to_dyn = 10.0

    dpdr = IMAS.gradient(rmin, press .* Pa_to_dyn)[gridpoint_cp]
    input_tglf.P_PRIME_LOC = abs.(q) ./ (rmin[gridpoint_cp] ./ a) .^ 2 .* rmin[gridpoint_cp] ./ bunit .^ 2 .* dpdr

    dqdr = IMAS.gradient(rmin, q_profile)[gridpoint_cp]
    s = rmin[gridpoint_cp] ./ q .* dqdr
    input_tglf.Q_PRIME_LOC = q .^ 2 .* a .^ 2 ./ rmin[gridpoint_cp] .^ 2 .* s

    # saturation rules
    input_tglf.ALPHA_ZF = 1.0 # 1 = default, -1 = low ky cutoff kypeak search
    input_tglf.USE_MHD_RULE = false
    input_tglf.NMODES = input_tglf.NS .+ 2 # capture main branches: ES each species + BPER + VPAR_SHEAR
    input_tglf.NKY = 12 # 12 is default, 16 for smoother spectrum
    input_tglf.ALPHA_QUENCH = 0 # 0 = spectral shift, 1 = quench
    input_tglf.SAT_RULE = parse(Int,split(string(sat),"sat")[end])
    if sat == :sat2 || sat == :sat3
        input_tglf.UNITS = "CGYRO"
        input_tglf.KYGRID_MODEL = 4
        input_tglf.NBASIS_MIN = 2
        input_tglf.NBASIS_MAX = 6
        input_tglf.USE_AVE_ION_GRID = true
        input_tglf.XNU_MODEL = 3
        input_tglf.WDIA_TRAPPED = 1.0
    else
        input_tglf.UNITS = "GYRO"
        if sat == :sat1
        elseif sat == :sat1geo
            input_tglf.UNITS = "CGYRO"
        elseif sat == :sat0quench
            input_tglf.ALPHA_QUENCH = 1
        end
        input_tglf.KYGRID_MODEL = 1
        input_tglf.NBASIS_MIN = 2
        input_tglf.NBASIS_MAX = 4
        input_tglf.USE_AVE_ION_GRID = false # default is false
        input_tglf.XNU_MODEL = 2
        input_tglf.WDIA_TRAPPED = 0.0
    end

    # electrostatic/electromagnetic
    if electromagnetic
        input_tglf.USE_BPER = true
        input_tglf.USE_BPAR = false # TGLF does not have enough moments to resolve BPAR flutter
    else
        input_tglf.USE_BPER = false
        input_tglf.USE_BPAR = false
    end

    input_tglf.ALPHA_MACH = 0.0

    # ===== Compute full-radial geometry for TJLFEP profile struct =====
    bunit_full = IMAS.interp1d(eq1d.rho_tor_norm, GACODE.bunit(eqt) .* T_to_Gauss).(cp1d.grid.rho_tor_norm)
    rmin_safe = max.(rmin, rmin[2])  # avoid division by zero at magnetic axis (rmin[1] ≈ 0)
    dqdr_full = IMAS.gradient(rmin, q_profile)
    shear_full = rmin_safe ./ q_profile .* dqdr_full
    q_prime_full = q_profile .^ 2 .* a .^ 2 ./ rmin_safe .^ 2 .* shear_full
    dpdr_full = IMAS.gradient(rmin, press .* Pa_to_dyn)
    p_prime_full = abs.(q_profile) ./ (rmin_safe ./ a) .^ 2 .* rmin_safe ./ bunit_full .^ 2 .* dpdr_full
    ir_rep = max(1, nr_full ÷ 2)

    # Build full-radial ZS matrix (nr_full × NS) from ep_species_data
    zs_full = hcat([ep_species_data["ZS_$j"] for j in 1:NS]...)
    mass_vec = vcat([getproperty(input_tglf, Symbol("MASS_$j"))[1] for j in 1:NS])


    # ===== Build extraEP dictionary with FULL radial profiles for TJLFEP =====
    extraEP = merge(ep_species_data, Dict(
        "gammaE" => gamma_e_full,
        "gammap" => gamma_p_full,
        "omegaGAM" => omegaGAM_full,
        "OMEGA_TAE" => OMEGA_TAE_full,
        "RHO_STAR" => rho_s_full,
        "CS" => c_s_full,
        "RMIN" => rmin,
        "NR" => nr_full,
        "NS" => ns_full,
        "SIGN_BT" => Float64(signb),
        "SIGN_IT" => Float64(signb * sign(q_profile[ir_rep])),
        "RMAJ" => Rmaj ./ a,
        "SHIFT" => drmaj,
        "Q" => q_profile,
        "SHEAR" => shear_full,
        "Q_PRIME" => q_prime_full,
        "P_PRIME" => p_prime_full,
        "KAPPA" => kappa,
        "S_KAPPA" => skappa,
        "DELTA" => delta,
        "S_DELTA" => sdelta,
        "ZETA" => zeta,
        "S_ZETA" => szeta,
        "BETAE" => betae_full,
        "ZEFF" => Float64.(cp1d.zeff),
        "N_ION" => ns_full - 4,
        "ZS" => zs_full,
        "MASS" => mass_vec
    ))

    return input_tglf, extraEP
end

function InputCGYRO(dd::IMAS.dd, gridpoint_cp::Integer, lump_ions::Bool)
    input_cgyro = TurbulentTransport.InputCGYRO()

    eq = dd.equilibrium
    eqt = eq.time_slice[]
    eq1d = eqt.profiles_1d
    cp1d = dd.core_profiles.profiles_1d[]

    if lump_ions
        ions = IMAS.lump_ions_as_bulk_and_impurity(cp1d)
    else
        ions = cp1d.ion
    end

    e = IMAS.cgs.e # statcoul
    k = IMAS.cgs.k # erg/eV
    mp = IMAS.cgs.mp # g
    me = IMAS.cgs.me # g
    md = 2 * mp # g
    m_to_cm = IMAS.cgs.m_to_cm
    m³_to_cm³ = IMAS.cgs.m³_to_cm³
    T_to_Gauss = IMAS.cgs.T_to_Gauss

    rmin = GACODE.r_min_core_profiles(eqt.profiles_1d, cp1d.grid.rho_tor_norm)
    a = rmin[end]

    Rmaj = IMAS.interp1d(eq1d.rho_tor_norm, m_to_cm * 0.5 * (eq1d.r_outboard .+ eq1d.r_inboard)).(cp1d.grid.rho_tor_norm)

    input_cgyro.RMIN = (rmin/a)[gridpoint_cp]
    input_cgyro.RMAJ = Rmaj[gridpoint_cp] / a

    dens_e = cp1d.electrons.density ./ m³_to_cm³
    ne = dens_e[gridpoint_cp]
    dlnnedr = -IMAS.calc_z(rmin ./ a, dens_e, :backward)
    dlnnedr = dlnnedr[gridpoint_cp]

    temp_e = cp1d.electrons.temperature
    Te = temp_e[gridpoint_cp]
    dlntedr = -IMAS.calc_z(rmin ./ a, temp_e, :backward)
    dlntedr = dlntedr[gridpoint_cp]

    n_norm = ne
    t_norm = Te
    for iion in eachindex(ions)
        species = iion

        Ti = ions[iion].temperature ./ t_norm
        dlntidr = -IMAS.calc_z(rmin ./ a, Ti, :backward)
        Ti = Ti[gridpoint_cp]
        dlntidr = dlntidr[gridpoint_cp]

        Zi = IMAS.avgZ(ions[iion].element[1].z_n, Ti*t_norm)
        setproperty!(input_cgyro, Symbol("Z_$species"), Zi)
        setproperty!(input_cgyro, Symbol("MASS_$species"), ions[iion].element[1].a .* mp / md)

        ni = ions[iion].density_thermal ./ m³_to_cm³ / n_norm
        dlnnidr = -IMAS.calc_z(rmin ./ a, ni, :backward)
        ni = ni[gridpoint_cp]
        dlnnidr = dlnnidr[gridpoint_cp]

        setproperty!(input_cgyro, Symbol("TEMP_$species"), Ti)
        setproperty!(input_cgyro, Symbol("DENS_$species"), ni)
        setproperty!(input_cgyro, Symbol("DLNNDR_$species"), dlnnidr)
        setproperty!(input_cgyro, Symbol("DLNTDR_$species"), dlntidr)
    end

    # electrons last for CGYRO
    i = length(ions) + 1
    setproperty!(input_cgyro, Symbol("DENS_$i"), ne / n_norm)
    setproperty!(input_cgyro, Symbol("TEMP_$i"), Te / t_norm)
    setproperty!(input_cgyro, Symbol("MASS_$i"), me / md)
    setproperty!(input_cgyro, Symbol("Z_$i"), -1.0)
    setproperty!(input_cgyro, Symbol("DLNNDR_$i"), dlnnedr)
    setproperty!(input_cgyro, Symbol("DLNTDR_$i"), dlntedr)

    input_cgyro.N_SPECIES = length(ions) + 1 # add 1 to include electrons

    c_s = GACODE.c_s(cp1d)[gridpoint_cp]
    loglam = 24.0 - log(sqrt(ne) / (Te))
    nu_ee = (a / c_s) * (loglam * 4 * pi * ne * e^4) / ((2 * k * Te)^(3 / 2) * me^(1 / 2))
    input_cgyro.NU_EE = nu_ee

    kappa = IMAS.interp1d(eq1d.rho_tor_norm, eq1d.elongation).(cp1d.grid.rho_tor_norm)
    input_cgyro.KAPPA = kappa[gridpoint_cp]

    skappa = rmin .* IMAS.gradient(rmin, kappa) ./ kappa
    input_cgyro.S_KAPPA = skappa[gridpoint_cp]

    drmaj = IMAS.gradient(rmin, Rmaj)
    input_cgyro.SHIFT = drmaj[gridpoint_cp]

    delta = IMAS.interp1d(eq1d.rho_tor_norm, 0.5 * (eq1d.triangularity_lower + eq1d.triangularity_upper)).(cp1d.grid.rho_tor_norm)
    input_cgyro.DELTA = delta[gridpoint_cp]
    sdelta = rmin .* IMAS.gradient(rmin, delta)
    input_cgyro.S_DELTA = sdelta[gridpoint_cp]

    zeta =
        IMAS.interp1d(
            eq1d.rho_tor_norm,
            0.25 * (eq1d.squareness_lower_inner .+ eq1d.squareness_lower_outer .+ eq1d.squareness_upper_inner .+ eq1d.squareness_upper_outer)
        ).(cp1d.grid.rho_tor_norm)
    input_cgyro.ZETA = zeta[gridpoint_cp]
    szeta = rmin .* IMAS.gradient(rmin, zeta)
    input_cgyro.S_ZETA = szeta[gridpoint_cp]

    Z0 = IMAS.interp1d(eq1d.rho_tor_norm, eq1d.geometric_axis.z * 1e2).(cp1d.grid.rho_tor_norm)
    input_cgyro.ZMAG = Z0[gridpoint_cp] / a
    sZ0 = IMAS.gradient(rmin, Z0)
    input_cgyro.DZMAG = sZ0[gridpoint_cp]

    q_profile = IMAS.interp1d(eq1d.rho_tor_norm, eq1d.q).(cp1d.grid.rho_tor_norm)
    q = q_profile[gridpoint_cp]

    input_cgyro.Q = q

    w0 = -1 * cp1d.rotation_frequency_tor_sonic
    w0p = IMAS.gradient(rmin, w0)
    gamma_p = -Rmaj[gridpoint_cp] * w0p[gridpoint_cp]
    gamma_e = -rmin[gridpoint_cp] / q * w0p[gridpoint_cp]
    mach = Rmaj[gridpoint_cp] * w0[gridpoint_cp] / c_s
    input_cgyro.GAMMA_P = (a / c_s) * gamma_p
    input_cgyro.GAMMA_E = (a / c_s) * gamma_e
    input_cgyro.MACH = mach

    bunit = IMAS.interp1d(eq1d.rho_tor_norm, GACODE.bunit(eqt) .* T_to_Gauss).(cp1d.grid.rho_tor_norm)[gridpoint_cp]

    input_cgyro.BETAE_UNIT = 8.0 * pi * ne * k * Te / bunit^2

    dqdr = IMAS.gradient(rmin, q_profile)[gridpoint_cp]
    s = rmin[gridpoint_cp] / q * dqdr
    input_cgyro.S = s

    Bt = eqt.global_quantities.vacuum_toroidal_field.b0

    input_cgyro.BTCCW = sign(Bt)
    input_cgyro.IPCCW = sign(Bt) * sign(q)

    return input_cgyro
end

"""
    load(input_tglf::InputTGLF, filename::AbstractString)

Reads filename (`input.tglf` or `input.tglf.gen` format) and populates input_tglf
"""
function load(input_tglf::InputTGLF, filename::String)
    multi_line = open(filename, "r") do file
        return read(filename, String)
    end

    filtered_lines = [line for line in split(multi_line, '\n') if occursin("=", line)]
    ip_dict = Dict(Symbol(replace(split(string(line), "=")[1], r"\s+" => "")) => replace(split(string(line), "=")[2], r"\s+" => "") for line in filtered_lines)
    field_types = fieldtypes(TurbulentTransport.InputTGLF)

    for (idx, name) in enumerate(fieldnames(typeof(input_tglf)))
        if typeof(field_types[idx]) <: Union
            type_of_item = field_types[idx].b
        else
            type_of_item = typeof(field_types[idx])
        end
        if name ∉ keys(ip_dict)
            continue
        end
        if ip_dict[name] == "T" || ip_dict[name] == ".true."
            setproperty!(input_tglf, name, true)
        elseif type_of_item == "F"
            ip_dict[name] == ".false."
            setproperty!(input_tglf, name, false)
        elseif type_of_item <: Float64
            setproperty!(input_tglf, name, parse(Float64, ip_dict[name]))
        elseif type_of_item <: Int64
            setproperty!(input_tglf, name, Int(parse(Float64, ip_dict[name])))
        elseif type_of_item <: String
            setproperty!(input_tglf, name, ip_dict[name])
        end
    end
    return input_tglf
end

"""
    run_tglf(input_tglf::InputTGLF)

Run TGLF starting from a InputTGLF.

Returns a `flux_solution` structure
"""
function run_tglf(input_tglf::InputTGLF)
    folder = mktempdir()

    save(input_tglf, joinpath(folder, "input.tglf"))

    open(joinpath(folder, "command.sh"), "w") do io
        return write(
            io,
            """
         if command -v timeout &> /dev/null; then
         	(time (timeout 120 tglf -n 1 -e .)) &> command.log
         else
         	(time (tglf -n 1 -e .)) &> command.log
         fi
         """)
    end

    fluxes = try
        run(Cmd(`bash command.sh`; dir=folder))

        tmp = open(joinpath(folder, "out.tglf.gbflux"), "r") do io
            return read(io, String)
        end

        parse_out_tglf_gbflux(tmp)

    catch e
        # show last 100 lines of  chease.output
        txt = open(joinpath(folder, "command.log"), "r") do io
            return split(read(io, String), "\n")
        end
        @error "ERROR running TGLF\n...\n" * join(txt[max(1, length(txt) - 100):end], "\n")
        rethrow(e)
    end

    sol = IMAS.flux_solution(
        fluxes["Q/Q_GB_elec"],
        fluxes["Q/Q_GB_ions"],
        fluxes["Gam/Gam_GB_elec"],
        fluxes["Gam/Gam_GB_all_ions"],
        fluxes["Pi/Pi_GB_ions"])

    rm(folder; recursive=true)

    return sol
end

"""
    run_tglf(input_tglf::InputTGLF)

Run TGLF starting from a vector of InputTGLFs.

NOTE: Each run is done asyncronously (ie. in separate parallel processes)

Returns a `flux_solution` structure
"""
function run_tglf(input_tglfs::Vector{InputTGLF})
    return collect(asyncmap(input_tglf -> TurbulentTransport.run_tglf(input_tglf), input_tglfs))
end

export run_tglf

"""
    save(input::Union{InputTGLF, InputCGYRO}, filename::AbstractString)

Write input_tglf/input_cgyro to file in input.tglf/input.cgyro/input.qlgyro format to be read by TGLF/CGYRO
"""
function save(input::Union{InputTGLF,InputCGYRO,InputQLGYRO}, filename::AbstractString)
    open(filename, "w") do io
        for key in fieldnames(typeof(input))
            if startswith(String(key), "_")
                continue
            end
            try
                value = getfield(input, key)
                if ismissing(value)
                    continue
                elseif isa(value, Int)
                    println(io, "$(key)=$(convert(Int, value))")
                elseif isa(value, String)
                    println(io, "$(key)='$value'")
                elseif isa(value, Bool)
                    println(io, "$(key)=.$value.")
                else
                    println(io, "$(key)=$(convert(Float64, value))")
                end
            catch e
                println("Error writing $key to input file")
                rethrow(e)
            end
        end
    end
end

"""
    compare_two_input_tglfs(itp_1::TurbulentTransport.InputTGLF, itp_2::TurbulentTransport.InputTGLF)

Compares two input_tglfs, prints the difference and stores the difference in a new InputTGLF
"""
function compare_two_input_tglfs(itp_1::TurbulentTransport.InputTGLF, itp_2::TurbulentTransport.InputTGLF)
    itp_diff = TurbulentTransport.InputTGLF()
    for field in fieldnames(typeof(itp_diff))
        if typeof(getproperty(itp_1, field)) <: String
            setproperty!(itp_diff, field, getproperty(itp_1, field) * "  " * getproperty(itp_2, field))
        else
            setproperty!(itp_diff, field, getproperty(itp_1, field) - getproperty(itp_2, field))
        end
    end

    for key in fieldnames(typeof(itp_diff))
        itp_1_value = getproperty(itp_1, key)
        itp_diff_value = getproperty(itp_diff, key)
        if typeof(itp_diff_value) <: String || typeof(itp_diff_value) <: Missing
            continue
        end
        println("Difference for $key = $(round(itp_diff_value/itp_1_value*100,digits=2)) % itp_2 = $(itp_diff_value+itp_1_value), itp_1 = $itp_1_value ")
    end

    return itp_diff
end

export compare_two_input_tglfs

"""
    parse_out_tglf_gbflux(lines::String; outnames=["Gam/Gam_GB", "Q/Q_GB", "Pi/Pi_GB", "S/S_GB"])

parse out.tglf.gbflux file into a dictionary with possibility of using custom names for outputs
"""
function parse_out_tglf_gbflux(lines::String; outnames::NTuple{4,String}=("Gam/Gam_GB", "Q/Q_GB", "Pi/Pi_GB", "S/S_GB"))
    vals = map(x -> parse(Float64, x), split.(lines))
    ns = Int(length(vals) / length(outnames))
    out = Dict()
    let k = 1
        for t in outnames
            out[t*"_elec"] = vals[k]
            k += 1
            out[t*"_all_ions"] = Float64[]
            for s in 1:ns-1
                out[t*"_ion$s"] = vals[k]
                push!(out[t*"_all_ions"], vals[k])
                k += 1
            end
            out[t*"_ions"] = sum(out[t*"_all_ions"])
        end
    end
    return out
end

"""
    input_qlgyro(input_qlgyro::InputQLGYRO, input_cgyro::InputCGYRO)

Run QLGYRO starting from a InputQLGYRO and InputCGYRO

Returns a `flux_solution` structure
"""
function run_qlgyro(input_qlgyro::InputQLGYRO, input_cgyro::InputCGYRO)
    folder = mktempdir()

    save(input_cgyro, joinpath(folder, "input.cgyro"))
    save(input_qlgyro, joinpath(folder, "input.qlgyro"))

    n_parallel = input_qlgyro.N_PARALLEL
    open(joinpath(folder, "command.sh"), "w") do io
        return write(
            io,
            """
         	(time (qlgyro -n $n_parallel -e .)) &> command.log
         """)
    end

    fluxes = try
        run(Cmd(`bash command.sh`; dir=folder))

        tmp = open(joinpath(folder, "out.qlgyro.gbflux"), "r") do io
            return read(io, String)
        end
        fluxes = parse_out_tglf_gbflux(tmp)

    catch e
        # show last 100 lines of  chease.output
        txt = open(joinpath(folder, "command.log"), "r") do io
            return split(read(io, String), "\n")
        end
        @error "ERROR running QLGYRO\n...\n" * join(txt[max(1, length(txt) - 100):end], "\n")
        rethrow(e)
    end

    sol = IMAS.flux_solution(
        fluxes["Q/Q_GB_elec"],
        fluxes["Q/Q_GB_ions"],
        fluxes["Gam/Gam_GB_elec"],
        fluxes["Gam/Gam_GB_all_ions"],
        fluxes["Pi/Pi_GB_ions"])

    rm(folder; recursive=true)

    return sol
end

"""
    run_qlgyro(input_qlgyros::Vector{InputQLGYRO}, input_cgyros::Vector{InputCGYRO})

Run QLGYRO starting from a vectors of InputQLGYRO and InputCGYRO

NOTE: Each run is done sequentially, one after the other

Returns a vector of `flux_solution` structures
"""
function run_qlgyro(input_qlgyros::Vector{InputQLGYRO}, input_cgyros::Vector{InputCGYRO})
    @assert length(input_qlgyros) == length(input_cgyros)
    return collect(map((input_qlgyro, input_cgyro) -> TurbulentTransport.run_qlgyro(input_qlgyro, input_cgyro), input_qlgyros, input_cgyros))
end

export run_qlgyro
