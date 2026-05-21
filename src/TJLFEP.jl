module TJLFEP
# File-based validation (dump.profile + runTHD string): set TJLFEP_FILE_ONLY=1 to skip IMAS/FUSE/TurbulentTransport.
const _FILE_ONLY = get(ENV, "TJLFEP_FILE_ONLY", "0") == "1"

using Base.Threads
using LinearAlgebra # Don't really need this for tjlf-ep
# using CUDA
using SparseArrays
using Printf
using StaticArrays
using TJLF
if !_FILE_ONLY
    using IMAS
    using TurbulentTransport
    using FUSE
end
using Plots

include("tjlfep_modules.jl")
include("tjlfep_ad_extensions.jl")
include("tjlfep_read_inputs.jl")
include("tjlfep_debug.jl")
include("EXPROconst.jl")
include("tjlfep_ky.jl")
include("tjlfep_kwscale_scan.jl")
include("mainsub.jl")
include("conv_input.jl")
include("tjlfep_complete_output.jl")
include("run_tjlfep_file.jl")
if !_FILE_ONLY
    include("run_tjlfep_imas.jl")
end

include("tjlfep_generate_input.jl")

if !_FILE_ONLY
    include("context.jl")
    include("actor_context.jl")
end

include("plotCritGrads.jl")

#sgould this line just be the big struct?
export profile, Options, InputTJLF
export readMTGLF, readTGLFEP, TJLF_map, readEXPRO, save_TGLFEP, save_MTGLF, save_EXPRO, save_all
export read_input_profile, readprofile, expro_dict_from_profile, ir_exp_from_scan, setup_fortran_file_inputs, setup_gacode_file_inputs, profile_from_gacode
export expro_vectors_from_gacode, preprocess_gacode_inputs
export expro_species_for_gacode_is_ep, read_gacode_scalar_field, read_gacode_ion_field
export expro_bound_deriv, expro_log_gradients, read_expro_for_alpha, compute_alpha_crit_profiles
export convert_input, revert_input
export tjlfep_complete_output
export runTHD, runTHD_from_gacode, runTHDs
export run_gacode_scan_task, finalize_gacode_scan, slurm_array_task_id

export TJLFEP_generate_input, readline_values

if !_FILE_ONLY
    export InputTGLFEP
    export populate_tjlfep_profile!, update_input_tjlfep!
    export preprocess_imas_inputs, save_imas_preprocessed_inputs, remap_extraEP_for_fortran_save!
end

export make_crit_grad_plots
export TJLF

end 