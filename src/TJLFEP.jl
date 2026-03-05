module TJLFEP
using Base.Threads
using LinearAlgebra # Don't really need this for tjlf-ep
using SparseArrays
using StaticArrays
using TJLF
using IMAS
using TurbulentTransport
using FUSE

include("tjlfep_modules.jl")
include("tjlfep_read_inputs.jl")
include("EXPROconst.jl")
include("tjlfep_ky.jl")
include("tjlfep_kwscale_scan.jl")
include("mainsub.jl")
include("conv_input.jl")
include("tjlfep_complete_output.jl")
include("run_tjlfep.jl")

include("tjlfep_generate_input.jl")

include("tglf_(context).jl")
include("tglf_actor_(context).jl")

#sgould this line just be the big struct?
export profile, Options, InputTJLF, makeStructs
export readMTGLF, readTGLFEP, TJLF_map, readEXPRO
export convert_input, revert_input
export tjlfep_complete_output
export runTHD, runTHDs

export TJLFEP_generate_input, readline_values

export InputTGLFEP
export populate_tjlfep_profile!

end 