# Full module load (IMAS + file paths). TJLFEP.jl uses run_tjlfep_file.jl alone when TJLFEP_FILE_ONLY=1.
include("run_tjlfep_file.jl")
include("run_tjlfep_imas.jl")
