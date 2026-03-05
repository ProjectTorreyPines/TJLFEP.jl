using Revise
using Pkg
using Plots
Pkg.activate("..")
# println("Post activate, pre instantiate")
# include("../../TJLF/src/TJLF.jl")
# println("Post instantiate, pre include")
include("TJLFEP.jl")
Pkg.instantiate()
println("Post include")
using .TJLFEP
using .TJLFEP: convert_input
using .TJLFEP: revert_input
using TJLF
using Base.Threads
using LinearAlgebra
using Dates
BLAS.set_num_threads(1)
begin
    homedirectory = pwd()
    mtglffilepath = homedirectory*"/../tests/tglfep_tests/input.MTGLF"
    exprofilepath = homedirectory*"/../tests/tglfep_tests/input.EXPRO"

    OptionsDict = Dict{String, Any}("nn" => 5, "nr" => 201, "jtscale_max" => 1, "nmodes" => 4,
    "PROCESS_IN" => 5, "THRESHOLD_FLAG" => 0, "N_BASIS" => 2,
    "SCAN_METHOD" => 1, "REJECT_I_PINCH_FLAG" => 0, "REJECT_E_PINCH_FLAG" => 0, "REJECT_TH_PINCH_FLAG" => 1, "REJECT_EP_PINCH_FLAG" => 0,
    "REJECT_TEARING_FLAG" => 1, "ROTATIONAL_SUPPRESSION_FLAG" => 1, "QL_RATIO_THRESH" => 0.001, "THETA_2_THRESH" => 100.0, "Q_SCALE" => 1.0,
    "WRITE_WAVEFUNCTION" => 1, "KY_MODEL" => 2, "SCAN_N" => 6, "IRS" => 2, "FACTOR_IN_PROFILE" => false, "FACTOR_IN" => 1.0,
    "WIDTH_IN_FLAG" => false, "WIDTH_MIN" => 1.0, "WIDTH_MAX" => 2.0, "INPUT_PROFILE_METHOD" => 2, "N_ION" => 3, "IS_EP" => 3, "REAL_FREQ" => 1)
    
    profile, inputTJLFEP, params, ir_exp = makeStructs(mtglffilepath, OptionsDict, exprofilepath)
    runTHD(profile, inputTJLFEP, params, ir_exp; printout = true)
end
# I will now run runTHDs on two examples:
#=
tglfepfilepath1 = homedirectory*"/outputs/tjlfeptests/isEP3v6/input.TGLFEP"
tglfepfilepath2 = homedirectory*"/outputs/tjlfeptests/isEP3v7/input.TGLFEP"

mtglffilepath1 = homedirectory*"/outputs/tjlfeptests/isEP3v6/input.MTGLF"
mtglffilepath2 = homedirectory*"/outputs/tjlfeptests/isEP3v7/input.MTGLF"

exprofilepath1 = homedirectory*"/outputs/tjlfeptests/isEP3v6/input.EXPRO"
exprofilepath2 = homedirectory*"/outputs/tjlfeptests/isEP3v7/input.EXPRO"

tglfeps = [tglfepfilepath1, tglfepfilepath2]
mtglfs = [mtglffilepath1, mtglffilepath2]
expros = [exprofilepath1, exprofilepath2]
=#
#widths, kymark_outs, SFmins, dpdr_crit_outs, dndr_crit_outs = runTHDs(tglfeps, mtglfs, expros)