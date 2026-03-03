using Revise
using Pkg
using Plots
Pkg.activate("..")
# println("Post activate, pre instantiate")
# include("../../TJLF/src/TJLF.jl")
# println("Post instantiate, pre include")
include("TJLFEP.jl")
# Pkg.instantiate()
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

    
    tglfepfilepath = homedirectory*"/../tests/tglfep_tests/input.TGLFEP"
    mtglffilepath = homedirectory*"/../tests/tglfep_tests/input.MTGLF"
    exprofilepath = homedirectory*"/../tests/tglfep_tests/input.EXPRO"

    # tjlf_ep_input = TJLF_EP_Input(dd)    

    runTHD(tglfepfilepath, mtglffilepath, exprofilepath, printout = true)
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
