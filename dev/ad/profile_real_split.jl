# Sample-profile ONE real gamma_dgamma_dfactor call at nb=32 to split the
# per-solve cost between the Dual-typed matrix assembly (get_matrix) and the
# eigensolve rule (_herm_eigen / geev). This decides what to optimize.
#
#   module load julia/1.11.7 && export JULIA_DEPOT_PATH=$PSCRATCH/.julia
#   julia --project=. dev/ad/profile_real_split.jl

using TJLFEP
using Profile
using Printf

const CASE = joinpath(@__DIR__, "..", "..", "examples", "DIIID_202017C42_500ms_v3.1")

function main()
    opts, prof, _ = preprocess_gacode_inputs(joinpath(CASE, "input.gacode"),
                                             joinpath(CASE, "input.TGLFEP"))
    opts.IR = 38; opts.N_BASIS = 32
    opts.KYHAT_IN = 0.25; opts.WIDTH_IN = 1.0; opts.FACTOR_IN = 2.5

    println("warm up...")
    gamma_dgamma_dfactor(opts, prof)

    Profile.clear()
    Profile.init(; n = 10_000_000, delay = 0.005)
    t = @elapsed (@profile gamma_dgamma_dfactor(opts, prof))
    @printf("one Dual solve = %.1f s\n", t)

    open(joinpath(@__DIR__, "prof_flat.txt"), "w") do io
        Profile.print(io; format = :flat, sortedby = :count, mincount = 20)
    end
    println("wrote prof_flat.txt")
end

main()
