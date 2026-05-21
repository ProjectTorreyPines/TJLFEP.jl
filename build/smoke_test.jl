# Single-process DIII-D smoke test (no SlurmClusterManager).
# Run from TJLFEP/build:
#   julia --project=.. --sysimage=noTJLF_TJLFEP_sysimage.so smoke_test.jl

using Pkg
Pkg.activate("..")

const TJLFEP_ROOT = normpath(@__DIR__, "..")
const SYSIMAGE = joinpath(@__DIR__, "noTJLF_TJLFEP_sysimage.so")
if isfile(SYSIMAGE)
    println("Using sysimage: ", SYSIMAGE)
else
    @warn "Sysimage missing; running with standard precompile" SYSIMAGE
end

using TJLFEP
using TJLF
using GACODE
using IMAS
using LinearAlgebra

BLAS.set_num_threads(1)
TJLF.pick_device(:auto)

const INPUT = joinpath(TJLFEP_ROOT, "src", "DIIIDfiles", "202017C42_500ms_v3.1", "input.gacode")
@assert isfile(INPUT) "Missing test input: $INPUT"

SCAN_N = 1
rho = [0.5]
N_BASIS = 8  # smaller basis for faster smoke test
use_gpu = (TJLF.pick_device(:auto) === :gpu)

OptionsDict = Dict{String, Any}(
    "nn" => 5, "nr" => 101, "jtscale_max" => 1, "nmodes" => 4,
    "PROCESS_IN" => 5, "THRESHOLD_FLAG" => 0, "N_BASIS" => N_BASIS,
    "SCAN_METHOD" => 2, "REJECT_I_PINCH_FLAG" => 0, "REJECT_E_PINCH_FLAG" => 0,
    "REJECT_TH_PINCH_FLAG" => 0, "REJECT_EP_PINCH_FLAG" => 0,
    "REJECT_TEARING_FLAG" => 1, "ROTATIONAL_SUPPRESSION_FLAG" => 0,
    "PPRIME_METHOD" => 3, "QL_RATIO_THRESH" => 10.0, "THETA_SQ_THRESH" => 100.0,
    "Q_SCALE" => 1.0, "WRITE_WAVEFUNCTION" => 0, "KY_MODEL" => 2,
    "SCAN_N" => SCAN_N, "IRS" => 2, "FACTOR_IN_PROFILE" => false, "FACTOR_IN" => 10.0,
    "WIDTH_IN_FLAG" => false, "WIDTH_MIN" => 1.0, "WIDTH_MAX" => 2.0,
    "INPUT_PROFILE_METHOD" => 2, "N_ION" => 2, "IS_EP" => 1, "REAL_FREQ" => 1,
)

println("device: ", use_gpu ? "GPU" : "CPU")
println("input: ", INPUT)
println("rho: ", rho)

inputGACODE = GACODE.load(INPUT)
dd = IMAS.dd(inputGACODE)

outdir = joinpath(@__DIR__, "smoke_out_$(get(ENV, "SLURM_JOB_ID", "local"))")
mkpath(outdir)

t0 = time()
cd(outdir) do
    width, kymark, SFmin, dpdr_crit, dndr_crit = runTHD(
        dd, rho, OptionsDict;
        printout = true,
        saveFiles = false,
        dir = joinpath(TJLFEP_ROOT, "src", "DIIIDfiles", "fileInput"),
        use_gpu = use_gpu,
    )
    println("SFmin = ", SFmin)
    println("width = ", width)
    println("kymark = ", kymark)
    @assert all(isfinite, SFmin) "SFmin not finite: $SFmin"
    @assert all(SFmin .< 9000) "No stable threshold found (sentinel 10000)"
end
dt = time() - t0
println("smoke_test OK in $(round(dt; digits=1)) s; outputs in $outdir")
