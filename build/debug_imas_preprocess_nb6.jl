# Preprocessing-only DIII-D check: IMAS input construction vs Fortran dump.profile path.
# No TJLF eigensolver. Writes ref + test input.MTGLF / input.EXPRO / input.TGLFEP and runs compare.
#
# Run:
#   TJLFEP_FILE_ONLY=0 julia --startup-file=no --project=.. build/debug_imas_preprocess_nb6.jl

# Must be set before TJLFEP loads (otherwise run_tjlfep_imas.jl is skipped).
ENV["TJLFEP_FILE_ONLY"] = "0"

using Pkg
Pkg.activate("..")

using TJLFEP
using TJLFEP: preprocess_imas_inputs, save_imas_preprocessed_inputs, setup_fortran_file_inputs
using IMAS
using GACODE
using Printf

const TJLFEP_ROOT = normpath(@__DIR__, "..")
const CASE_DIR = get(ENV, "CASE_DIR",
    joinpath(TJLFEP_ROOT, "src", "DIIIDfiles", "202017C42_500ms_v3.1"))
const TGLFEP_FILE = get(ENV, "TGLFEP_FILE",
    joinpath(TJLFEP_ROOT, "build", "debug_nb6", "input.TGLFEP"))
const REF_DIR = get(ENV, "REF_DIR",
    joinpath(TJLFEP_ROOT, "build", "imas_ref_local"))
const TEST_DIR = get(ENV, "TEST_DIR",
    joinpath(TJLFEP_ROOT, "build", "imas_test_local"))

const rho = [0.01]
const IS_EP = 1

function logmsg(args...)
    println(args...)
    flush(stdout)
    flush(stderr)
end

@assert isfile(joinpath(CASE_DIR, "dump.profile"))
@assert isfile(joinpath(CASE_DIR, "input.gacode"))
@assert isfile(TGLFEP_FILE)

OptionsDict = Dict{String, Any}(
    "nn" => 5,
    "nr" => 101,
    "jtscale_max" => 1,
    "nmodes" => 4,
    "PROCESS_IN" => 5,
    "THRESHOLD_FLAG" => 0,
    "N_BASIS" => 6,
    "SCAN_METHOD" => 2,
    "REJECT_I_PINCH_FLAG" => 0,
    "REJECT_E_PINCH_FLAG" => 0,
    "REJECT_TH_PINCH_FLAG" => 0,
    "REJECT_EP_PINCH_FLAG" => 0,
    "REJECT_TEARING_FLAG" => 1,
    "ROTATIONAL_SUPPRESSION_FLAG" => 0,
    "PPRIME_METHOD" => 3,
    "QL_RATIO_THRESH" => 10.0,
    "THETA_SQ_THRESH" => 100.0,
    "Q_SCALE" => 1.0,
    "WRITE_WAVEFUNCTION" => 1,
    "KY_MODEL" => 2,
    "SCAN_N" => 1,
    "IRS" => 2,
    "FACTOR_IN_PROFILE" => false,
    "FACTOR_IN" => 10.0,
    "WIDTH_IN_FLAG" => false,
    "WIDTH_MIN" => 1.0,
    "WIDTH_MAX" => 2.0,
    "INPUT_PROFILE_METHOD" => 2,
    "N_ION" => 2,
    "IS_EP" => IS_EP,
    "REAL_FREQ" => 1,
)

logmsg("=== IMAS preprocessing-only (DIII-D nb6, SCAN_N=1) ===")
logmsg("CASE_DIR: ", CASE_DIR)
logmsg("TGLFEP_FILE: ", TGLFEP_FILE)
logmsg("REF_DIR: ", REF_DIR)
logmsg("TEST_DIR: ", TEST_DIR)
logmsg("rho: ", rho)

logmsg("--- reference: setup_fortran_file_inputs (dump.profile) ---")
setup_fortran_file_inputs(CASE_DIR, REF_DIR; tglfep_file=TGLFEP_FILE)

logmsg("--- test: IMAS preprocess_imas_inputs + save ---")
# Same IMAS entry as DIIID_juliaValidation.jl (colleague workflow):
inputFile = joinpath(CASE_DIR, "input.gacode")
@time inputGACODE = GACODE.load(inputFile)
@time dd = IMAS.dd(inputGACODE)
Options, profile, extraEP, _ = preprocess_imas_inputs(dd, rho, OptionsDict; verbose=true)
save_imas_preprocessed_inputs(Options, profile, extraEP, TEST_DIR)

logmsg("IR_EXP: ", Options.IR_EXP)
logmsg("profile NR/NS: ", profile.NR, " / ", profile.NS)
logmsg("ref NS (from dump.profile): 3 ; IMAS NS: ", extraEP["NS"])

compare_script = joinpath(@__DIR__, "compare_preprocess_inputs.jl")
if isfile(compare_script)
    logmsg("--- compare ---")
    ENV["REF_DIR"] = REF_DIR
    ENV["TEST_DIR"] = TEST_DIR
    include(compare_script)
else
    logmsg("compare script missing: ", compare_script)
end

logmsg("=== done ===")
