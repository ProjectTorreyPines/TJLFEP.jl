# Single-radius, file-based run (N_BASIS=16, SCAN_N=1, ir=2).

get(ENV, "TJLFEP_FILE_ONLY", "0") == "1" || (ENV["TJLFEP_FILE_ONLY"] = "1")

using Pkg
Pkg.activate("..")

function logmsg(args...)
    println(args...)
    flush(stdout)
    flush(stderr)
end

using TJLFEP
using TJLF
using LinearAlgebra

BLAS.set_num_threads(1)

const TJLFEP_ROOT = normpath(@__DIR__, "..")
const DEBUG_DIR = normpath(@__DIR__, "debug_nb16")
const CASE_DIR = get(ENV, "CASE_DIR",
    joinpath(TJLFEP_ROOT, "src", "DIIIDfiles", "202017C42_500ms_v3.1"))
const TGLFEP_FILE = get(ENV, "TGLFEP_FILE", joinpath(DEBUG_DIR, "input.TGLFEP"))
const FILE_DIR = get(ENV, "FILE_DIR",
    joinpath(DEBUG_DIR, "fileInput_$(get(ENV, "SLURM_JOB_ID", "local"))"))

@assert isfile(joinpath(CASE_DIR, "dump.profile"))
@assert isfile(TGLFEP_FILE)

use_gpu = TJLF.pick_device(:auto) === :gpu
logmsg("TJLFEP_DEBUG=", get(ENV, "TJLFEP_DEBUG", "0"))
logmsg("N_BASIS=16 SCAN_N=1 (single radius)")
logmsg("device: ", use_gpu ? "GPU" : "CPU")
logmsg("CASE_DIR: ", CASE_DIR)
logmsg("TGLFEP_FILE: ", TGLFEP_FILE)
logmsg("FILE_DIR: ", FILE_DIR)

prof, ir_exp = setup_fortran_file_inputs(CASE_DIR, FILE_DIR; tglfep_file=TGLFEP_FILE)
logmsg("IR_EXP: ", ir_exp)

tglfep = abspath(joinpath(FILE_DIR, "input.TGLFEP"))
mtglf = abspath(joinpath(FILE_DIR, "input.MTGLF"))
expro = abspath(joinpath(FILE_DIR, "input.EXPRO"))
outdir = abspath(joinpath(@__DIR__, "debug_out_nb16_$(get(ENV, "SLURM_JOB_ID", "local"))"))
mkpath(outdir)

t0 = time()
cd(outdir) do
    logmsg("--- runTHD start (threads, 1 radius) ---")
    @time runTHD(tglfep, mtglf, expro; printout=true, use_gpu=use_gpu, parallel=:threads)
    make_crit_grad_plots("neither"; dir=".", code="julia")
end
logmsg("debug_compare OK in $(round(time() - t0; digits=1)) s -> $outdir")
