# N_BASIS=32, SCAN_N=20, multi-node (SlurmClusterManager + pmap).
# Run: sbatch batch_debug_nb32_julia_scan20_10n.sh

get(ENV, "TJLFEP_FILE_ONLY", "0") == "1" || (ENV["TJLFEP_FILE_ONLY"] = "1")

using Pkg
Pkg.activate("..")

function logmsg(args...)
    println(args...)
    flush(stdout)
    flush(stderr)
    return nothing
end

const TJLFEP_ROOT = normpath(@__DIR__, "..")
const DEBUG_DIR = normpath(@__DIR__, "debug_nb32")
const CASE_DIR = get(ENV, "CASE_DIR",
    joinpath(TJLFEP_ROOT, "src", "DIIIDfiles", "202017C42_500ms_v3.1"))
const FILE_DIR = get(ENV, "FILE_DIR",
    joinpath(DEBUG_DIR, "fileInput_scan20_10n_$(get(ENV, "SLURM_JOB_ID", "local"))"))
const SCAN_N = parse(Int, get(ENV, "SCAN_N", "20"))
const THREADS_PER_WORKER = parse(Int, get(ENV, "JULIA_WORKER_THREADS",
    get(ENV, "SLURM_CPUS_PER_TASK", "64")))
const TGLFEP_FILE = get(ENV, "TGLFEP_FILE", joinpath(DEBUG_DIR, "input_scan20.TGLFEP"))
const GACODE_DUMP = get(ENV, "GACODE_DUMP", joinpath(CASE_DIR, "input.gacode"))
if !haskey(ENV, "GACODE_DUMP")
    ENV["GACODE_DUMP"] = GACODE_DUMP
end

@assert isfile(joinpath(CASE_DIR, "dump.profile"))
@assert isfile(TGLFEP_FILE)
@assert isfile(GACODE_DUMP)

logmsg("=== nb32 scan20 distributed N_BASIS=32 SCAN_N=$SCAN_N ===")
Pkg.instantiate()
t0 = time()
Pkg.precompile()
logmsg("Pkg.precompile done in ", round(time() - t0; digits=1), " s")

using Distributed
using SlurmClusterManager
using TJLFEP
using TJLF
using LinearAlgebra

BLAS.set_num_threads(1)

project_path = TJLFEP_ROOT
_sysimage = get(ENV, "TJLFEP_SYSIMAGE", joinpath(TJLFEP_ROOT, "build", "noTJLF_TJLFEP_sysimage.so"))
exeflags = if isfile(_sysimage)
    `--project=$(project_path) --sysimage=$(_sysimage) -t $(THREADS_PER_WORKER)`
else
    `--project=$(project_path) -t $(THREADS_PER_WORKER)`
end
logmsg("worker exeflags sysimage=", isfile(_sysimage) ? _sysimage : "none")

worker_env = Dict{String,String}()
haskey(ENV, "JULIA_DEPOT_PATH") && (worker_env["JULIA_DEPOT_PATH"] = ENV["JULIA_DEPOT_PATH"])
worker_env["JULIA_PKG_PRECOMPILE_AUTO"] = "0"
worker_env["TJLFEP_FILE_ONLY"] = "1"
worker_env["GACODE_DUMP"] = GACODE_DUMP

if haskey(ENV, "SLURM_JOB_ID") || haskey(ENV, "SLURM_JOBID")
    ntasks = parse(Int, get(ENV, "SLURM_NTASKS", string(SCAN_N)))
    @assert ntasks == SCAN_N
    logmsg("SlurmClusterManager: ntasks=$ntasks threads=$THREADS_PER_WORKER")
    addprocs(SlurmManager(); exeflags=exeflags, env=worker_env)
else
    addprocs(SCAN_N; exeflags=exeflags, env=worker_env)
end

@everywhere begin
    get(ENV, "TJLFEP_FILE_ONLY", "0") == "1" || (ENV["TJLFEP_FILE_ONLY"] = "1")
    using TJLFEP
    using TJLF
    using LinearAlgebra
    BLAS.set_num_threads(1)
end

if !isfile(joinpath(FILE_DIR, "input.MTGLF"))
    setup_fortran_file_inputs(CASE_DIR, FILE_DIR; tglfep_file=TGLFEP_FILE)
end

use_gpu = false
tglfep = abspath(joinpath(FILE_DIR, "input.TGLFEP"))
mtglf = abspath(joinpath(FILE_DIR, "input.MTGLF"))
expro = abspath(joinpath(FILE_DIR, "input.EXPRO"))
outdir = abspath(joinpath(@__DIR__, "debug_out_nb32_scan20_$(get(ENV, "SLURM_JOB_ID", "local"))_dist"))
mkpath(outdir)

logmsg("workers=", nworkers(), " FILE_DIR=", FILE_DIR, " output -> ", outdir)

t0 = time()
cd(outdir) do
    @time runTHD(tglfep, mtglf, expro; printout=false, use_gpu=use_gpu, parallel=:distributed)
    make_crit_grad_plots("neither"; dir=".", code="julia")
end
logmsg("OK in ", round(time() - t0; digits=1), " s")
