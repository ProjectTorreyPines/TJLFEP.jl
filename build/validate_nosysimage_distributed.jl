# File-based DIII-D validation with multi-node parallelism (SlurmClusterManager + pmap).
# One Julia worker per radius (SLURM_NTASKS=SCAN_N). Threads are only used inside each
# worker (kwscale_scan); runTHD(..., parallel=:distributed) does not nest @threads over radii.
#
# Run inside: sbatch batch_prod_nb32_julia.sh

get(ENV, "TJLFEP_FILE_ONLY", "0") == "1" || (ENV["TJLFEP_FILE_ONLY"] = "1")

using Pkg
Pkg.activate("..")

"""Print and flush (Slurm .out is not a TTY; otherwise logs appear only at job end)."""
function logmsg(args...)
    println(args...)
    flush(stdout)
    flush(stderr)
    return nothing
end

const TJLFEP_ROOT = normpath(@__DIR__, "..")
const CASE_DIR = get(ENV, "CASE_DIR",
    joinpath(TJLFEP_ROOT, "src", "DIIIDfiles", "202017C42_500ms_v3.1"))
const FILE_DIR = get(ENV, "FILE_DIR",
    joinpath(@__DIR__, "fileInput_$(get(ENV, "SLURM_JOB_ID", "local"))"))
const SCAN_N = parse(Int, get(ENV, "SCAN_N", "20"))
const THREADS_PER_WORKER = parse(Int, get(ENV, "JULIA_WORKER_THREADS",
    get(ENV, "SLURM_CPUS_PER_TASK", "64")))

const TGLFEP_FILE = get(ENV, "TGLFEP_FILE",
    joinpath(TJLFEP_ROOT, "build", "debug_prod", "input.TGLFEP"))
const GACODE_DUMP = get(ENV, "GACODE_DUMP", joinpath(CASE_DIR, "input.gacode"))
if !haskey(ENV, "GACODE_DUMP")
    ENV["GACODE_DUMP"] = GACODE_DUMP
end

@assert isfile(joinpath(CASE_DIR, "dump.profile"))
@assert isfile(TGLFEP_FILE)
@assert isfile(GACODE_DUMP)

logmsg("=== precompile on manager (depot=", get(ENV, "JULIA_DEPOT_PATH", "<default>"), ") ===")
logmsg("Pkg.instantiate...")
Pkg.instantiate()
logmsg("Pkg.precompile...")
t0 = time()
Pkg.precompile()
logmsg("Pkg.precompile done in ", round(time() - t0; digits=1), " s")

using Distributed
using SlurmClusterManager
using TJLFEP
using TJLF
using LinearAlgebra

BLAS.set_num_threads(1)
@assert occursin("/dev/TJLF", pathof(TJLF)) "expected dev TJLF, got $(pathof(TJLF))"
logmsg("TJLF: ", pathof(TJLF))
logmsg("TJLFEP: ", pathof(TJLFEP))
logmsg("=== manager packages loaded ===")

project_path = TJLFEP_ROOT
exeflags = `--project=$(project_path) -t $(THREADS_PER_WORKER)`

worker_env = Dict{String,String}()
if haskey(ENV, "JULIA_DEPOT_PATH")
    worker_env["JULIA_DEPOT_PATH"] = ENV["JULIA_DEPOT_PATH"]
end
worker_env["JULIA_PKG_PRECOMPILE_AUTO"] = "0"
worker_env["TJLFEP_FILE_ONLY"] = get(ENV, "TJLFEP_FILE_ONLY", "1")
worker_env["GACODE_DUMP"] = GACODE_DUMP

if haskey(ENV, "SLURM_JOB_ID") || haskey(ENV, "SLURM_JOBID")
    ntasks = parse(Int, get(ENV, "SLURM_NTASKS", string(SCAN_N)))
    @assert ntasks == SCAN_N "SLURM_NTASKS ($ntasks) should equal SCAN_N ($SCAN_N) for one radius per worker"
    logmsg("Adding SlurmClusterManager workers (SLURM_NTASKS=$ntasks, threads=$THREADS_PER_WORKER)...")
    addprocs(SlurmManager(); exeflags=exeflags, env=worker_env)
else
    logmsg("No SLURM allocation — adding local workers for testing")
    addprocs(SCAN_N; exeflags=exeflags, env=worker_env)
end

logmsg("Loading packages on $(nworkers()) workers...")
t_load = time()
@everywhere begin
    get(ENV, "TJLFEP_FILE_ONLY", "0") == "1" || (ENV["TJLFEP_FILE_ONLY"] = "1")
    using TJLFEP
    using TJLF
    using LinearAlgebra
    BLAS.set_num_threads(1)
    function logmsg(args...)
        println(args...)
        flush(stdout)
        flush(stderr)
        return nothing
    end
end
logmsg("Worker package load finished in ", round(time() - t_load; digits=1), " s")

@everywhere logmsg("worker $(myid()) on $(gethostname()) — TJLF ", pathof(TJLF))
logmsg("Launched $(nworkers()) workers")

if !isfile(joinpath(FILE_DIR, "input.MTGLF")) || !isfile(joinpath(FILE_DIR, "input.EXPRO"))
    logmsg("--- setup_fortran_file_inputs (main) ---")
    setup_fortran_file_inputs(CASE_DIR, FILE_DIR; tglfep_file=TGLFEP_FILE)
    logmsg("--- setup_fortran_file_inputs done ---")
end

use_gpu = TJLF.pick_device(:auto) === :gpu
logmsg("device: ", use_gpu ? "GPU" : "CPU")
logmsg("workers: ", nworkers(), "  SCAN_N: ", SCAN_N)
logmsg("CASE_DIR: ", CASE_DIR)
logmsg("FILE_DIR: ", FILE_DIR)

tglfep = abspath(joinpath(FILE_DIR, "input.TGLFEP"))
mtglf = abspath(joinpath(FILE_DIR, "input.MTGLF"))
expro = abspath(joinpath(FILE_DIR, "input.EXPRO"))
outdir = abspath(joinpath(@__DIR__, "validate_out_$(get(ENV, "SLURM_JOB_ID", "local"))_files_dist"))
mkpath(outdir)

logmsg("output dir: ", outdir)
t0 = time()
cd(outdir) do
    logmsg("--- runTHD (distributed pmap) start ---")
    @time runTHD(tglfep, mtglf, expro; printout=true, use_gpu=use_gpu, parallel=:distributed)
    logmsg("--- runTHD done; writing plots ---")
    make_crit_grad_plots("neither"; dir=".", code="julia")
end
logmsg("validate (distributed) OK in $(round(time() - t0; digits=1)) s -> $outdir")
