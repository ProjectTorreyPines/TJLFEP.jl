using Pkg
Pkg.activate("..")
using TJLFEP, TJLF, LinearAlgebra
BLAS.set_num_threads(1)

const CASE_DIR = ENV["CASE_DIR"]
const TGLFEP_FILE = ENV["TGLFEP_FILE"]
const FILE_DIR = ENV["FILE_DIR"]

setup_fortran_file_inputs(CASE_DIR, FILE_DIR; tglfep_file=TGLFEP_FILE)
outdir = joinpath(@__DIR__, "debug_out_nb6_scan20_$(get(ENV, "SLURM_JOB_ID", "local"))")
mkpath(outdir)

cd(outdir) do
    @time runTHD(
        abspath(joinpath(FILE_DIR, "input.TGLFEP")),
        abspath(joinpath(FILE_DIR, "input.MTGLF")),
        abspath(joinpath(FILE_DIR, "input.EXPRO"));
        printout=true, use_gpu=false, parallel=:threads)
    make_crit_grad_plots("FILE"; dir=outdir)
end
println("OK -> ", outdir)
