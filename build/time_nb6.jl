#!/usr/bin/env julia
# Wall time for N_BASIS=6, SCAN_N=1 via runTHD_from_gacode (matches Fortran gacode driver path).
# Prints a single parseable line: TIMING_RESULT backend=... device=... seconds=...

ENV["TJLFEP_FILE_ONLY"] = "1"
ENV["TJLFEP_DEBUG"] = get(ENV, "TJLFEP_DEBUG", "0")

using Pkg
Pkg.activate(normpath(@__DIR__, ".."))

use_gpu = get(ENV, "USE_GPU", "") == "1"
if use_gpu
    using CUDA
end

using Printf
using TJLFEP
using TJLF

const ROOT = normpath(@__DIR__, "..")
const GACODE = joinpath(ROOT, "src", "DIIIDfiles", "202017C42_500ms_v3.1", "input.gacode")
const TGLFEP = joinpath(ROOT, "build", "debug_nb6", "input.TGLFEP")

@assert isfile(GACODE)
@assert isfile(TGLFEP)

opts, _, _ = preprocess_gacode_inputs(GACODE, TGLFEP)
@assert opts.N_BASIS == 6
@assert opts.SCAN_N == 1

device = use_gpu ? "gpu" : "cpu"
println("TIMING_START backend=julia device=$device SCAN_N=$(opts.SCAN_N) N_BASIS=$(opts.N_BASIS) threads=$(Threads.nthreads())")
flush(stdout)

wall_s = @elapsed runTHD_from_gacode(
    GACODE, TGLFEP; printout=false, use_gpu=use_gpu, parallel=:threads)

println(@sprintf("TIMING_RESULT backend=julia device=%s seconds=%.3f SCAN_N=%d N_BASIS=%d threads=%d",
    device, wall_s, opts.SCAN_N, opts.N_BASIS, Threads.nthreads()))
flush(stdout)
