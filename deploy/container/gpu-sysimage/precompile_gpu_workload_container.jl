# Precompile workload for the container GPU sysimage (CUDA + TJLF + TJLFEP
# standalone — the lean image has no IMAS/GACODE/TurbulentTransport, so the
# TJLFEPIMASExt extension stays dormant and is not baked).
#
# Adapted from build/sysimage/precompile_gpu_workload_fileonly.jl; the only
# delta is that example paths come from pkgdir(TJLFEP) (the package lives in
# the image depot, not a repo checkout). One nb6 radius per solver on the GPU
# bakes the tjlf_LS eigenvalue + eigenvector CUSOLVER paths, the batched
# shift-invert driver, and the AD/robust/truth critical-factor solvers.
using CUDA
using TJLF
using TJLFEP
using LinearAlgebra
BLAS.set_num_threads(1)

if !CUDA.functional()
    error("precompile_gpu_workload_container.jl must run on a GPU node (CUDA.functional() == false)")
end
@info "container precompile workload GPU" name = CUDA.name(first(CUDA.devices()))

# NB: keep these as LOCALS inside a `let` (not top-level `const`) — a baked
# `const GACODE` in Main would collide with the MPS task script's own
# `const GACODE = ENV[...]` at runtime.
let
    EX     = joinpath(pkgdir(TJLFEP), "examples", "DIIID_202017C42_500ms_v3.1")
    GACODE = joinpath(EX, "input.gacode")
    TGLFEP = joinpath(EX, "input_singleradius_nb6.TGLFEP")   # N_BASIS=6, SCAN_N=1

    @assert isfile(GACODE) "missing $GACODE"
    @assert isfile(TGLFEP) "missing $TGLFEP"

    # Plain (grid) GPU per-combo path.
    mktempdir() do tmp
        res = run_gacode_scan_task(GACODE, TGLFEP, 1;
            out_dir=tmp, use_gpu=true, printout=false, inner=:threads, team=nothing)
        @info "container precompile workload done" scan_index=res.scan_index ir=res.ir sfmin=res.sfmin
    end

    # Hybrid batched shift-invert grid path (inner=:batched_si): bakes the
    # batched cuBLAS SI eigensolver + the collect/replay two-phase driver so
    # the first real batched_si call in a run does not JIT the GPU kernels.
    mktempdir() do tmp
        res = run_gacode_scan_task(GACODE, TGLFEP, 1;
            out_dir=tmp, use_gpu=true, printout=false, inner=:batched_si, team=nothing)
        @info "container precompile batched_si workload done" scan_index=res.scan_index ir=res.ir sfmin=res.sfmin
    end

    # solver=:ad path (critical_factor_optimize): the fast-turnaround :only /
    # :wide / :locate modes all reach this.
    mktempdir() do tmp
        res = run_gacode_scan_task(GACODE, TGLFEP, 1;
            out_dir=tmp, use_gpu=true, printout=false, inner=:threads, team=nothing, solver=:ad)
        @info "container precompile AD workload done" scan_index=res.scan_index ir=res.ir sfmin=res.sfmin
    end

    # solver=:robust_ad path (critical_factor_robust + adaptive (ky,w) refinement).
    mktempdir() do tmp
        res = run_gacode_scan_task(GACODE, TGLFEP, 1;
            out_dir=tmp, use_gpu=true, printout=false, inner=:threads, team=nothing,
            solver=:robust_ad, refine_rounds=1)
        @info "container precompile robust_ad workload done" scan_index=res.scan_index ir=res.ir sfmin=res.sfmin
    end

    # solver=:truth path (critical_factor_truth: extended log-width locate +
    # nbasis convergence).
    mktempdir() do tmp
        res = run_gacode_scan_task(GACODE, TGLFEP, 1;
            out_dir=tmp, use_gpu=true, printout=false, inner=:threads, team=nothing, solver=:truth)
        @info "container precompile truth workload done" scan_index=res.scan_index ir=res.ir sfmin=res.sfmin
    end
end
