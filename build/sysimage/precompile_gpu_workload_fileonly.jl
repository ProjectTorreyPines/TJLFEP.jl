# Precompile workload for the FILE-ONLY GPU sysimage (CUDA + TJLF + TJLFEP, no FUSE/IMAS).
# Same GPU compute trace as precompile_gpu_workload_generic.jl (one nb6 radius per solver on the
# GPU -> tjlf_LS eigenvalue + eigenvector CUSOLVER paths), but loads TJLFEP STANDALONE: only
# `using CUDA, TJLF, TJLFEP`, never `using FUSE/IMAS`, so the TJLFEPIMASExt extension stays
# dormant and the IMAS/FUSE stack is NOT baked. This is what a TGLF-EP user running the
# file-based scan path (run_gacode_scan_task) gets -- a leaner, faster-loading image than the
# generic one, with the same solver compute paths baked in.
#
# Built against the registry-resolved "lean" build env (see build_gpu_sysimage_fileonly.jl);
# CUDA/TJLF resolve as direct deps. Do NOT Pkg.activate a different project here -- that
# would fight create_sysimage's build project.
using CUDA
using TJLF
using TJLFEP
using LinearAlgebra
BLAS.set_num_threads(1)

if !CUDA.functional()
    error("precompile_gpu_workload_fileonly.jl must run on a GPU node (CUDA.functional() == false)")
end
@info "file-only precompile workload GPU" name=CUDA.name(first(CUDA.devices()))

# NB: keep these as LOCALS inside a `let` (not top-level `const`) -- a baked `const GACODE`
# in Main would collide with the MPS task script's own `const GACODE = ENV[...]` at runtime.
let
    # Example inputs come from the installed TJLFEP package (registry copy), not a repo
    # checkout -- the traced code and its inputs stay from the same source.
    EX     = joinpath(pkgdir(TJLFEP), "examples", "DIIID_202017C42_500ms_v3.1")
    GACODE = joinpath(EX, "input.gacode")
    TGLFEP = joinpath(EX, "input_singleradius_nb6.TGLFEP")   # N_BASIS=6, SCAN_N=1

    @assert isfile(GACODE) "missing $GACODE"
    @assert isfile(TGLFEP) "missing $TGLFEP"

    # Plain (grid) GPU per-combo path.
    mktempdir() do tmp
        res = run_gacode_scan_task(GACODE, TGLFEP, 1;
            out_dir=tmp, use_gpu=true, printout=false, inner=:threads, team=nothing)
        @info "file-only precompile workload done" scan_index=res.scan_index ir=res.ir sfmin=res.sfmin
    end

    # Hybrid batched shift-invert grid path (inner=:batched_si): bakes the batched cuBLAS SI
    # eigensolver + the collect/replay two-phase driver so the first real batched_si call in a run
    # does not JIT the GPU kernels.
    mktempdir() do tmp
        res = run_gacode_scan_task(GACODE, TGLFEP, 1;
            out_dir=tmp, use_gpu=true, printout=false, inner=:batched_si, team=nothing)
        @info "file-only precompile batched_si workload done" scan_index=res.scan_index ir=res.ir sfmin=res.sfmin
    end

    # solver=:ad path (critical_factor_optimize): the fast-turnaround :only / :wide / :locate
    # modes all reach this, so baking it keeps the AD runs from JIT-ing per task.
    mktempdir() do tmp
        res = run_gacode_scan_task(GACODE, TGLFEP, 1;
            out_dir=tmp, use_gpu=true, printout=false, inner=:threads, team=nothing, solver=:ad)
        @info "file-only precompile AD workload done" scan_index=res.scan_index ir=res.ir sfmin=res.sfmin
    end

    # solver=:robust_ad path (critical_factor_robust + adaptive (ky,w) refinement).
    mktempdir() do tmp
        res = run_gacode_scan_task(GACODE, TGLFEP, 1;
            out_dir=tmp, use_gpu=true, printout=false, inner=:threads, team=nothing,
            solver=:robust_ad, refine_rounds=1)
        @info "file-only precompile robust_ad workload done" scan_index=res.scan_index ir=res.ir sfmin=res.sfmin
    end

    # solver=:truth path (critical_factor_truth: extended log-width locate + nbasis convergence).
    mktempdir() do tmp
        res = run_gacode_scan_task(GACODE, TGLFEP, 1;
            out_dir=tmp, use_gpu=true, printout=false, inner=:threads, team=nothing, solver=:truth)
        @info "file-only precompile truth workload done" scan_index=res.scan_index ir=res.ir sfmin=res.sfmin
    end
end
