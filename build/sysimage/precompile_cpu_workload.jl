# Precompile workload for the full CPU sysimage: run ONE nb6 radius on the CPU so that
# PackageCompiler traces and bakes in the hot per-combo path -- tjlf_LS -> tjlf_eigensolver ->
# _standard_eigenvalues_via_solve (LAPACK gesv!/geev!) and the CPU eigenvector inverse
# iteration ldiv!(lu!(zmat), v). Method specializations depend on TYPES (Float64 / ComplexF64),
# not matrix size, so nb6 (small, fast) covers the nb8/nb16/nb32 production cases.
#
# CPU analogue of precompile_gpu_workload.jl (use_gpu=false, no CUDA). The shared tjlf_LS
# compute kernel is what dominates per-radius JIT, so tracing it here removes the same cold
# compilation every distributed CPU worker would otherwise pay.

using Pkg
Pkg.activate(normpath(@__DIR__, "..", ".."))

using TJLF
using TJLFEP
using LinearAlgebra
BLAS.set_num_threads(1)

# NB: keep these as LOCALS inside a `let` (not top-level `const`) so PackageCompiler does not
# bake a `Main.GACODE` constant into the image (which would collide with task scripts that
# define their own `const GACODE = ENV[...]`). The `let` preserves the compile trace exactly.
let
    ROOT   = normpath(@__DIR__, "..", "..")
    GACODE = joinpath(ROOT, "examples", "DIIID_202017C42_500ms_v3.1", "input.gacode")
    TGLFEP = joinpath(ROOT, "examples", "DIIID_202017C42_500ms_v3.1", "input_singleradius_nb6.TGLFEP")   # N_BASIS=6, SCAN_N=1

    @assert isfile(GACODE) "missing $GACODE"
    @assert isfile(TGLFEP) "missing $TGLFEP"

    # Exercise the full CPU per-combo path once (SFmin value irrelevant; we only need the trace).
    mktempdir() do tmp
        res = run_gacode_scan_task(GACODE, TGLFEP, 1;
            out_dir=tmp, use_gpu=false, printout=false, inner=:threads, team=nothing)
        @info "precompile workload done" scan_index=res.scan_index ir=res.ir sfmin=res.sfmin
    end

    # Also bake the autodiff (solver=:ad) path: critical_factor_optimize -> Dual eigensolves +
    # IFT descent, reached via mainsub's solver=:ad branch. Without this trace the AD timing runs
    # would pay a one-time per-task JIT that unfairly inflates the wallclock-vs-grid comparison.
    mktempdir() do tmp
        res = run_gacode_scan_task(GACODE, TGLFEP, 1;
            out_dir=tmp, use_gpu=false, printout=false, inner=:threads, team=nothing, solver=:ad)
        @info "precompile AD workload done" scan_index=res.scan_index ir=res.ir sfmin=res.sfmin
    end
end
