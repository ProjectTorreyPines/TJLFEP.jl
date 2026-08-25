# Build the full CPU sysimage: TJLF + TJLFEP (file-only), with the CPU eigensolve path
# precompiled via precompile_cpu_workload.jl. This bakes in the per-radius JIT that cold
# distributed CPU workers otherwise pay on every spawn. CUDA and the FUSE/IMAS stack are
# excluded by construction: the TJLFEP project manifest never resolves the IMAS/GACODE/
# TurbulentTransport weak deps, so TJLFEPIMASExt stays dormant.
#
# CPU analogue of build_gpu_sysimage_generic.jl, baking TJLF + TJLFEP for the file-based path. This
# DOES bake TJLF (it was previously left out so TJLF source could be iterated during
# development without rebuilding the image). Run on a CPU node.

using Pkg
Pkg.activate(normpath(@__DIR__, "..", ".."))
Pkg.instantiate()
using PackageCompiler

create_sysimage(
    [:TJLF, :TJLFEP];
    sysimage_path = get(ENV, "CPU_SYSIMAGE_OUT", normpath(@__DIR__, "..", "TJLFEP_cpu_sysimage.so")),
    precompile_execution_file = normpath(@__DIR__, "precompile_cpu_workload.jl"),
    cpu_target = PackageCompiler.default_app_cpu_target(),
)
