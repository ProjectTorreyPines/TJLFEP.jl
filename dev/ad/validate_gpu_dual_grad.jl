# Validate the L1 Complex{Dual} GPU eigensolve+IFT kernel (TJLF._gpu_solve_eig_grad!)
# against the CPU Dual path and central finite differences, at N_BASIS=32 on the DIII-D
# verification case. The GPU and CPU Dual paths share the identical math (M=B⁻¹A, eigen,
# ∂λ=(R⁻¹∂M R)ᵢᵢ); only the BLAS/LAPACK-vs-CUSOLVER backend differs, so they should agree
# to ~roundoff (~1e-8 relative), and both should match FD to the usual FD truncation error.
#
# Run on a GPU node (JIT — do NOT use a sysimage built before the kernel was added):
#   module load cudatoolkit/12.9 julia/1.11.7 && export JULIA_DEPOT_PATH=$PSCRATCH/.julia
#   julia --project=. dev/ad/validate_gpu_dual_grad.jl

using CUDA
using TJLF
using TJLFEP
using Printf

const CASE = joinpath(@__DIR__, "..", "..", "examples", "DIIID_202017C42_500ms_v3.1")

const IR        = 38
const KYHAT_IN  = 0.25
const WIDTH_IN  = 1.0
const FACTOR_IN = 2.5
const N_BASIS   = 32

function gamma_at(ep0, pr, use_gpu; kyhat = KYHAT_IN, width = WIDTH_IN, factor = FACTOR_IN)
    ep = deepcopy(ep0)
    ep.KYHAT_IN  = kyhat
    ep.WIDTH_IN  = width
    ep.FACTOR_IN = factor
    gamma_dgamma_dfactor(ep, pr; use_gpu = use_gpu).gamma
end

function main()
    @printf("CUDA.functional() = %s  (%s)\n", CUDA.functional(),
            CUDA.functional() ? CUDA.name(first(CUDA.devices())) : "no GPU")
    @assert CUDA.functional() "this validation must run on a functional GPU node"
    @assert TJLF._cuda_functional() "TJLF does not see a functional CUDA runtime"

    gacode = joinpath(CASE, "input.gacode")
    tglfep = joinpath(CASE, "input.TGLFEP")
    opts, prof, _ = preprocess_gacode_inputs(gacode, tglfep)
    opts.IR        = IR
    opts.N_BASIS   = N_BASIS
    opts.KYHAT_IN  = KYHAT_IN
    opts.WIDTH_IN  = WIDTH_IN
    opts.FACTOR_IN = FACTOR_IN

    @printf("DIII-D point: IR=%d  N_BASIS=%d  kyhat=%.4f  width=%.4f  factor=%.4f\n\n",
            IR, N_BASIS, KYHAT_IN, WIDTH_IN, FACTOR_IN)

    vars = (:FACTOR_IN, :KYHAT_IN, :WIDTH_IN)
    # warm up / compile both paths
    gamma_grad(opts, prof; vars = vars, use_gpu = false)
    gamma_grad(opts, prof; vars = vars, use_gpu = true)

    gc = gamma_grad(opts, prof; vars = vars, use_gpu = false)
    gg = gamma_grad(opts, prof; vars = vars, use_gpu = true)
    nm = length(gc.gamma)

    maxval = 0.0; maxgrad = 0.0
    println("(A) GPU vs CPU Dual: γ/freq values and ∂γ/∂(factor,kyhat,width)")
    println("  mode |   γ_cpu        γ_gpu       relΔ   |  freq relΔ |  max grad relΔ over 3 cols")
    for n in 1:nm
        rv = abs(gc.gamma[n] - gg.gamma[n]) / max(abs(gc.gamma[n]), 1e-12)
        rf = abs(gc.freq[n] - gg.freq[n]) / max(abs(gc.freq[n]), 1e-12)
        rg = 0.0
        for k in 1:length(vars)
            rg = max(rg, abs(gc.dgamma[n,k] - gg.dgamma[n,k]) / max(abs(gc.dgamma[n,k]), 1e-12))
        end
        maxval = max(maxval, rv); maxgrad = max(maxgrad, rg)
        @printf("   %2d  | %12.5e %12.5e  %8.1e | %8.1e  |  %8.1e\n",
                n, gc.gamma[n], gg.gamma[n], rv, rf, rg)
    end
    @printf("\n  MAX rel Δ:  γ value = %.2e   gradient = %.2e\n", maxval, maxgrad)

    # (B) both AD paths vs finite differences for the kyhat/width columns
    col = Dict(v => k for (k, v) in enumerate(vars))
    for (label, var, base) in (("KYHAT_IN", :KYHAT_IN, KYHAT_IN), ("WIDTH_IN", :WIDTH_IN, WIDTH_IN))
        h = 1e-4
        gp = label == "KYHAT_IN" ? gamma_at(opts, prof, false; kyhat = base + h) : gamma_at(opts, prof, false; width = base + h)
        gm = label == "KYHAT_IN" ? gamma_at(opts, prof, false; kyhat = base - h) : gamma_at(opts, prof, false; width = base - h)
        println("\n(B) ∂γ/∂", label, " — GPU-AD vs FD(h=$h)  [mode: AD_gpu  FD  relΔ]")
        for n in 1:nm
            fd = (gp[n] - gm[n]) / (2h)
            ad = gg.dgamma[n, col[var]]
            @printf("   mode %2d:  %12.5e  %12.5e   %8.2e\n", n, ad, fd, abs(ad-fd)/max(abs(fd),1e-12))
        end
    end

    ok = (maxval < 1e-6) && (maxgrad < 1e-5)
    println("\nRESULT: ", ok ? "PASS" : "FAIL",
            @sprintf("  (γ relΔ=%.1e < 1e-6, grad relΔ=%.1e < 1e-5)", maxval, maxgrad))
    ok || error("GPU Dual derivatives disagree with CPU beyond tolerance")
end

main()
