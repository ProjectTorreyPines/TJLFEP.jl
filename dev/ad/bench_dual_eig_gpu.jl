# Benchmark the Complex{Dual} standard-eigenvalue solve (the :ad hot kernel) on GPU vs CPU,
# at the DIII-D (iur=1440) and UCP (iur=2400) nb32 sizes, for np=1 (gamma_dgamma_dfactor) and
# np=3 (gamma_grad). This isolates whether fix #1 (cuSOLVER Xgeev + IFT) actually accelerates the
# reactor-size Dual eigensolve, or whether cuSOLVER Xgeev is itself the bottleneck.
#   USE_GPU=1 julia --project=<root> --sysimage=<gpu.so> ad/bench_dual_eig_gpu.jl
using Pkg
Pkg.activate(normpath(@__DIR__, "..", ".."))
using CUDA, TJLF, ForwardDiff, LinearAlgebra, Printf
BLAS.set_num_threads(1)

println("CUDA functional: ", CUDA.functional(), "  device: ", CUDA.functional() ? CUDA.name(first(CUDA.devices())) : "n/a")
println("_cuda_functional(): ", TJLF._cuda_functional())
flush(stdout)

function make_dual_pencil(n::Int, np::Int)
    Tag = ForwardDiff.Tag(:bench, Float64)
    D = ForwardDiff.Dual{typeof(Tag), Float64, np}
    mk() = begin
        M = Matrix{Complex{D}}(undef, n, n)
        for j in 1:n, i in 1:n
            re = ForwardDiff.Dual{typeof(Tag)}(randn(), ntuple(_ -> randn(), np)...)
            im = ForwardDiff.Dual{typeof(Tag)}(randn(), ntuple(_ -> randn(), np)...)
            M[i, j] = Complex(re, im)
        end
        M
    end
    A = mk()
    # B diagonally dominant so B\A is well-conditioned
    B = mk()
    for i in 1:n
        B[i, i] += Complex(ForwardDiff.Dual{typeof(Tag)}(Float64(n), ntuple(_ -> 0.0, np)...),
                           ForwardDiff.Dual{typeof(Tag)}(0.0, ntuple(_ -> 0.0, np)...))
    end
    return A, B
end

function timeit(f; reps=3)
    f()  # warm
    best = Inf
    for _ in 1:reps
        t = @elapsed f()
        best = min(best, t)
    end
    best
end

for n in (1440, 2400), np in (1, 3)
    A, B = make_dual_pencil(n, np)
    tcpu = timeit(() -> TJLF._standard_eigenvalues_via_solve(A, B; use_gpu=false))
    tgpu = try
        timeit(() -> TJLF._standard_eigenvalues_via_solve(A, B; use_gpu=true))
    catch e
        println("  GPU path errored: ", sprint(showerror, e)); NaN
    end
    # correctness: compare eigenvalue values (sorted) + a partial
    vc = TJLF._standard_eigenvalues_via_solve(A, B; use_gpu=false)
    vg = TJLF._standard_eigenvalues_via_solve(A, B; use_gpu=true)
    valc = sort(ForwardDiff.value.(real.(vc)))
    valg = sort(ForwardDiff.value.(real.(vg)))
    maxval = maximum(abs, valc .- valg)
    @printf("n=%-5d np=%d  CPU=%7.3fs  GPU=%7.3fs  speedup=%5.1fx  max|Δeig_val|=%.2e\n",
            n, np, tcpu, tgpu, tcpu/tgpu, maxval)
    flush(stdout)
end
println("=== bench done ===")
