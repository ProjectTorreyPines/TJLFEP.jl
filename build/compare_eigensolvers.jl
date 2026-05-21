# Compare TJLFEP CPU (ggev) vs GPU (gesv+Xgeev) eigenvalue paths.
using LinearAlgebra
import LinearAlgebra.LAPACK: ggev!, gesv!, geev!
using CUDA
ENV["TJLFEP_FILE_ONLY"] = "1"
using TJLFEP
using TJLF

function λ_ggev(A, B)
    (α, β, _, _) = ggev!('N', 'N', copy(A), copy(B))
    return α ./ β
end

function λ_gesv_geev(A, B)
    Ac, Bc = copy(A), copy(B)
    gesv!(Bc, Ac)
    return geev!('N', 'N', Ac)[1]
end

# Random well-conditioned test
n = 120
A = randn(ComplexF64, n, n)
B = randn(ComplexF64, n, n)
B = B * B' + 10.0 * I
λ1 = λ_ggev(A, B)
λ2 = λ_gesv_geev(A, B)
λ3 = TJLF._gpu_solve!(copy(A), copy(B))
sortmatch(a, b) = begin
    as = sort(a; by=x -> (real(x), imag(x)))
    bs = sort(b; by=x -> (real(x), imag(x)))
    maximum(abs.(as .- bs) ./ max.(abs.(as), 1e-30))
end

println("random ", n, "x", n, " (sorted eigenvalues):")
println("  ggev vs gesv+geev max rel: ", sortmatch(λ1, λ2))
println("  ggev vs gpu max rel:       ", sortmatch(λ1, λ3))

# Unstable-mode growth proxy (max positive real part, as in tjlf_LS)
println("  max real(ggev): ", maximum(real.(λ1)))
println("  max real(gesv): ", maximum(real.(λ2)))
println("  max real(gpu):  ", maximum(real.(λ3)))
