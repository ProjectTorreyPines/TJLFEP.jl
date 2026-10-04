# The EP single-ky Dual path dispatches to _standard_eigenvalues_via_solve →
# _herm_eigen(B \ A). Profile that ACTUAL rule at production size, and the
# underlying Float64 eigen() calls it makes, to localize the per-solve overhead.
#
#   module load julia/1.11.7 && export JULIA_DEPOT_PATH=$PSCRATCH/.julia
#   julia --project=. dev/ad/profile_ad_eigen.jl

using TJLF
using ForwardDiff
using LinearAlgebra
using Printf

function make_dual_nonherm(n)
    Tag = typeof(ForwardDiff.Tag(:bench, Float64))
    Mr = randn(n, n); Mi = randn(n, n); Pr = randn(n, n); Pi = randn(n, n)
    [Complex(ForwardDiff.Dual{Tag}(Mr[i, j], Pr[i, j]),
             ForwardDiff.Dual{Tag}(Mi[i, j], Pi[i, j])) for i in 1:n, j in 1:n]
end

run1(f) = (f(); @elapsed f())

function profile_n(n)
    M = make_dual_nonherm(n)
    Af = map(a -> ComplexF64(ForwardDiff.value(real(a)), ForwardDiff.value(imag(a))), M)

    t_herm = run1(() -> TJLF._herm_eigen(M))
    t_eigV = run1(() -> eigen(Af))            # right eigenvectors (geev 'V')
    t_eigA = run1(() -> eigen(Af'))           # left  (second decomposition)
    t_eigN = run1(() -> eigvals(copy(Af)))    # eigenvalues only (geev 'N')

    @printf("n=%4d  _herm_eigen(Dual)=%7.2f s | eigen(A)=%6.2f  eigen(A')=%6.2f  eigvals=%6.2f s\n",
            n, t_herm, t_eigV, t_eigA, t_eigN)
end

function main()
    profile_n(200)    # warm up / small reference
    profile_n(1440)   # production size
end

main()
