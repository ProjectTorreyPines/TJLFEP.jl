# Isolate TJLF's Dual `_herm_eigen` eigenvalue-derivative rule from the TGLF
# physics, to confirm whether the dγ/dfactor mismatch comes from the eigen rule
# itself (non-normal, clustered spectrum) rather than the matrix assembly.
#
# We build a parameterized complex matrix M(p), compute dλ/dp two ways:
#   - AD : TJLF._herm_eigen on a Dual matrix
#   - FD : central difference of eigenvalues of the Float64 matrix
# and report the worst per-eigenvalue mismatch for a few matrix regimes.
#
#   module load julia/1.11.7 && export JULIA_DEPOT_PATH=$PSCRATCH/.julia
#   julia --project=. dev/ad/check_herm_eigen_ad.jl

using TJLF
using ForwardDiff
using LinearAlgebra
using Printf
using Random

# Eigenvalues of M(p) sorted by a stable key (descending real part) so AD/FD
# refer to the same eigenvalue. Returns Vector{Complex}.
sorteig(vals) = sort(vals; by = z -> (-real(z), -imag(z)))

function ad_deriv(Mfun, p0::Float64)
    Tag = typeof(ForwardDiff.Tag(Mfun, Float64))
    D   = ForwardDiff.Dual{Tag, Float64, 1}
    Md  = Mfun(ForwardDiff.Dual{Tag}(p0, 1.0))   # Matrix{Complex{D}}
    vals = TJLF._herm_eigen(Md).values
    vals = sorteig(vals)
    return [Complex(ForwardDiff.partials(real(z), 1), ForwardDiff.partials(imag(z), 1)) for z in vals]
end

function fd_deriv(Mfun, p0::Float64; h = 1e-6)
    vp = sorteig(eigen(Mfun(p0 + h)).values)
    vm = sorteig(eigen(Mfun(p0 - h)).values)
    return (vp .- vm) ./ (2h)
end

function report(name, Mfun, p0)
    ad = ad_deriv(Mfun, p0)
    fd = fd_deriv(Mfun, p0)
    n  = length(ad)
    errs = abs.(ad .- fd)
    rels = errs ./ max.(abs.(fd), 1e-12)
    worst = argmax(rels)
    @printf("%-28s n=%-4d  max|Δ|=%.2e  max rel=%.2e  (worst eig: AD=%s FD=%s)\n",
            name, n, maximum(errs), maximum(rels),
            string(round(ad[worst]; digits=4)), string(round(fd[worst]; digits=4)))
end

function main()
    Random.seed!(1)

    # 1) Small, well-separated non-normal matrix: eigen rule should be ~exact.
    let A = randn(ComplexF64, 6, 6), Bp = randn(ComplexF64, 6, 6)
        report("well-separated 6x6", p -> A .+ p .* Bp, 0.3)
    end

    # 2) Matrix with a near-degenerate pair (clustered eigenvalues).
    let
        base = diagm(ComplexF64[1.0, 1.0 + 1e-3, 2.0, 3.0, -1.0])
        N    = 1e-2 .* randn(ComplexF64, 5, 5)   # non-normal coupling
        Bp   = randn(ComplexF64, 5, 5)
        report("near-degenerate 5x5", p -> base .+ N .+ p .* Bp, 0.1)
    end

    # 3) Large non-normal matrix sized like the EP pencil after B\A.
    let n = 200
        A  = randn(ComplexF64, n, n)
        Bp = randn(ComplexF64, n, n)
        report("large non-normal 200", p -> A .+ p .* Bp, 0.05)
    end

    # 4) Same large matrix but mildly non-normal + clustered (closer to TGLF).
    let n = 200
        D0 = diagm(ComplexF64[(0.5 * (i % 7)) + 0.01im * i for i in 1:n])  # many repeats
        N  = 5e-2 .* randn(ComplexF64, n, n)
        Bp = randn(ComplexF64, n, n)
        report("clustered non-normal 200", p -> D0 .+ N .+ p .* Bp, 0.05)
    end
end

main()
