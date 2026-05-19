# ─────────────────────────────────────────────────────────────────────────────
# AD-compatible dispatch for the generalized complex eigenvalue problem.
#
# TJLF's existing tjlf_ad_extensions.jl covers:
#   - eigen(Symmetric{Dual})          — real symmetric
#   - eigen(AbstractMatrix{Complex{Dual}})  — standard complex
#
# Missing (needed by EP_TJLF): eigen(A, B) where both are Complex{Dual}.
# This arises in tjlf_eigensolver._generalized_eigenvalues, which calls
# eigen(A, B) for non-ComplexF64 element types.
#
# Generalized problem:  A rᵢ = λᵢ B rᵢ,   lᵢᴴ A = λᵢ lᵢᴴ B
# B-biorthogonality (after normalisation): lᵢᴴ B rⱼ = δᵢⱼ
#
# IFT derivatives:
#   ∂λᵢ/∂p  =  lᵢᴴ (∂A/∂p − λᵢ ∂B/∂p) rᵢ
#   Sₖ[j,i] =  [Lᴴ(∂A/∂p − λᵢ ∂B/∂p)R]_{j,i} / (λᵢ − λⱼ)   (j ≠ i)
#   ∂R/∂p   =  R Sₖ
# ─────────────────────────────────────────────────────────────────────────────

import ForwardDiff
import LinearAlgebra
import LinearAlgebra.LAPACK.ggev!
import TJLF

const _GEN_DEGEN_THRESHOLD = 1e-12

function LinearAlgebra.eigen(A::AbstractMatrix{Complex{D}}, B::AbstractMatrix{Complex{D}}; kwargs...) where {D <: ForwardDiff.Dual}
    np  = ForwardDiff.npartials(D)
    Tag = ForwardDiff.tagtype(D)

    # Extract Float64 value matrices
    Af = map(a -> Complex{Float64}(ForwardDiff.value(real(a)), ForwardDiff.value(imag(a))), A)
    Bf = map(b -> Complex{Float64}(ForwardDiff.value(real(b)), ForwardDiff.value(imag(b))), B)

    # LAPACK ggev! returns left (vl) and right (vr) eigenvectors.
    # It overwrites its inputs, so pass copies.
    (alpha, beta, L, R) = ggev!('V', 'V', copy(Af), copy(Bf))
    λf = alpha ./ beta    # Vector{ComplexF64}
    n  = length(λf)

    # B-biorthogonal normalisation: rescale L so that lᵢᴴ B rᵢ = 1
    BfR = Bf * R
    for i in 1:n
        s = LinearAlgebra.dot(L[:, i], BfR[:, i])
        if abs(s) > 1e-30
            L[:, i] ./= conj(s)
        end
    end

    # Precompute derivatives for all partial directions
    dλ_re = zeros(n, np)
    dλ_im = zeros(n, np)
    dR_re = zeros(n, n, np)
    dR_im = zeros(n, n, np)

    for k in 1:np
        dAk = map(a -> Complex{Float64}(ForwardDiff.partials(real(a), k),
                                        ForwardDiff.partials(imag(a), k)), A)
        dBk = map(b -> Complex{Float64}(ForwardDiff.partials(real(b), k),
                                        ForwardDiff.partials(imag(b), k)), B)

        # Ck[j,i] = [Lᴴ dAk R - λᵢ Lᴴ dBk R]_{j,i}
        # λᵢ varies per column → broadcast transpose(λf) across rows
        LH_dAk_R = L' * (dAk * R)
        LH_dBk_R = L' * (dBk * R)
        Ck = LH_dAk_R .- LH_dBk_R .* transpose(λf)

        for i in 1:n
            dλ_re[i, k] = real(Ck[i, i])
            dλ_im[i, k] = imag(Ck[i, i])
        end

        Sk = zeros(ComplexF64, n, n)
        for i in 1:n, j in 1:n
            if i != j
                gap = λf[i] - λf[j]
                if abs(gap) > _GEN_DEGEN_THRESHOLD
                    Sk[j, i] = Ck[j, i] / gap
                end
            end
        end
        dRk = R * Sk
        dR_re[:, :, k] .= real.(dRk)
        dR_im[:, :, k] .= imag.(dRk)
    end

    # Construct Dual eigenvalues
    λ = map(1:n) do i
        Complex(ForwardDiff.Dual{Tag}(real(λf[i]), ntuple(k -> dλ_re[i, k], Val(np))...),
                ForwardDiff.Dual{Tag}(imag(λf[i]), ntuple(k -> dλ_im[i, k], Val(np))...))
    end

    # Construct Dual right eigenvectors
    Rd = Matrix{Complex{D}}(undef, n, n)
    for j in 1:n, i in 1:n
        Rd[i, j] = Complex(ForwardDiff.Dual{Tag}(real(R[i, j]), ntuple(k -> dR_re[i, j, k], Val(np))...),
                           ForwardDiff.Dual{Tag}(imag(R[i, j]), ntuple(k -> dR_im[i, j, k], Val(np))...))
    end

    return (values = λ, vectors = Rd)
end

# ── AD-safe run for non-Float64 element types (Dual, etc.) ───────────────────
#
# TJLF.run's single-ky branch hardcodes zeros(Float64,...)/zeros(ComplexF64,...)
# for result arrays, which fails for Dual. This function is identical to that
# branch but allocates with zeros(T,...). All internal TJLF calls are unchanged.
# Called from TJLFEP_ky via:  T === Float64 ? TJLF.run(x) : _tjlf_run_dual(x, T)
function _tjlf_run_dual(inputTJLF, ::Type{T}; use_gpu::Bool = false) where {T}
    TJLF.checkInput(inputTJLF)
    outputHermite = TJLF.gauss_hermite(inputTJLF)
    satParams     = TJLF.get_sat_params(inputTJLF)
    inputTJLF.KY_SPECTRUM .= TJLF.get_ky_spectrum(inputTJLF, satParams.grad_r0)

    ns     = inputTJLF.NS
    nmodes = inputTJLF.NMODES
    nbasis = inputTJLF.NBASIS_MAX

    nmodes_out, gamma_out, freq_out,
        particle_QL_out, energy_QL_out, stress_tor_QL_out, stress_par_QL_out, exchange_QL_out,
        _ft, field_weight_out_3d, _phi = TJLF.tjlf_LS(inputTJLF, satParams, outputHermite,
                                                       inputTJLF.KY, nbasis, inputTJLF.VEXB_SHEAR, 1;
                                                       use_gpu=use_gpu)

    eigenvalue       = zeros(T,          nmodes, 1, 2)
    QL_weights       = zeros(T,          3, ns, nmodes, 1, 5)
    field_weight_out = zeros(Complex{T}, 3, nbasis, nmodes, 1)

    eigenvalue[:, 1, 1]                      .= gamma_out
    eigenvalue[:, 1, 2]                      .= freq_out
    QL_weights[:, :, 1:nmodes_out, 1, 1]    .= particle_QL_out[:, :, 1:nmodes_out]
    QL_weights[:, :, 1:nmodes_out, 1, 2]    .= energy_QL_out[:, :, 1:nmodes_out]
    QL_weights[:, :, 1:nmodes_out, 1, 3]    .= stress_tor_QL_out[:, :, 1:nmodes_out]
    QL_weights[:, :, 1:nmodes_out, 1, 4]    .= stress_par_QL_out[:, :, 1:nmodes_out]
    QL_weights[:, :, 1:nmodes_out, 1, 5]    .= exchange_QL_out[:, :, 1:nmodes_out]
    field_weight_out[:, :, 1:nmodes_out, 1] .= field_weight_out_3d[:, :, 1:nmodes_out]

    QL_flux_out   = Array{T}(undef, 0, 0, 0)
    flux_spectrum = Array{T}(undef, 0, 0, 0)
    return (QL_weights=QL_weights, eigenvalue=eigenvalue, QL_flux_out=QL_flux_out,
            flux_spectrum=flux_spectrum, field_weight_out=field_weight_out)
end

# ── Strip Dual partials from InputTJLF → InputTJLF{Float64} ─────────────────
# Used before calling get_wavefunction, which has Float64-typed geometry arrays
# internally. The wavefunction shape is independent of the differentiated
# parameters, so extracting .value is exact (no information lost).
#
# NOTE: `src` may be either TJLFEP.InputTJLF or TJLF.InputTJLF; no type
# annotation so both are accepted. Runtime field copy works for either.
function _to_float64_input(src)
    dst = TJLF.InputTJLF{Float64}(src.NS, src.NKY)
    for fn in fieldnames(typeof(src))
        v = getfield(src, fn)
        if v isa ForwardDiff.Dual
            setfield!(dst, fn, ForwardDiff.value(v))
        elseif v isa AbstractArray && eltype(v) <: ForwardDiff.Dual
            setfield!(dst, fn, ForwardDiff.value.(v))
        elseif v isa AbstractArray && eltype(v) <: Complex && real(eltype(v)) <: ForwardDiff.Dual
            setfield!(dst, fn, map(x -> ComplexF64(ForwardDiff.value(real(x)), ForwardDiff.value(imag(x))), v))
        else
            try; setfield!(dst, fn, v); catch; end
        end
    end
    return dst
end
