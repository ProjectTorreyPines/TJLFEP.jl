# Is the narrow-width edge AE mode a resolved eigenmode or a basis artifact, and
# is the alpha_MHD^(-1/4) "Gaussian self-consistency" floor the right validity
# bound for it?
#
# Two tests, per edge radius (scan 18/19/20 = IR90/95/101), at the AD-located
# (ky, width, sfmin) read from the archived nb=32 :locate task files:
#
#   (1) N_BASIS convergence at FIXED narrow WIDTH: sweep the Hermite basis size
#       and record the leading in-AE-band growth rate gamma_AE, its frequency,
#       and the |phi|-weighted ballooning moment <theta^2>. If gamma_AE and
#       <theta^2> converge as N_BASIS grows, the narrow-width mode is a genuine,
#       resolved eigenmode (the Hermite basis at small WIDTH captures it), NOT a
#       basis artifact -> the alpha^(-1/4) floor is not a numerical necessity.
#
#   (2) Width-scale comparison: the actual mode width <theta^2> vs the two
#       candidate "expected" widths -- the thermal ideal-MHD ballooning scale
#       W_char = alpha_MHD^(-1/4), and the Alfven/geometry scale ~epsilon = r/R.
#       alpha_MHD = q^2 (R/a) sum_s beta_s (a/Ln_s + a/LT_s), beta_s = betae*as*taus.
#
#   module load julia/1.11.7 && export JULIA_DEPOT_PATH=$PSCRATCH/.julia
#   julia -t 8 --project=. dev/ad/narrow_mode_validity.jl

using TJLFEP
using TJLF
using LinearAlgebra
using Serialization
using Printf
using Plots

const CASE   = joinpath(@__DIR__, "..", "..", "examples", "DIIID_202017C42_500ms_v3.1")
const GACODE = joinpath(CASE, "input.gacode")
const TGLFEP = joinpath(CASE, "input_scan20_nb32.TGLFEP")
const AD_TASKS = joinpath(@__DIR__, "ad_threads_sfmin_nb32_ad_tasks")
const OUTDIR = @__DIR__
const USE_GPU = get(ENV, "NMV_GPU", "0") != "0"
const RHO_OF_SCAN = Dict(18=>0.89, 19=>0.94, 20=>1.00, 3=>0.11, 8=>0.37, 13=>0.63)

load_task(dir, ix) = (f = joinpath(dir, "task_$(ix).jls"); isfile(f) ? deserialize(f) : nothing)

# Single eigenvalue solve at (ir, ky, w, factor) with an explicit Hermite basis
# size nb. iflux=true also returns the eigenvectors (wavefunction) for <theta^2>.
function eval_point(opts0, prof, ir, ky, w, factor, nb; iflux::Bool=true)
    ep = deepcopy(opts0)
    ep.IR = ir; ep.KYHAT_IN = ky; ep.WIDTH_IN = w; ep.FACTOR_IN = factor
    ep.N_BASIS = nb
    ep.FREQ_AE_UPPER = -abs(prof.omegaGAM[ir])
    ep.GAMMA_THRESH = 1.0e-7; ep.GAMMA_THRESH_MAX = 1.0e-7

    it = TJLFEP.TJLF_map(ep, prof)
    it.USE_TRANSPORT_MODEL = false
    it.KYGRID_MODEL = 0
    it.NMODES = ep.NMODES
    it.NBASIS_MIN = nb
    it.NBASIS_MAX = nb
    it.NXGRID = 32
    it.WIDTH = w
    it.FIND_WIDTH = false
    it.USE_AVE_ION_GRID = false
    it.WIDTH_SPECTRUM .= w
    it.FIND_EIGEN = true
    it.RLNP_CUTOFF = 18.0
    it.BETA_LOC = 0.0
    it.DAMP_PSI = 0.0; it.DAMP_SIG = 0.0; it.WDIA_TRAPPED = 0.0
    if it.SAT_RULE == 2 || it.SAT_RULE == 3
        it.UNITS = "CGYRO"; it.XNU_MODEL = 3; it.WDIA_TRAPPED = 1.0
    end
    it.KX0_LOC = 0.0
    it.IFLUX = iflux

    local result
    try
        result = TJLF.run(it; use_gpu=USE_GPU)
    catch err
        err isa LinearAlgebra.SingularException || rethrow()
        return nothing
    end
    g = result.eigenvalue[:, 1, 1]
    f = result.eigenvalue[:, 1, 2]
    fu = -abs(prof.omegaGAM[ir])
    th2 = NaN; gband = 0.0; fband = NaN
    band = [n for n in eachindex(g) if f[n] < fu && g[n] > 1e-9]
    if !isempty(band)
        nb_lead = band[argmax(g[band])]
        gband = g[nb_lead]; fband = f[nb_lead]
        if iflux
            fw = result.field_weight_out[:, :, :, 1]
            sp = TJLF.get_sat_params(it)
            wf, angle, _, _ = TJLF.get_wavefunction(it, sp, fw)
            phi = abs.(wf[nb_lead, 1, :]); s = sum(phi)
            th2 = s > 0 ? sum((angle .^ 2) .* phi) / s : NaN
        end
    end
    return (; gband, fband, fu, th2)
end

# s-alpha thermal ballooning parameter alpha = q^2 (R/a) sum_s beta_s (a/Ln+a/LT).
function alpha_mhd(prof, ir)
    q = abs(prof.Q[ir]); Roa = prof.RMAJ[ir]; betae = prof.BETAE[ir]
    a = 0.0
    for s in 1:prof.NS
        beta_s = betae * prof.AS[ir, s] * prof.TAUS[ir, s]
        a += beta_s * (prof.RLNS[ir, s] + prof.RLTS[ir, s])
    end
    return q^2 * Roa * a
end

function main()
    @printf("threads=%d  GPU=%s\n", Threads.nthreads(), USE_GPU)
    opts0, prof, _ = preprocess_gacode_inputs(GACODE, TGLFEP)
    @printf("input N_BASIS=%d NMODES=%d\n", opts0.N_BASIS, opts0.NMODES)

    scans = parse.(Int, split(get(ENV, "NMV_SCANS", "18,19,20"), ","))
    nbs = parse.(Int, split(get(ENV, "NMV_NBASIS", "2,4,6,8,12,16,24,32"), ","))
    # Convergence is evaluated at a FIXED super-critical factor (default nominal
    # EP density 1.0, ~4-6x marginal here) so gamma_AE is robustly positive and
    # not pinned to ~0 by the marginal definition of sfmin.
    fconv = parse(Float64, get(ENV, "NMV_FACTOR", "1.0"))

    conv = Dict{Int,Any}()
    scaletab = NamedTuple[]

    for ix in scans
        at = load_task(AD_TASKS, ix); at === nothing && continue
        ir = at.ir; rho = get(RHO_OF_SCAN, ix, NaN)
        ky = at.kymark; w = at.width; sf = at.sfmin
        @printf("\n=== scan %d  IR=%d  rho=%.2f | ky=%.3f w=%.3f sf=%.4g  (conv @ factor=%.3g) ===\n", ix, ir, rho, ky, w, sf, fconv)
        @printf("  %-6s %-12s %-12s %-10s\n", "NB", "gamma_AE", "freq", "<theta^2>")
        gs = Float64[]; f0 = Float64[]; t2 = Float64[]
        for nb in nbs
            r = eval_point(opts0, prof, ir, ky, w, fconv, nb)
            g = r === nothing ? NaN : r.gband
            fr = r === nothing ? NaN : r.fband
            th2 = r === nothing ? NaN : r.th2
            push!(gs, g); push!(f0, fr); push!(t2, th2)
            @printf("  %-6d %-12.5g %-12.5g %-10.4g\n", nb, g, fr, th2)
        end
        conv[ix] = (; ir, rho, ky, w, sf, nbs, gs, f0, t2)

        # width scales
        al = alpha_mhd(prof, ir)
        wchar = al > 0 ? al^(-0.25) : NaN
        eps = prof.RMIN[ir] / prof.RMAJ[ir]
        th2_conv = t2[end]                    # <theta^2> at the largest N_BASIS
        push!(scaletab, (; ix, ir, rho, w, alpha=al, wchar, eps, th2=th2_conv,
                          gconv=gs[end], gnb6=gs[findfirst(==(6), nbs)]))
    end

    open(joinpath(OUTDIR, "narrow_mode_validity.txt"), "w") do io
        @printf(io, "N_BASIS convergence at fixed narrow WIDTH (AD-located edge AE modes)\n")
        for ix in scans
            haskey(conv, ix) || continue
            c = conv[ix]
            @printf(io, "\nscan %d IR=%d rho=%.2f ky=%.3f w=%.3f sf=%.4g\n", ix, c.ir, c.rho, c.ky, c.w, c.sf)
            @printf(io, "  %-6s %-12s %-12s %-10s\n", "NB", "gamma_AE", "freq", "<theta^2>")
            for (j, nb) in enumerate(c.nbs)
                @printf(io, "  %-6d %-12.5g %-12.5g %-10.4g\n", nb, c.gs[j], c.f0[j], c.t2[j])
            end
        end
        @printf(io, "\nWidth-scale comparison (per radius)\n")
        @printf(io, "%-4s %-5s %-5s | %-8s %-9s %-11s %-8s %-9s | %-10s\n",
            "ix","IR","rho","w_used","alpha","alpha^-1/4","eps=r/R","<th^2>","gamma_conv")
        for r in scaletab
            @printf(io, "%-4d %-5d %-5.2f | %-8.3f %-9.4g %-11.4g %-8.4g %-9.4g | %-10.4g\n",
                r.ix, r.ir, r.rho, r.w, r.alpha, r.wchar, r.eps, r.th2, r.gconv)
        end
    end
    println("\nwrote narrow_mode_validity.txt")

    # convergence plot: gamma_AE(N_BASIS)/gamma_AE(nb=32) per radius
    p = plot(title="Narrow-WIDTH edge AE: gamma_AE convergence vs N_BASIS",
             xlabel="N_BASIS", ylabel="gamma_AE / gamma_AE(N=32)", legend=:bottomright,
             size=(760, 520), left_margin=4Plots.mm, bottom_margin=4Plots.mm)
    for ix in scans
        haskey(conv, ix) || continue
        c = conv[ix]
        gref = c.gs[end]
        isfinite(gref) && gref > 0 || continue
        plot!(p, c.nbs, c.gs ./ gref; marker=:circle, linewidth=2,
              label=@sprintf("IR%d rho=%.2f w=%.2f", c.ir, c.rho, c.w))
    end
    hline!(p, [1.0]; color=:black, linestyle=:dot, label="converged")
    savefig(p, joinpath(OUTDIR, "narrow_mode_nbasis_convergence.png"))
    println("wrote narrow_mode_nbasis_convergence.png")

    serialize(joinpath(OUTDIR, "narrow_mode_validity.jls"), (; conv, scaletab))
end

main()
