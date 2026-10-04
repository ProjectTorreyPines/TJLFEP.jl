# Edge-vs-core AE mode physics in the DIII-D 202017C42 benchmark.
#
# Where TGLF-EP Fortran (and the w>=1 grid box) and TJLFEP disagree, the
# width-extended AD solver locates *narrow*-width modes (w ~ 0.1-0.3) that the
# canonical WIDTH_MIN=1 box cannot represent. We have previously confirmed these
# are AE modes (frequency in the AE band). This script looks past the mode
# *identification* at the mode *structure / physics* that the ALPHA model may be
# sensitive to:
#
#   (1) eigenfunction phi(theta) (ballooning structure) of the binding AE mode,
#       at the AD-located (ky, width) vs the wide grid (ky, width);
#   (2) the growth-rate "bump" gamma_AE(scale factor): onset + stiffness;
#   (3) scalar mode metrics: ballooning width <theta^2>, peak |phi| location,
#       tearing/parity, electromagnetic ratio |A_par|/|phi|, frequency vs the
#       AE-band upper bound.
#
# The binding (ky, width, sfmin) per radius are read from the archived nb=32 task
# files (grid = w>=1 box; ad = width-extended :locate), so no re-solve is needed.
#
#   module load julia/1.11.7 && export JULIA_DEPOT_PATH=$PSCRATCH/.julia
#   julia -t 8 --project=. dev/ad/edge_vs_core_modes.jl
#
# Optional env:
#   EVC_SCANS="3,8,13,18,19,20"   # scan indices to analyze (default these 6)
#   EVC_NSWEEP=21                 # factor-sweep points
#   EVC_GPU=1                     # use GPU eigensolves

using TJLFEP
using TJLF
using LinearAlgebra
using Serialization
using Printf
using Plots

const CASE   = joinpath(@__DIR__, "..", "..", "examples", "DIIID_202017C42_500ms_v3.1")
const GACODE = joinpath(CASE, "input.gacode")
const TGLFEP = joinpath(CASE, "input_scan20_nb32.TGLFEP")
const GRID_TASKS = joinpath(@__DIR__, "ad_threads_sfmin_nb32_tasks")      # w>=1 box (Fortran-like)
const AD_TASKS   = joinpath(@__DIR__, "ad_threads_sfmin_nb32_ad_tasks")   # width-extended :locate
const OUTDIR = @__DIR__
const USE_GPU = get(ENV, "EVC_GPU", "0") != "0"

# scan index -> rho (documented SCAN_N=20 grid, IRS=2, NR=101)
const RHO_OF_SCAN = Dict(1=>0.01,2=>0.06,3=>0.11,4=>0.16,5=>0.21,6=>0.27,7=>0.32,
    8=>0.37,9=>0.42,10=>0.47,11=>0.53,12=>0.58,13=>0.63,14=>0.68,15=>0.73,
    16=>0.79,17=>0.84,18=>0.89,19=>0.94,20=>1.00)

load_task(dir, ix) = (f = joinpath(dir, "task_$(ix).jls"); isfile(f) ? deserialize(f) : nothing)

# Replicate TJLFEP_ky's inputTJLF configuration (tjlfep_ky.jl lines 38-77).
# iflux=true populates field_weight_out (eigenvectors) for the wavefunction; the
# Float64 eigenvalue solve (iflux=false) is ~1-2s vs ~33s for the AD/Dual path,
# so the factor sweep uses iflux=false (plain gamma, no derivative).
function eval_point(opts0, prof, ir, ky, w, factor; iflux::Bool = true)
    ep = deepcopy(opts0)
    ep.IR = ir; ep.KYHAT_IN = ky; ep.WIDTH_IN = w; ep.FACTOR_IN = factor
    fu = -abs(prof.omegaGAM[ir])
    ep.FREQ_AE_UPPER = fu
    ep.GAMMA_THRESH = 1.0e-7; ep.GAMMA_THRESH_MAX = 1.0e-7

    inputTJLF = TJLFEP.TJLF_map(ep, prof)
    inputTJLF.USE_TRANSPORT_MODEL = false
    inputTJLF.KYGRID_MODEL = 0
    inputTJLF.NMODES = ep.NMODES
    inputTJLF.NBASIS_MIN = ep.N_BASIS
    inputTJLF.NBASIS_MAX = ep.N_BASIS
    inputTJLF.NXGRID = 32
    inputTJLF.WIDTH = w
    inputTJLF.FIND_WIDTH = false
    inputTJLF.USE_AVE_ION_GRID = false
    inputTJLF.WIDTH_SPECTRUM .= w
    inputTJLF.FIND_EIGEN = true
    inputTJLF.RLNP_CUTOFF = 18.0
    inputTJLF.BETA_LOC = 0.0
    inputTJLF.DAMP_PSI = 0.0
    inputTJLF.DAMP_SIG = 0.0
    inputTJLF.WDIA_TRAPPED = 0.0
    if inputTJLF.SAT_RULE == 2 || inputTJLF.SAT_RULE == 3
        inputTJLF.UNITS = "CGYRO"; inputTJLF.XNU_MODEL = 3; inputTJLF.WDIA_TRAPPED = 1.0
    end
    inputTJLF.KX0_LOC = 0.0
    inputTJLF.IFLUX = iflux   # eigenvectors (wavefunction) only when needed

    local result
    try
        result = TJLF.run(inputTJLF; use_gpu = USE_GPU)
    catch err
        err isa LinearAlgebra.SingularException || rethrow()
        return nothing
    end
    g  = result.eigenvalue[:, 1, 1]
    f  = result.eigenvalue[:, 1, 2]
    wf = nothing; angle = nothing; nmout = 0
    if iflux
        fw = result.field_weight_out[:, :, :, 1]
        sp = TJLF.get_sat_params(inputTJLF)
        wf, angle, nplot, nmout = TJLF.get_wavefunction(inputTJLF, sp, fw)
    end
    return (; g, f, fu, wf, angle, nmout)
end

# AE-band leading growth rate: max gamma among modes with freq < AE-upper bound.
function gamma_ae_lead(g, f, fu)
    band = [n for n in eachindex(g) if f[n] < fu]
    isempty(band) && return 0.0
    return max(0.0, maximum(g[band]))
end

# Pick the binding AE mode: in-band (freq < fu) with the largest growth rate.
function binding_mode(res)
    NM = length(res.g)
    inband = [n for n in 1:NM if res.f[n] < res.fu && res.g[n] > 1e-9]
    isempty(inband) && return argmax(res.g)
    return inband[argmax(res.g[inband])]
end

# Ballooning second moment <theta^2> (|phi|-weighted) and peak |phi| location.
function ballooning_metrics(res, n)
    phi = abs.(res.wf[n, 1, :])
    th  = res.angle
    s   = sum(phi)
    th2 = s > 0 ? sum((th .^ 2) .* phi) / s : NaN
    pk  = th[argmax(phi)]
    # electromagnetic ratio |A_par|/|phi| at the peak (field index 2)
    apar = abs.(res.wf[n, 2, :])
    emrat = maximum(phi) > 0 ? maximum(apar) / maximum(phi) : NaN
    # tearing/parity metric: max |phi(theta) - phi(-theta)| / max|phi|
    L = length(phi)
    tear = maximum(abs.(res.wf[n,1,:] .- res.wf[n,1,end:-1:1])) / (maximum(phi) + 1e-12)
    return (; th2, pk, emrat, tear)
end

# gamma_AE(factor) at fixed (ir, ky, w): leading in-AE-band growth rate, via the
# fast Float64 eigenvalue solve (no AD derivative).
function gamma_bump(opts0, prof, ir, ky, w, sfmin; nsweep = 21)
    facs = exp10.(range(log10(max(sfmin, 1e-4)) - 1.0, log10(max(sfmin, 1e-4)) + 1.0; length = nsweep))
    gam  = similar(facs)
    for (i, fc) in enumerate(facs)
        r = eval_point(opts0, prof, ir, ky, w, fc; iflux = false)
        gam[i] = r === nothing ? 0.0 : gamma_ae_lead(r.g, r.f, r.fu)
    end
    return facs, gam
end

function main()
    @printf("threads=%d  GPU=%s\n", Threads.nthreads(), USE_GPU)
    scans = parse.(Int, split(get(ENV, "EVC_SCANS", "3,8,13,18,19,20"), ","))
    nsweep = parse(Int, get(ENV, "EVC_NSWEEP", "21"))

    opts0, prof, _ = preprocess_gacode_inputs(GACODE, TGLFEP)
    @printf("N_BASIS=%d  NMODES=%d\n", opts0.N_BASIS, opts0.NMODES)

    rows = NamedTuple[]
    wfs  = Dict{Int,Any}()   # scan_ix => (ad_res, ad_n, grid_res, grid_n)
    bumps = Dict{Int,Any}()  # scan_ix => (facs, gam_ad, gam_grid, sf_ad, sf_grid)

    for ix in scans
        gt = load_task(GRID_TASKS, ix); at = load_task(AD_TASKS, ix)
        at === nothing && continue
        ir = at.ir; rho = get(RHO_OF_SCAN, ix, NaN)
        @printf("\n=== scan %d  IR=%d  rho=%.2f ===\n", ix, ir, rho)

        # AD-located narrow mode
        ad_res = eval_point(opts0, prof, ir, at.kymark, at.width, at.sfmin)
        ad_n = ad_res === nothing ? 0 : binding_mode(ad_res)
        adm = ad_res === nothing ? (; th2=NaN, pk=NaN, emrat=NaN, tear=NaN) : ballooning_metrics(ad_res, ad_n)
        ad_g = ad_res === nothing ? NaN : ad_res.g[ad_n]
        ad_f = ad_res === nothing ? NaN : ad_res.f[ad_n]
        ad_fu = ad_res === nothing ? NaN : ad_res.fu

        # Wide grid mode (where the w>=1 box found one: width>0, sfmin<cap)
        grid_ok = gt !== nothing && gt.width > 0 && isfinite(gt.kymark) && gt.sfmin < 9.999
        grid_res = grid_ok ? eval_point(opts0, prof, ir, gt.kymark, gt.width, gt.sfmin) : nothing
        grid_n = grid_res === nothing ? 0 : binding_mode(grid_res)
        gm = grid_res === nothing ? (; th2=NaN, pk=NaN, emrat=NaN, tear=NaN) : ballooning_metrics(grid_res, grid_n)

        wfs[ix] = (ad_res, ad_n, grid_res, grid_n)

        # gamma(factor) bumps
        f_ad, g_ad = gamma_bump(opts0, prof, ir, at.kymark, at.width, at.sfmin; nsweep)
        if grid_ok
            f_gr, g_gr = gamma_bump(opts0, prof, ir, gt.kymark, gt.width, gt.sfmin; nsweep)
        else
            f_gr, g_gr = f_ad, fill(NaN, length(f_ad))
        end
        bumps[ix] = (f_ad, g_ad, f_gr, g_gr, at.sfmin, grid_ok ? gt.sfmin : NaN)

        push!(rows, (; ix, ir, rho,
            ad_w = at.width, ad_ky = at.kymark, ad_sf = at.sfmin,
            ad_freq = ad_f, ad_fu = ad_fu, ad_gam = ad_g,
            ad_th2 = adm.th2, ad_pk = adm.pk, ad_em = adm.emrat, ad_tear = adm.tear,
            gr_w = grid_ok ? gt.width : NaN, gr_ky = grid_ok ? gt.kymark : NaN,
            gr_sf = grid_ok ? gt.sfmin : NaN, gr_th2 = gm.th2, gr_em = gm.emrat))
        @printf("  AD:   ky=%.3f w=%.3f sf=%.4g | freq=%.4g (AEupper=%.4g) gam=%.3g <th^2>=%.3g pk=%.2f emrat=%.3g tear=%.2g\n",
            at.kymark, at.width, at.sfmin, ad_f, ad_fu, ad_g, adm.th2, adm.pk, adm.emrat, adm.tear)
        if grid_ok
            @printf("  grid: ky=%.3f w=%.3f sf=%.4g | <th^2>=%.3g emrat=%.3g\n",
                gt.kymark, gt.width, gt.sfmin, gm.th2, gm.emrat)
        else
            @printf("  grid: w>=1 box empty / capped (no wide mode)\n")
        end
    end

    serialize(joinpath(OUTDIR, "edge_vs_core_modes.jls"), (; rows, wfs, bumps))

    # ----- summary table -----
    open(joinpath(OUTDIR, "edge_vs_core_modes.txt"), "w") do io
        @printf(io, "%-4s %-5s %-5s | %-7s %-7s %-9s %-9s %-8s %-9s %-7s | %-7s %-9s %-9s\n",
            "ix","IR","rho","AD_w","AD_ky","AD_sf","AD_freq","AD_<th2>","AD_emrat","AD_tear","gr_w","gr_sf","gr_<th2>")
        for r in rows
            @printf(io, "%-4d %-5d %-5.2f | %-7.3f %-7.3f %-9.4g %-9.4g %-8.3g %-9.3g %-7.2g | %-7.3f %-9.4g %-9.3g\n",
                r.ix, r.ir, r.rho, r.ad_w, r.ad_ky, r.ad_sf, r.ad_freq, r.ad_th2, r.ad_em, r.ad_tear,
                r.gr_w, r.gr_sf, r.gr_th2)
        end
    end
    println("\nwrote edge_vs_core_modes.txt / .jls")

    make_plots(rows, wfs, bumps)
end

function make_plots(rows, wfs, bumps)
    # classify core vs edge by rho
    iscore(r) = r.rho <= 0.45
    # ---- Panel 1: eigenfunctions |phi(theta)| (AD narrow modes), core vs edge
    p1 = plot(title = "AD-located AE eigenfunction |phi(theta)|", xlabel = "ballooning angle theta/pi",
              ylabel = "|phi| (norm)", legend = :topright, size = (760, 520))
    for r in rows
        ent = get(wfs, r.ix, nothing); ent === nothing && continue
        ad_res, ad_n = ent[1], ent[2]
        ad_res === nothing && continue
        phi = abs.(ad_res.wf[ad_n, 1, :]); phi ./= (maximum(phi) + 1e-30)
        th = ad_res.angle ./ pi
        lab = @sprintf("IR%d rho=%.2f w=%.2f%s", r.ir, r.rho, r.ad_w, iscore(r) ? " (core)" : " (edge)")
        plot!(p1, th, phi; label = lab, linewidth = 2,
              linestyle = iscore(r) ? :solid : :dash)
    end
    savefig(p1, joinpath(OUTDIR, "edge_vs_core_eigenfunction.png"))

    # ---- Panel 2: gamma(factor) bumps, normalized factor/sfmin
    p2 = plot(title = "AE growth-rate bump gamma(scale factor)", xlabel = "scale factor / sfmin",
              ylabel = "gamma_AE  [c_s/a]", legend = :topleft, xscale = :log10, size = (760, 520))
    for r in rows
        b = get(bumps, r.ix, nothing); b === nothing && continue
        f_ad, g_ad, _, _, sf_ad, _ = b
        x = f_ad ./ sf_ad
        lab = @sprintf("IR%d rho=%.2f w=%.2f%s", r.ir, r.rho, r.ad_w, iscore(r) ? " (core)" : " (edge)")
        plot!(p2, x, max.(g_ad, 1e-8); label = lab, linewidth = 2,
              linestyle = iscore(r) ? :solid : :dash)
    end
    vline!(p2, [1.0]; label = "marginal (sfmin)", color = :black, linestyle = :dot)
    savefig(p2, joinpath(OUTDIR, "edge_vs_core_gamma_bump.png"))

    # ---- Panel 3: ballooning width <theta^2> vs rho (AD narrow vs wide grid)
    rr = sort(rows, by = r -> r.rho)
    p3 = plot(title = "ballooning width <theta^2> vs radius", xlabel = "rho",
              ylabel = "<theta^2>  (|phi|-weighted)", legend = :topleft, size = (760, 520))
    plot!(p3, [r.rho for r in rr], [r.ad_th2 for r in rr]; label = "AD-located (narrow)",
          marker = :circle, linewidth = 2, color = :firebrick)
    grm = [(r.rho, r.gr_th2) for r in rr if isfinite(r.gr_th2)]
    if !isempty(grm)
        plot!(p3, first.(grm), last.(grm); label = "grid w>=1 (wide)",
              marker = :utriangle, linewidth = 2, color = :seagreen, linestyle = :dash)
    end
    savefig(p3, joinpath(OUTDIR, "edge_vs_core_width.png"))

    println("wrote edge_vs_core_eigenfunction.png / _gamma_bump.png / _width.png")
end

main()
