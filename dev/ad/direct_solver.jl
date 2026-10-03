# Shared prototype solvers for the (ky,w) critical-factor search, used by the offline experiments.
# Both reuse TJLFEP internals (_ae_unstable_window, marginal_factor_faithful, critical_factor_optimize)
# and the team/threads fan-out (_ad_pmap). Kept here (not in core) while we benchmark; promote once a
# winner is chosen.
#
#   critical_factor_direct      — NLopt GN_DIRECT_L global search on the cheap AE-onset surface +
#                                 early-stop few-confirm (the accuracy ceiling; cost set by max_evals).
#   critical_factor_ad_f1seed   — pinned-aware f1 seed grid (confirm's front-end, so floor-pinned
#                                 basins are SEEN) → :ad gradient descent on INTERIOR basins (slides
#                                 off-node) → early-stop few-confirm. Aims for :ad-class speed with
#                                 confirm/DIRECT robustness.

using NLopt

# ── Option 2 reference: DIRECT global search on cheap AE-edge + early-stop confirm ──
function critical_factor_direct(ep0, prof; gamma_thresh=nothing,
                                scan_lo=nothing, scan_hi=nothing,
                                n_eig::Int=24, max_evals::Int=40, ky_lo::Float64=0.25,
                                inner::Symbol=:threads, team=nothing,
                                use_gpu::Bool=false, verbose::Bool=false)
    gth = gamma_thresh === nothing ? TJLFEP._gamma_thresh_for(ep0, prof) : gamma_thresh
    shi = scan_hi === nothing ? Float64(ep0.FACTOR_IN) : scan_hi
    slo = scan_lo === nothing ? shi / 512.0 : scan_lo
    kylo, kyhi = ky_lo, 1.0   # ky_lo=0.25 matches the canonical kwscale_scan kyhat floor {0.25,..,1.0}
    wlo, whi = Float64(ep0.WIDTH_MIN), Float64(ep0.WIDTH_MAX)
    PENALTY = shi

    samples = NamedTuple[]
    eig_evals = Ref(0)
    cheap_onset = function (ky, w)
        ep = deepcopy(ep0); ep.KYHAT_IN = ky; ep.WIDTH_IN = w
        win = TJLFEP._ae_unstable_window(ep, prof, gth; scan_lo=slo, scan_hi=shi,
                  n_eig=n_eig, threaded=true, use_gpu=use_gpu, inner=inner, team=team)
        eig_evals[] += win.evals
        f = win.unstable ? win.f1 : PENALTY
        push!(samples, (; ky=ky, w=w, f=f, unstable=win.unstable, pinned=win.pinned_lo))
        return f
    end

    opt = NLopt.Opt(:GN_DIRECT_L, 2)
    NLopt.lower_bounds!(opt, [kylo, wlo])
    NLopt.upper_bounds!(opt, [kyhi, whi])
    NLopt.maxeval!(opt, max_evals)
    NLopt.min_objective!(opt, (x, grad) -> cheap_onset(x[1], x[2]))
    (_minf, _minx, ret) = NLopt.optimize(opt, [0.5*(kylo+kyhi), 0.5*(wlo+whi)])

    order = sortperm([Float64(s.f) for s in samples])
    best_f = Inf; best_ky=NaN; best_w=NaN; best_bind=:none; total_full=0; n_confirm=0
    for i in order
        s = samples[i]
        s.unstable || break
        s.f >= best_f && break
        r = TJLFEP.marginal_factor_faithful(ep0, prof; kyhat=s.ky, width=s.w,
                gamma_thresh=gth, scan_lo=slo, scan_hi=shi, threaded=true,
                inner=inner, team=team, use_gpu=use_gpu)
        n_confirm += 1; total_full += r.evals_full; eig_evals[] += r.evals_eig
        if r.binding != :none && isfinite(r.factor_faithful) && r.factor_faithful < best_f
            best_f = r.factor_faithful; best_ky=s.ky; best_w=s.w; best_bind=r.binding
        end
    end
    status = (best_bind===:none || !isfinite(best_f)) ? :no_onset : (best_f >= 0.999*shi ? :cap : :ok)
    return (; sfmin=best_f, kyhat=best_ky, width=best_w, binding=best_bind, status=status,
            n_samples=length(samples), n_confirm=n_confirm, nlopt_ret=ret,
            total_evals_full=total_full, total_evals_eig=eig_evals[])
end

# NOTE: critical_factor_ad_f1seed was promoted to core (src/tjlfep_ad_extensions.jl) and is exported
# by TJLFEP — the escalate harness below calls the core version via `using TJLFEP`.

# ── Production policy: fast adf1 default, escalate to a robust global search only on flagged radii ──
# adf1 wins clean/pinned radii cheaply, but on multimodal/keep-filter-divergent surfaces it can land
# above the true min (IR=48 cheap/faithful gap; IR=95 sparse). Cheap trust gate decides per radius:
#   cheap_gap > gap_thresh     → winner's faithful ≫ its cheap f1: descent chased a keep-filter basin.
#   feasible_frac < feas_thresh→ sparse unstable seed grid: surface multimodal/under-bracketed.
#   adf1 no_onset / on cap     → nothing trustworthy found.
# Flagged → run the escalation target (:direct = canonical-bounds DIRECT-40, or :grid = grid-zoom
# robust_ad), then keep the LOWER faithful sfmin (both are faithful-confirmed lower-is-better globals).
function critical_factor_ad_escalate(ep0, prof; escalate_to::Symbol=:direct,
                                     gap_thresh::Float64=1.5, feas_thresh::Float64=0.25,
                                     direct_evals::Int=40, direct_neig::Int=24,
                                     ky_lo::Float64=0.25, nseed_ky::Int=4, nseed_w::Int=8, n_eig_seed::Int=12,
                                     gamma_thresh=nothing, scan_lo=nothing, scan_hi=nothing,
                                     inner::Symbol=:threads, team=nothing,
                                     use_gpu::Bool=false, verbose::Bool=false)
    @assert escalate_to in (:direct, :grid) "escalate_to must be :direct or :grid"
    r = critical_factor_ad_f1seed(ep0, prof; ky_lo=ky_lo, nseed_ky=nseed_ky, nseed_w=nseed_w,
            n_eig_seed=n_eig_seed, gamma_thresh=gamma_thresh, scan_lo=scan_lo, scan_hi=scan_hi,
            inner=inner, team=team, use_gpu=use_gpu)
    reasons = Symbol[]
    (r.status === :no_onset)            && push!(reasons, :no_onset)
    (r.status === :cap)                 && push!(reasons, :cap)
    (r.cheap_gap > gap_thresh)          && push!(reasons, :cheap_gap)
    (r.feasible_frac < feas_thresh)     && push!(reasons, :sparse)
    flagged = !isempty(reasons)
    full = r.total_evals_full; eig = r.total_evals_eig
    if !flagged
        return (; sfmin=r.sfmin, kyhat=r.kyhat, width=r.width, binding=r.binding, status=r.status,
                escalated=false, reasons=reasons, source=:adf1,
                adf1_sfmin=r.sfmin, esc_sfmin=NaN, cheap_gap=r.cheap_gap, feasible_frac=r.feasible_frac,
                total_evals_full=full, total_evals_eig=eig)
    end
    if escalate_to === :direct
        e = critical_factor_direct(ep0, prof; max_evals=direct_evals, n_eig=direct_neig, ky_lo=ky_lo,
                gamma_thresh=gamma_thresh, scan_lo=scan_lo, scan_hi=scan_hi,
                inner=inner, team=team, use_gpu=use_gpu)
    else
        e = TJLFEP.critical_factor_robust(ep0, prof; gamma_thresh=gamma_thresh,
                scan_lo=scan_lo, scan_hi=scan_hi, inner=inner, team=team, use_gpu=use_gpu)
    end
    full += get(e, :total_evals_full, 0); eig += get(e, :total_evals_eig, 0)
    # Keep the better (lower) faithful onset; both are faithful-confirmed globals.
    use_esc = isfinite(e.sfmin) && (!isfinite(r.sfmin) || e.sfmin < r.sfmin)
    win = use_esc ? e : r
    verbose && @info "escalate" reasons=reasons adf1=r.sfmin esc=e.sfmin chosen=(use_esc ? escalate_to : :adf1)
    return (; sfmin=win.sfmin, kyhat=getfield(win, :kyhat), width=getfield(win, :width),
            binding=win.binding, status=win.status,
            escalated=true, reasons=reasons, source=(use_esc ? escalate_to : :adf1),
            adf1_sfmin=r.sfmin, esc_sfmin=e.sfmin, cheap_gap=r.cheap_gap, feasible_frac=r.feasible_frac,
            total_evals_full=full, total_evals_eig=eig)
end
