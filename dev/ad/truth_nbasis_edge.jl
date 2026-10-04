# Tier-3 of the :truth protocol in isolation: PIN the located (ky,width) optimum
# and walk the nbasis ladder with marginal_factor_faithful (NOT a per-nb
# leading-mode reselection, which is what made IR95 mode-hop). This is the exact
# +nbasis correction :truth applies on top of robust_ad. Writes each (radius,nb)
# result to dev/ad/truth_nbasis_edge.txt immediately so progress is visible.
#
#   module load julia/1.11.7 && export JULIA_DEPOT_PATH=$PSCRATCH/.julia
#   julia -t 8 --project=. dev/ad/truth_nbasis_edge.jl

using TJLFEP, Printf, Serialization
const CASE   = joinpath(@__DIR__, "..", "..", "examples", "DIIID_202017C42_500ms_v3.1")
const GACODE = joinpath(CASE, "input.gacode")
const TGLFEP = joinpath(CASE, "input_scan20_nb32.TGLFEP")
const AD_TASKS = joinpath(@__DIR__, "ad_threads_sfmin_nb32_ad_tasks")
const OUT = joinpath(@__DIR__, "truth_nbasis_edge.txt")

opts0, prof, _ = preprocess_gacode_inputs(GACODE, TGLFEP)
gth = TJLFEP._gamma_thresh_for(opts0, prof)
shi = Float64(opts0.FACTOR_IN); slo = shi/512.0
nb_steps = parse.(Int, split(get(ENV,"NB_STEPS","32,40,48,56"), ","))
scan_ir  = [(18,90,0.89),(19,95,0.94),(20,101,1.00)]

logline(s) = (open(OUT,"a") do io; println(io, s); end; println(s); flush(stdout))
rm(OUT, force=true)
logline(@sprintf("pinned nbasis ladder  gamma_thresh=%.3g  nb_steps=%s", gth, string(nb_steps)))

for (ix,ir,rho) in scan_ir
    at = deserialize(joinpath(AD_TASKS, "task_$(ix).jls"))
    ky = at.kymark; w = at.width
    logline(@sprintf("\n=== ix=%d IR=%d rho=%.2f  ky=%.3f w=%.3f (archived sfmin=%.4g) ===", ix,ir,rho,ky,w,at.sfmin))
    nbs = Int[]; vals = Float64[]
    for nb in nb_steps
        epn = deepcopy(opts0); epn.IR = ir; epn.N_BASIS = nb
        t0 = time()
        v = try
            r = TJLFEP.marginal_factor_faithful(epn, prof; kyhat=ky, width=w, gamma_thresh=gth,
                    scan_lo=slo, scan_hi=shi, threaded=true, inner=:threads, use_gpu=false)
            (r.binding !== :none && isfinite(r.factor_faithful)) ? r.factor_faithful : NaN
        catch e
            NaN
        end
        push!(nbs, nb); push!(vals, v)
        logline(@sprintf("  nb=%-3d  sfmin=%-10.5g  (%.1fs)", nb, v, time()-t0))
    end
    ex = TJLFEP._nbasis_extrapolate(nbs, vals)
    logline(@sprintf("  -> ratio=%.3g  converged=%s  limit=%.5g", ex.ratio, string(ex.converged), ex.limit))
end
