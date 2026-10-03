# Formalized `:wide` width-extension (critical_factor_optimize extend_mode=:wide) vs the original
# w≥1-only :ad and robust_ad, over the 20-radius DIII-D scan. Sweeps wide_kdesc = 1,2,3 to show how
# the small multistart closes the residual over-prediction gap to robust_ad and at what cost.
#
#   module load julia/1.11.7 && export JULIA_DEPOT_PATH=$PSCRATCH/.julia
#   USE_GPU=1 NB=8 julia --project=. dev/ad/profile_ad_wide_kdesc.jl
# env: NB (default 8), KDESC (default 1,2,3), RADII (default all 20), WITH_ROBUST (default 1)

using TJLF, TJLFEP, Printf

const USE_GPU = get(ENV, "USE_GPU", "0") == "1"
if USE_GPU
    using CUDA
    @assert CUDA.functional() "USE_GPU=1 but no functional GPU"
    CUDA.device!(first(CUDA.devices()))
end

const AX_CASE   = joinpath(@__DIR__, "..", "..", "examples", "DIIID_202017C42_500ms_v3.1")
const AX_GACODE = joinpath(AX_CASE, "input.gacode")
const AX_NB     = parse(Int, get(ENV, "NB", "8"))
const KDESC     = parse.(Int, split(get(ENV, "KDESC", "1,2,3"), ','))
const WITH_ROB  = get(ENV, "WITH_ROBUST", "1") == "1"

faithful_of(r) = (r.faithful !== nothing && r.faithful.binding != :none) ? r.faithful.factor_faithful : r.sfmin

function run_case(opts, prof, ir)
    ep = deepcopy(opts); ep.IR = ir
    gth = TJLFEP._gamma_thresh_for(ep, prof); shi = Float64(ep.FACTOR_IN); slo = shi / 512.0
    kw = (; gamma_thresh = gth, scan_lo = slo, scan_hi = shi, inner = :threads, team = nothing, use_gpu = USE_GPU)
    @printf("\n===== IR=%d  NB=%d  (FACTOR_IN=%.4g) =====\n", ir, ep.N_BASIS, shi); flush(stdout)

    row = Dict{Symbol,Any}(:ir => ir)

    t = @elapsed o = critical_factor_optimize(ep, prof; faithful_confirm = true, extend_width = false, kw...)
    row[:orig] = faithful_of(o); row[:orig_t] = t; row[:orig_ev] = o.evals
    @printf("  orig            sfmin=%.5g  (ky=%.3f w=%.3f)  ev=%d  %.1fs\n", row[:orig], o.kyhat, o.width, o.evals, t); flush(stdout)

    for k in KDESC
        t = @elapsed w = critical_factor_optimize(ep, prof; faithful_confirm = true, extend_width = true,
                extend_mode = :wide, wide_kdesc = k, kw...)
        row[Symbol("w$k")] = faithful_of(w); row[Symbol("w$(k)_t")] = t; row[Symbol("w$(k)_ev")] = w.evals
        @printf("  wide k=%d        sfmin=%.5g  (ky=%.3f w=%.3f)  ev=%d  nx=%d  %.1fs\n",
                k, row[Symbol("w$k")], w.kyhat, w.width, w.evals, w.n_ext_confirm, t); flush(stdout)
    end

    if WITH_ROB
        t = @elapsed rob = critical_factor_robust(ep, prof; extend_width = true, kw...)
        row[:robust] = rob.sfmin; row[:robust_t] = t; row[:robust_ev] = rob.total_evals_full + rob.total_evals_eig
        @printf("  robust_ad (ref) sfmin=%.5g  (ky=%.3f w=%.3f)  ev=%d  %.1fs\n",
                rob.sfmin, rob.kyhat, rob.width, row[:robust_ev], t); flush(stdout)
    end
    return row
end

function main()
    @printf("formalized :wide kdesc sweep (%s)  NB=%d  KDESC=%s  threads=%d\n",
            USE_GPU ? "GPU" : "CPU", AX_NB, string(KDESC), Threads.nthreads()); flush(stdout)
    tglfep = joinpath(AX_CASE, "input_scan20_nb$(AX_NB).TGLFEP")
    opts, prof, _ = preprocess_gacode_inputs(AX_GACODE, tglfep)
    opts.N_BASIS = AX_NB

    ir_exp = (!isempty(opts.IR_EXP) && !all(iszero, opts.IR_EXP)) ? Int.(opts.IR_EXP) :
             TJLFEP.ir_exp_from_scan(prof.NR, prof.IRS, opts.SCAN_N)
    radii = haskey(ENV, "RADII") ? parse.(Int, split(ENV["RADII"], ',')) : Int.(ir_exp[1:opts.SCAN_N])
    @printf("radii (%d): %s\n", length(radii), string(radii)); flush(stdout)

    run_case(opts, prof, first(radii))   # warm up compilation (excluded from timings below)
    rows = [run_case(opts, prof, ir) for ir in radii]

    kcols = ["w$k" for k in KDESC]
    println("\n===================== SUMMARY (sfmin: ratio to robust) =====================")
    @printf("  %-5s %-11s %s %-11s\n", "IR", "orig", join([@sprintf("%-11s", "wide k=$k") for k in KDESC]), "robust")
    for r in rows
        ro = get(r, :robust, NaN)
        kvals = join([@sprintf("%-11.5g", get(r, Symbol(c), NaN)) for c in kcols])
        @printf("  %-5d %-11.5g %s %-11.5g\n", r[:ir], get(r, :orig, NaN), kvals, ro)
    end
    # mean |ratio-1| of each wide-k vs robust (over radii where both finite)
    if WITH_ROB
        println("\n  accuracy vs robust_ad (mean |wide/robust - 1|, max ratio):")
        for k in KDESC
            rels = Float64[]; ratios = Float64[]
            for r in rows
                v = get(r, Symbol("w$k"), NaN); ro = get(r, :robust, NaN)
                (isfinite(v) && isfinite(ro) && ro > 0) || continue
                push!(rels, abs(v/ro - 1)); push!(ratios, v/ro)
            end
            @printf("    wide k=%d : mean|rel|=%.1f%%   max ratio=%.2fx\n", k, 100*sum(rels)/length(rels), maximum(ratios))
        end
    end

    println("\n===================== SUMMARY (cost: sum wallclock s / mean evals) =====================")
    to = sum(get(r, :orig_t, 0.0) for r in rows); tr = sum(get(r, :robust_t, 0.0) for r in rows)
    @printf("  orig:   sum=%.1fs  mean_ev=%d\n", to, round(Int, sum(get(r,:orig_ev,0) for r in rows)/length(rows)))
    for k in KDESC
        tk = sum(get(r, Symbol("w$(k)_t"), 0.0) for r in rows)
        ev = round(Int, sum(get(r, Symbol("w$(k)_ev"), 0) for r in rows)/length(rows))
        @printf("  wide k=%d: sum=%.1fs  mean_ev=%d   (%.2fx orig)\n", k, tk, ev, tk/max(to,eps()))
    end
    WITH_ROB && @printf("  robust: sum=%.1fs  mean_ev=%d\n", tr, round(Int, sum(get(r,:robust_ev,0) for r in rows)/length(rows)))
end

main()
