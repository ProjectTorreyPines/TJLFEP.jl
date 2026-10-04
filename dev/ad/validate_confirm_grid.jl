# §7.1 A/B exactness for the cheap-rank→few-confirm `confirm_grid` path in `critical_factor_robust`.
#
# For a core radius and a plasma-edge radius on DIII-D n_scan=20, at small N_BASIS, assert that the
# new cheap-rank path returns the SAME grid minimum as the old brute-faithful path while paying fewer
# IFLUX=true (faithful) evals:
#
#   (A) adaptive=false, refine_rounds=1 (identical node set):
#         confirm_grid=false (old brute) vs confirm_grid=true (new) must match BITWISE on
#         sfmin, sfmin_w1, kyhat, width, binding; total_evals_full must be LOWER with confirm_grid=true.
#   (B) adaptive=true (default): sfmin must match to ≤1e-9 relative (the `sparse` zoom trigger now
#         keys off the cheap AE-unstable count, which can change which refine boxes are explored).
#   (C) critical_factor_truth with/without confirm_grid must match sfmin/sfmin_w1 to ≤1e-9 relative.
#
# CPU only (no GPU), fast. Env: RADII(22,95) NB_LIST(8,16).

using TJLF, TJLFEP, Printf

const USE_GPU = get(ENV, "USE_GPU", "0") == "1"
if USE_GPU
    using CUDA
    @assert CUDA.functional() "USE_GPU=1 but no functional GPU"
    CUDA.device!(first(CUDA.devices()))
end
const KW_INNER = :threads
const KW_TEAM  = nothing

# Prefixed names avoid colliding with constants (e.g. GACODE) baked into Main by the GPU sysimage.
const CG_CASE   = joinpath(@__DIR__, "..", "..", "examples", "DIIID_202017C42_500ms_v3.1")
const CG_GACODE = joinpath(CG_CASE, "input.gacode")
const CG_RADII  = parse.(Int, split(get(ENV, "RADII", "22,95"), ','))
const CG_NB_LIST = parse.(Int, split(get(ENV, "NB_LIST", "8,16"), ','))

reldiff(a, b) = (isfinite(a) && isfinite(b)) ? (a == b ? 0.0 : abs(a - b) / max(abs(a), abs(b), eps())) :
                (a === b || (isnan(a) && isnan(b)) ? 0.0 : Inf)

# track pass/fail across all asserts so the script exits nonzero on any failure
const FAILS = Ref(0)
function check(cond, msg)
    if cond
        @printf("    [ok]   %s\n", msg)
    else
        @printf("    [FAIL] %s\n", msg); FAILS[] += 1
    end
    flush(stdout)
end

function run_case(opts, prof, ir, nb)
    ep = deepcopy(opts); ep.IR = ir; ep.N_BASIS = nb
    gth = TJLFEP._gamma_thresh_for(ep, prof); shi = Float64(ep.FACTOR_IN); slo = shi / 512.0
    kw = (; gamma_thresh = gth, scan_lo = slo, scan_hi = shi, inner = KW_INNER, team = KW_TEAM, use_gpu = USE_GPU)

    @printf("\n===== IR=%d  NB=%d =====\n", ir, nb); flush(stdout)

    # (A) identical node set: adaptive=false, refine_rounds=1
    old = critical_factor_robust(ep, prof; confirm_grid = false, adaptive = false, refine_rounds = 1, kw...)
    new = critical_factor_robust(ep, prof; confirm_grid = true,  adaptive = false, refine_rounds = 1, kw...)
    @printf("  (A) adaptive=false refine=1\n")
    @printf("      old : sfmin=%.10g sfmin_w1=%.10g ky=%.6g w=%.6g bind=%s full=%d eig=%d\n",
            old.sfmin, old.sfmin_w1, old.kyhat, old.width, String(old.binding), old.total_evals_full, old.total_evals_eig)
    @printf("      new : sfmin=%.10g sfmin_w1=%.10g ky=%.6g w=%.6g bind=%s full=%d eig=%d\n",
            new.sfmin, new.sfmin_w1, new.kyhat, new.width, String(new.binding), new.total_evals_full, new.total_evals_eig)
    check(old.sfmin == new.sfmin,         @sprintf("sfmin bitwise equal (%.17g)", new.sfmin))
    check(old.sfmin_w1 == new.sfmin_w1,   @sprintf("sfmin_w1 bitwise equal (%.17g)", new.sfmin_w1))
    check(old.kyhat == new.kyhat,         @sprintf("kyhat bitwise equal (%.17g)", new.kyhat))
    check(old.width == new.width,         @sprintf("width bitwise equal (%.17g)", new.width))
    check(old.binding === new.binding,    @sprintf("binding equal (%s)", String(new.binding)))
    check(new.total_evals_full <= old.total_evals_full,
          @sprintf("confirm_grid full evals not higher: new=%d old=%d (save %d)",
                   new.total_evals_full, old.total_evals_full, old.total_evals_full - new.total_evals_full))

    # (B) default adaptive=true
    olda = critical_factor_robust(ep, prof; confirm_grid = false, adaptive = true, kw...)
    newa = critical_factor_robust(ep, prof; confirm_grid = true,  adaptive = true, kw...)
    @printf("  (B) adaptive=true (default)\n")
    @printf("      old : sfmin=%.10g sfmin_w1=%.10g full=%d\n", olda.sfmin, olda.sfmin_w1, olda.total_evals_full)
    @printf("      new : sfmin=%.10g sfmin_w1=%.10g full=%d\n", newa.sfmin, newa.sfmin_w1, newa.total_evals_full)
    check(reldiff(olda.sfmin, newa.sfmin) <= 1e-9,
          @sprintf("sfmin rel diff %.2e <= 1e-9", reldiff(olda.sfmin, newa.sfmin)))
    check(reldiff(olda.sfmin_w1, newa.sfmin_w1) <= 1e-9,
          @sprintf("sfmin_w1 rel diff %.2e <= 1e-9", reldiff(olda.sfmin_w1, newa.sfmin_w1)))

    # (C) truth parity
    nb_steps = [nb, nb + 8]
    to = critical_factor_truth(ep, prof; confirm_grid = false, nb_work = nb, nb_steps = nb_steps, kw...)
    tn = critical_factor_truth(ep, prof; confirm_grid = true,  nb_work = nb, nb_steps = nb_steps, kw...)
    @printf("  (C) truth (nb_work=%d nb_steps=%s)\n", nb, string(nb_steps))
    @printf("      old : sfmin=%.10g sfmin_w1=%.10g\n", to.sfmin, to.sfmin_w1)
    @printf("      new : sfmin=%.10g sfmin_w1=%.10g\n", tn.sfmin, tn.sfmin_w1)
    check(reldiff(to.sfmin, tn.sfmin) <= 1e-9,
          @sprintf("truth sfmin rel diff %.2e <= 1e-9", reldiff(to.sfmin, tn.sfmin)))
    check(reldiff(to.sfmin_w1, tn.sfmin_w1) <= 1e-9,
          @sprintf("truth sfmin_w1 rel diff %.2e <= 1e-9", reldiff(to.sfmin_w1, tn.sfmin_w1)))

    return (; ir, nb, save = old.total_evals_full - new.total_evals_full,
            old_full = old.total_evals_full, new_full = new.total_evals_full)
end

function main()
    @printf("DIII-D confirm_grid A/B exactness (%s)  radii=%s  NB=%s\n",
            USE_GPU ? "GPU" : "CPU", string(CG_RADII), string(CG_NB_LIST))
    flush(stdout)
    rows = NamedTuple[]
    warmed = false
    for nb in CG_NB_LIST
        # build options at the FIRST nb's TGLFEP, reuse for all radii at this nb
        tglfep = joinpath(CG_CASE, "input_scan20_nb$(nb).TGLFEP")
        opts, prof, _ = preprocess_gacode_inputs(CG_GACODE, tglfep)
        opts.N_BASIS = nb
        if !warmed
            let ep = deepcopy(opts); ep.IR = CG_RADII[1]
                critical_factor_robust(ep, prof; confirm_grid = true, adaptive = false, refine_rounds = 0,
                    inner = KW_INNER, team = KW_TEAM, use_gpu = USE_GPU)
            end
            warmed = true
        end
        for ir in CG_RADII
            push!(rows, run_case(opts, prof, ir, nb))
        end
    end

    println("\n===================== SUMMARY =====================")
    @printf("  %-5s %-4s %-10s %-10s %-8s\n", "IR", "NB", "old_full", "new_full", "saved")
    for r in rows
        @printf("  %-5d %-4d %-10d %-10d %-8d\n", r.ir, r.nb, r.old_full, r.new_full, r.save)
    end
    if FAILS[] == 0
        println("\n=== ALL CHECKS PASSED ===")
    else
        @printf("\n=== %d CHECK(S) FAILED ===\n", FAILS[])
        exit(1)
    end
end

main()
