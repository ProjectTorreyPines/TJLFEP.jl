> **Historical paths.** This note predates the repo restructure: `build/ad/` is now `dev/ad/`, `build/timing/` is `dev/timing/`, `build/common/` is `slurm/common/`, and production submit scripts live in `slurm/`.

# Plan: make `:robust_ad` faster at the same accuracy (cheap-rank → few-confirm on the `w≥1` core)

Status: IMPLEMENTED, but the speed premise did **not** hold on the GPU/MPS production route, so
`confirm_grid` ships **off by default** (opt-in). See the **OUTCOME** section directly below.
Repo: `TJLFEP` (Julia 1.11.7). Solver code: `src/tjlfep_ad_extensions.jl`.

---

## OUTCOME (executed 2026-06-18)

What was built (per §4): a shared `_rank_confirm(ep, prof, pts; …)` helper (cheap AE-band rank →
early-stop few-confirm); wired into `critical_factor_robust.eval_grid!` behind a `confirm_grid` kwarg;
the `sparse` zoom trigger fed the cheap AE-unstable count in that mode; `confirm_grid` forwarded
through `critical_factor_truth`; and `critical_factor_confirm` refactored to reuse `_rank_confirm`
(dedup). Confirm faithful calls are wrapped in `with_blas_threads(1)` to stay bitwise-identical to the
brute `_ad_pmap` path.

Results:
- **Exactness (§7.1): CONFIRMED.** With `adaptive=false` (identical node set), `confirm_grid` true vs
  false is **bitwise identical** (sfmin/sfmin_w1/ky/w/binding) at DIII-D IR=2/17/22/95 × `nb=6/8/16`,
  and `total_evals_full` drops sharply (e.g. 780→29, 824→28, 703→162). The §5 proof holds.
- **Speed (§7.2 / primary goal): NOT ACHIEVED — it got *slower*.** Same 1-node GPU/MPS harness,
  20-radius scan: brute→confirm scan time went **621→735s (nb6), 710→780 (nb8), 1499→1602 (nb16),
  4523→4966 (nb32)** — i.e. **+7…18%**. Root cause: the early-stop confirm is **serial** (each confirm
  depends on the running incumbent) so it runs on **one** GPU while the brute path fans all 32 nodes
  across the 4-GPU MPS team; plus the cheap rank adds a whole extra eigen-scan pass. The plan's premise
  — that the `IFLUX=true` keep filter dominates per-node cost — is false on this hardware (the GPU
  **eigensolve** dominates, and brute already parallelizes it). At `nb=32`, confirm reproduces the brute
  `sfmin` **exactly** (0 rel diff, all radii) yet is ~10% slower: a pure loss.
- **Adaptive parity (§4.3 risk): MATERIALIZED.** With `adaptive=true`, the cheap (over-counted)
  feasibility count shifts the `sparse` trigger, so `sfmin` differs from the brute baseline at
  near-degenerate radii (IR=2/17, ~3–10%, only `nb<32`). Reproduced **in a single job** via the A/B
  (B)-case (brute=0.16402 vs confirm=0.15875 at IR=17/nb6), ruling out GPU non-determinism. (Confirm
  there actually matched the *more-refined* `adaptive=false` value; brute's adaptive stopped early — so
  it is not a correctness regression, but it does break strict baseline parity.)

Decision: keep the implementation (correct, useful on CPU/serial and as an exact opt-in with
`adaptive=false`), but **default `confirm_grid=false`** in `critical_factor_robust` and
`critical_factor_truth` so the production GPU `:robust_ad`/`:truth` route is byte-for-byte unchanged
(exact baseline parity, no slowdown). No `mainsub` plumbing was added. The node-hours re-collect/plot
step was skipped because the measured direction is a regression, not a win. Validation harnesses:
`build/ad/validate_confirm_grid.jl` (+ `batch_validate_confirm_grid{,_gpu}.sh`),
`build/timing/check_sfmin_parity.jl`.

A future win would require **batched-parallel** confirm (confirm the top-k cheap candidates across the
MPS team per round, then early-stop) to avoid the serial single-GPU tail — out of scope here.

---

## 1. Goal

`:robust_ad` is currently ~4× the node-hours of `:grid` MPS at `N_BASIS=32`
(1.26 vs 0.31 node-hours on the 20-radius DIII-D scan) while being more accurate.
Most of that extra cost is **not** the autodiff narrow-width extension — it is that
the `w≥1` canonical-box stage brute-forces the **faithful** (`IFLUX=true`, full
keep-filter) onset at *every* grid node.

We already proved (and built) the fix as the standalone `:confirm` solver
(`critical_factor_confirm`): rank the grid on the **cheap eigenvalue-only AE-band
onset** (`IFLUX=false`) and faithful-confirm only the few nodes that can still hold
the minimum. This is **provably exact** (faithful onset ≥ cheap AE-band onset
node-wise, because keep ⊆ AE-unstable), so it returns the same grid minimum while
paying `IFLUX=true` on ~1–few nodes instead of all 32.

**This plan folds that cheap-rank→few-confirm step into `critical_factor_robust`'s
`w≥1` passes** (`eval_grid!` + refine), so `:robust_ad` keeps its narrow-width
extension (`_locate_extended`) but stops over-paying on the core.

Expected outcome: `:robust_ad` `sfmin` per radius **unchanged**; node-hours drop
toward the `:grid` MPS line, especially at large `N_BASIS`.

---

## 2. Why this is the right change (and what was already tried)

- `critical_factor_confirm` (`:confirm`) **exists, is exported, validated, exact**,
  but (a) it is **not** wired into `mainsub`'s `solver` dispatch (only
  `:grid,:ad,:robust_ad,:truth`), so no scan/timing/DB-gen run ever used it; and
  (b) it covers **only `w≥1`** (it equals `critical_factor_robust(refine=0)` over
  `w∈[WIDTH_MIN,WIDTH_MAX]`, with **no** narrow-width extension). So `:confirm` is a
  faster `:grid`, not a faster `:robust_ad`.
- `_locate_extended` (the `w<1` narrow path inside `robust_ad`) **already** uses the
  cheap-rank→confirm pattern. So the asymmetry is purely on the `w≥1` core, which
  still confirms all nodes.
- Net: the algorithm is built and trusted; the only missing work is **reusing it
  inside `critical_factor_robust`'s `w≥1` grid passes**. No new numerics invented.

---

## 3. Code map (current state, `src/tjlfep_ad_extensions.jl`)

- `critical_factor_robust(...)` ≈ lines **1082–1282**.
  - `eval_grid!(kya,kyb,wa,wb)` ≈ **1119–1147**: builds a `nkyhat×nefwid` (4×8)
    **linear** `(ky,w)` grid and calls `marginal_factor_faithful` at **every** node
    via `_ad_pmap`; folds the min into the closure vars `best_f/best_ky/best_w/best_bind`;
    returns `(; nfeasible, np, kymin, kymax, wmin, wmax)` used by the zoom trigger.
  - `zoom_decision(bky,bw,st)` ≈ **1160–1166**: `boundary || sparse || nearcap`.
    `sparse` uses `st.nfeasible` (count of faithful binding≠none nodes).
  - coarse pass **1169**, hybrid polish **1184–1203** (`outer=:hybrid` AD-IFT polish),
    refine loop **1215–1233**, width extension **1251–1268** (`_locate_extended`),
    snapshot `sfmin_w1` at **1249**, return **1276–1281**.
- `critical_factor_confirm(...)` ≈ **1313–1377**: the cheap-rank→confirm-few logic to extract.
  - cheap AE-band onset per node via `_ae_unstable_window(...; n_eig=24)` **1336–1342**;
  - `order = sortperm([c.f ...])`, confirm in increasing cheap order, early stop
    `c.f >= best_f && break` **1351–1369**.
- `_locate_extended(...)` ≈ **1654–1690**: reference for the same pattern already in use
  (cheap rank `1663–1668`, sorted candidate confirm with `c[3] >= best_f && break`).
- `mainsub` solver dispatch: `src/mainsub.jl` **37–39** (`:grid,:ad,:robust_ad,:truth`).
- Export list: `src/TJLFEP.jl` line ~**53**.

---

## 4. The change

### 4.1 Extract a shared helper

Add (near `critical_factor_confirm`) a helper that ranks a fixed point list on the
cheap AE-band onset and faithful-confirms in increasing cheap order with early-stop
against a caller-supplied incumbent. Reuse it from BOTH `critical_factor_confirm`
and `critical_factor_robust`'s `eval_grid!`.

Proposed signature (keyword-compatible with existing call sites):

```julia
# Cheap-rank → few-confirm over a fixed (ky,w) point list. Returns the updated
# incumbent and stats. EXACT: faithful onset ≥ cheap AE-band onset node-wise.
function _rank_confirm(ep0, prof, pts::Vector{Tuple{Float64,Float64}};
                       gth, slo, shi, incumbent_f::Float64,
                       n_eig::Int=24, inner::Symbol=:threads, team=nothing,
                       use_gpu::Bool=false)
    # 1) cheap AE-band lower edge f1 at every node (IFLUX=false), parallel over team/threads
    cheap = _ad_pmap(idx -> begin
            ky, w = pts[idx]; ep = deepcopy(ep0); ep.KYHAT_IN = ky; ep.WIDTH_IN = w
            win = _ae_unstable_window(ep, prof, gth; scan_lo=slo, scan_hi=shi,
                      n_eig=n_eig, threaded=false, use_gpu=use_gpu)
            (; ky, w, f=(win.unstable ? win.f1 : Inf), evals=win.evals)
        end, length(pts); inner=inner, team=team)
    total_eig = sum(c.evals for c in cheap; init=0)
    # 2) confirm in increasing cheap order, early stop when cheap ≥ best faithful
    best_f = incumbent_f; best_ky = NaN; best_w = NaN; best_bind = :none
    total_full = 0; n_confirm = 0; n_cheap_feasible = 0
    for i in sortperm([c.f for c in cheap])
        c = cheap[i]
        isfinite(c.f) || break
        n_cheap_feasible += 1
        c.f >= best_f && break
        r = marginal_factor_faithful(ep0, prof; kyhat=c.ky, width=c.w, gamma_thresh=gth,
                scan_lo=slo, scan_hi=shi, threaded=false, use_gpu=use_gpu)
        n_confirm += 1; total_full += r.evals_full; total_eig += r.evals_eig
        if r.binding != :none && isfinite(r.factor_faithful) && r.factor_faithful < best_f
            best_f = r.factor_faithful; best_ky = c.ky; best_w = c.w; best_bind = r.binding
        end
    end
    return (; best_f, best_ky, best_w, best_bind,
            n_confirm, n_cheap_feasible, total_evals_full=total_full, total_evals_eig=total_eig)
end
```

Notes:
- `n_eig`: `critical_factor_confirm` uses `n_eig=24`; `_locate_extended` uses
  `n_eig_seed=12`. Use **24** for the `w≥1` core (matches the existing exact `:confirm`).
- `_ad_pmap`, `_ae_unstable_window`, `marginal_factor_faithful` already exist and are
  used exactly this way in both reference functions.

### 4.2 Wire it into `critical_factor_robust.eval_grid!`

Replace the "confirm every node" body so it (1) builds the same `pts`, then (2) calls
`_rank_confirm` with the **current global incumbent `best_f`** as `incumbent_f`,
folding any win back into the closure's `best_f/best_ky/best_w/best_bind`. Accumulate
`total_full/total_eig`. Gate behind a kwarg:

```julia
confirm_grid::Bool = true   # add to critical_factor_robust signature
```

When `confirm_grid=false`, keep the old brute-faithful path (for A/B exactness checks
and the filter-heavy ITER-rotational regime where many nodes bind anyway).

### 4.3 Keep the zoom trigger working

`zoom_decision`'s `boundary` and `nearcap` use `best_f/best_ky/best_w` — still **exact**
(from confirmed faithful). Only `sparse` uses a feasibility count. With early-stop we
no longer confirm every node, so:

- Redefine the trigger's feasibility count as `n_cheap_feasible` (cheap AE-unstable
  count) returned by `_rank_confirm`. This is a **heuristic** (it over-counts vs
  faithful feasibility), so document it. It changes *which* refine boxes are explored,
  not the minimum over the nodes actually evaluated.
- **Risk:** cheap over-count makes `sparse` fire *less*, so a refine the faithful-based
  logic would have done might be skipped → could miss a lower off-node min. Mitigate by
  validating with `adaptive=false` first (fixed `refine_rounds`, identical node set),
  then re-enabling `adaptive` and confirming `sfmin` still matches (it should: the
  winning basin is bracketed by the coarse cheap rank either way).

### 4.4 Plumb `confirm_grid` through `:truth` and `mainsub` (optional but recommended)

- `critical_factor_truth` calls `critical_factor_robust` (line ~**1732**) — forward
  `confirm_grid` so `:truth` benefits too (its tier-2 locate is the same object).
- No `mainsub` signature change needed; `confirm_grid=true` is the new default inside
  `critical_factor_robust`. Keep `solver=:robust_ad` unchanged on the scan path.

---

## 5. Correctness argument (why `sfmin` is unchanged)

For any fixed `(ky,w)` node, the faithful keep onset `f★` satisfies `f★ ≥ f1` (the
cheap AE-band lower edge), because the kept set is a subset of the AE-unstable hull.
Confirming nodes in increasing `f1` order and stopping when `f1 ≥ best_f★` cannot skip
any node whose `f★ < best_f★`. Therefore the minimum over a grid pass equals the
minimum you'd get by confirming **all** nodes. The only behavioral change across the
adaptive zoom is the `sparse` trigger heuristic (§4.3), validated empirically.

---

## 6. Files to touch

- `src/tjlfep_ad_extensions.jl`
  - add `_rank_confirm` helper;
  - refactor `eval_grid!` in `critical_factor_robust` to use it; add `confirm_grid` kwarg;
  - swap `sparse` count to `n_cheap_feasible`;
  - (optional) refactor `critical_factor_confirm` to call `_rank_confirm` (dedupe);
  - forward `confirm_grid` from `critical_factor_truth`.
- `docs/AD_SOLVERS_AND_SEARCH_BOUNDS.md` — note `robust_ad`'s `w≥1` stage now uses the
  exact cheap-rank→confirm-few scheme (same minimum, fewer `IFLUX=true` evals).
- `build/README.md` — one-line update to the `SOLVER=robust_ad` description if needed.

---

## 7. Validation (correctness first, then speed)

All Julia runs need the environment:

```bash
module load julia/1.11.7 && export JULIA_DEPOT_PATH=$PSCRATCH/.julia
cd /pscratch/sd/t/tneiser/.julia/dev/TJLFEP
```

### 7.1 Unit A/B exactness (single radius, CPU, fast)

Pick 2–3 representative radii (a core radius and a plasma-edge radius, e.g. IR=2 and
IR=95 from the DIII-D `ir_exp` list). For each, at `N_BASIS=8` and `16`:

- Run `critical_factor_robust(... confirm_grid=false, adaptive=false, refine_rounds=1)`
  (current behavior) and `confirm_grid=true` (new). Assert `sfmin`, `kyhat`, `width`,
  `binding` are **identical** (bitwise / `==`), and `total_evals_full` is **lower** with
  `confirm_grid=true`.
- Then with `adaptive=true` (default): assert `sfmin` matches to ≤1e-9 relative.

Write this as a small script under `build/ad/` (e.g. `validate_confirm_grid.jl`) or a
`@testset`. Confirm `:truth` with/without `confirm_grid` matches too.

### 7.2 Full 20-radius scan parity + node-hours (GPU)

Use the existing single-node backfill timing harness (now `BACKFILL_MODE`):

```bash
cd build
# robust_ad node-hours sweep (regular QoS by default; premium for speed)
SOLVER=robust_ad NB_LIST="6 8 16 32" QOS=premium ./timing/submit_nodehours_vs_nbasis.sh
# after they finish:
julia --project=. timing/collect_scan20_timing.jl   # -> timing_runs/scan20_timing.csv
./timing/plot_scan20_timing.sh                       # -> timing_runs/scan20_timing_*.png
```

- **Parity check:** compare the new `robust_ad` `SFmin = [...]` array (in each job log /
  merged `sfmin_scan.txt`) against the pre-change baseline already on disk:
  - baseline nb=32 robust_ad: `build/gacode_nb32_scan20_1node_robust_ad_54638544_tasks/sfmin_scan.txt`
    and the `SFmin = [...]` line in `build/time_scan20_nb32_julia_gpu_robust_ad_54638544.out`.
  - baselines nb=6/8/16: `build/time_scan20_nb{6,8,16}_julia_gpu_robust_ad_546351{90,91,92}.out`.
  - Require max relative diff ≤ 1e-6 per radius (ideally exact).
- **Speed check:** node-hours for `robust_ad` should drop vs the current line
  (1.26 nh @ nb=32), trending toward `:grid` MPS (~0.31 nh). Re-render the node-hours
  line plot and compare.

Baseline node-hours to beat (from `timing_runs/scan20_timing.csv`, current):

| N_BASIS | grid MPS | robust_ad (now) | truth (now) |
|---------|----------|-----------------|-------------|
| 6  | 0.20 | 0.18 | 0.21 |
| 8  | 0.21 | 0.20 | 0.24 |
| 16 | 0.22 | 0.42 | 0.48 |
| 32 | 0.31 | 1.26 | 1.38 |

(`robust_ad` already ≈ grid at nb=6/8; the big wins should appear at nb=16/32 where the
faithful core evals dominate.)

---

## 8. Acceptance criteria

1. Per-radius `robust_ad` `sfmin` unchanged vs baseline (≤1e-6 rel, exact preferred).
2. `:truth` `sfmin` unchanged (it floors on `robust_ad`'s `w≥1` min, so the snapshot
   `sfmin_w1` must remain the true `w≥1` grid min — verify `_rank_confirm` is used for
   the `w≥1` min that feeds `sfmin_w1`).
3. `robust_ad` node-hours at nb=16 and nb=32 measurably lower than the baseline table.
4. `confirm_grid=false` still reproduces the old brute-faithful behavior exactly.
5. Docs updated; commit + push.

---

## 9. Risks / fallbacks

- **Filter-heavy regimes (ITER rotational):** when secondary keep filters bind at many
  nodes, more confirms are needed and the speedup shrinks (degrades gracefully to the
  full faithful grid). `confirm_grid=false` is the explicit fallback.
- **`sparse` trigger heuristic** (§4.3): validate adaptive on/off parity; if a radius
  regresses, fall back to counting faithful feasibility only among confirmed nodes plus
  always-zoom-on-boundary, or set `adaptive=false` for that path.
- **`n_eig` choice:** keep `n_eig=24` for the core to match the exact `:confirm`; do not
  reduce without re-checking the bound holds (under-resolved cheap onset could over-prune).

---

## 10. How to start (fresh chat prompt)

> Read `docs/PLAN_robust_ad_cheap_confirm.md` and implement it. Start with §4
> (extract `_rank_confirm`, wire into `critical_factor_robust.eval_grid!` behind a
> `confirm_grid=true` kwarg, swap the `sparse` count to cheap feasibility, forward
> `confirm_grid` from `critical_factor_truth`). Then run §7.1 single-radius A/B
> exactness on CPU before any GPU scan. Only once `sfmin` parity holds, submit the §7.2
> 20-radius `robust_ad` node-hours sweep, check per-radius parity vs the on-disk
> baselines, re-collect/plot, and update the docs. Commit + push at the end.
