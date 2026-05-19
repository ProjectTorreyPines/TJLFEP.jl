# TGLF-EP: Julia (TJLFEP) vs Fortran comparison results

This document records verification outcomes from the detailed comparison plan. Fortran sources live in `TJLFEP/TGLF-EP/`; Julia in `TJLFEP/src/`.

## 1. RMIN_LOC normalization (`verify-rmin-norm`)

| Code | Assignment |
|------|------------|
| Fortran `TGLFEP_read_EXPRO` | `rmin(:) = EXPRO_rmin(:) / a_meters` then `tglf_rmin_loc_in = rmin(ir)` |
| Julia `TJLF_map` | `RMIN_LOC = RMIN[ir] / RMIN[end]` |
| Julia FUSE `context.jl` | `RMIN` stored in metres; `RMIN_LOC = rmin[ir]/a` |

**Conclusion:** Equivalent for the IMAS/FUSE path when `RMIN` is in physical units. No change required for `RMIN_LOC`.

**Additional fix:** `GAMMA_THRESH_MAX` used raw `RMIN[ir]` in metres; now uses `r_over_a = RMIN[ir]/RMIN[end]` to match Fortran’s dimensionless `rmin(ir)`.

## 2. QL_flux_ratio for `n > nmodes_out` (`verify-ql-ratio-guard`)

| Code | Behaviour |
|------|-----------|
| Fortran `TGLFEP_ky` | Accumulates QL for all `nmodes`; no `1e20` guard |
| Julia (before) | Skipped QL accumulation when `n > nmodes_out`; set `QL_flux_ratio[n] = 1e20` |
| Julia (after) | QL accumulation matches Fortran; no `1e20` branch |

**Conclusion:** Julia aligned with Fortran. For `nmodes_out < NMODES`, empty QL arrays yield `NaN` ratios; Fortran and Julia both treat `NaN < thresh` as false (no QL-ratio rejection).

## 3. TGLF defaults (`verify-tglf-defaults`)

Julia `TJLFEP_ky` sets explicitly:

- `USE_AVE_ION_GRID = false`
- `FIND_EIGEN = true`
- `RLNP_CUTOFF = 18.0`
- `BETA_LOC = 0.0`, `KX0_LOC = 0.0`

Fortran `TGLFEP_ky` relies on GACODE `tglf_run` defaults after `TGLFEP_tglf_map`. TJLF.jl module defaults: `BETA_LOC=0`, `KX0_LOC=0`, `FIND_EIGEN=true` via `tjlf_read_input.jl`.

**Conclusion:** Consistent with TJLF defaults used by the Fortran driver path.

## 4. FACTOR_IN rounding (`verify-factor-rounding`)

Julia had `FACTOR_IN = round(100000*FACTOR_IN)/100000` in `TJLF_map`; Fortran does not round.

**Change:** Rounding removed for Fortran parity.

## 5. QL weight axis ordering (`verify-ql-axis-order`)

| Accessor | Index order |
|----------|-------------|
| Fortran `get_QL_particle_flux(n, species, jfields)` | mode, species, field |
| Julia `particle_QL_out[jfields, species, n]` | field, species, mode |

Julia indexing `particle_QL_out[jfields, inputsEP.IS_EP + 1, n]` matches Fortran `get_QL_*(n, is_EP+1, jfields)`.

**Conclusion:** Axis order is consistent.

## 6. IR_EXP mapping (`diagnose-crit-grad-diff`)

| Code | `INPUT_PROFILE_METHOD == 2` |
|------|------------------------------|
| Fortran | `ir_exp(i) = IRS + floor((i-1)*(NR-IRS)/(SCAN_N-1))` |
| Julia (before) | `argmin(|grid - rho[i]|)` |
| Julia (after) | Fortran linear formula restored in `run_tjlfep.jl` |

**Change:** `IR_EXP` now matches Fortran for DIII-D validation.

**Note:** For `NR=101`, `IRS=2`, `SCAN_N=20`, and uniform `ρ` grid, `argmin` and the Fortran formula yield the **same index vector** `[2, 7, 12, …, 101]`. Plot shifts from IR_EXP alone are unlikely for this case; remaining differences are dominated by `SFmin` (TJLF vs GACODE TGLF) and the 3- vs 4-species layout.

## 7. IS_EP / species layout

| | Fortran (DIII-D reference) | Julia IMAS path |
|--|---------------------------|-----------------|
| `input.TGLFEP` | `IS_EP=2`, `N_ION=2` → `ns=3` | `InputTGLFEP(...; is_ep=1)` → fast D on 1st ion |
| TJLF EP index | `is = IS_EP+1 = 3` | `Options.IS_EP = ep_slot-1 = 3`, `is = 4` with `ns=4` |

Julia keeps an extra thermal ion slot (`ns = N_ION + 2`). `is_ep=1` in `DIIID_juliaValidation.jl` selects fast D in IMAS; this corresponds to Fortran’s 2nd EXPRO ion (`IS_EP=2`) in the 3-species layout.

Post-processing uses `DENS_$ep_slot` (last species), consistent with the EP slot in the 4-species TJLF grid.

## 8. DIII-D setup (`setup-diiid-run`)

- Paths: `EP_TJLF` → `TJLFEP` in `DIIID_juliaValidation.jl` and `batchRun.sh`.
- TJLF branch: use `origin/gpu_new` at `/pscratch/sd/t/tneiser/.julia/dev/TJLF`.
- Sysimage: `TJLFEP/build/noTJLF_TJLFEP_sysimage.so` (optional).
- Comparison script: `src/DIIIDfiles/compare_fortran_julia.jl`.

## 10. Debug workflow (N_BASIS=6)

| Slurm job name | Code | Script |
|----------------|------|--------|
| `TGLFEP_nb6` | Fortran TGLF-EP | `build/batch_debug_nb6_fortran.sh` |
| `TJLFEP_nb6` | Julia TJLFEP | `build/batch_debug_nb6_julia.sh` |

Set `TGLFEP_DEBUG=1` (Fortran) and `TJLFEP_DEBUG=1` (Julia) for matched debug lines tagged `[TGLFEP_DBG]` and `[TJLFEP_DBG]`. Julia file runs should set `TJLFEP_FILE_ONLY=1` so `using TJLFEP` does not load IMAS/FUSE (unused on the string `runTHD` path). No sysimage is recommended while TJLFEP/TJLF are under active development.

## 11. Single-radius basis scans (SCAN_N=1, ir=2)

Same case as nb6; only `N_BASIS` changes in `build/debug_nb{N}/input.TGLFEP`.

| N_BASIS | Fortran job | Julia job | Compare script |
|---------|-------------|-----------|----------------|
| 6 | `TGLFEP_nb6` | `TJLFEP_nb6` | `compare_debug_nb6.sh` |
| 16 | `TGLFEP_nb16` | `TJLFEP_nb16` | `compare_debug_nb16.sh` |
| 32 | `TGLFEP_nb32` | `TJLFEP_nb32` | `compare_debug_nb32.sh` |

Nb6 single-radius match (2026-05-18): SFmin 2.8125, n_EP/n_e 96.1%, β_crit 4.46% (jobs 53156424 / 53156921).

Multi-radius production (`SCAN_N=20`) uses `build/debug_prod/input.TGLFEP` and `batch_prod_nb32_*` when needed.

## 12. N_BASIS=6 SCAN_N=20 on 10 nodes (2026-05-19)

| Job | Code | Script |
|-----|------|--------|
| `TGLFEP_nb6_s20_10n` | Fortran | `build/batch_debug_nb6_fortran_scan20_10n.sh` |
| `TJLFEP_nb6_s20_10n` | Julia | `build/batch_debug_nb6_julia_scan20_10n.sh` |

Compare: `build/compare_debug_nb6_scan20.sh` → `build/compare_nb6_scan20_plots/`.

At scan radii (`IR_EXP`): SFmin max rel err ~0.03%; α(dn/dr) ~0.5%; α(dp/dr) ~0.4%. Full reproduction steps: `docs/REPRODUCE_FORTRAN_MATCH.md`.

## 9. Running comparison (`test-agreement`)

```bash
cd /pscratch/sd/t/tneiser/.julia/dev/TJLFEP/src/DIIIDfiles
module load julia   # Perlmutter
julia --project=../.. compare_fortran_julia.jl path/to/GPU_n32_*_20_*
julia --project=../.. diagnose_crit_grad.jl path/to/GPU_n32_*_20_* 4
```

Fortran reference outputs: `202017C42_500ms_v3.1/out.scalefactor_r*`, `alpha_*_crit.input`, `out.TGLFEP`.

The comparison parser was smoke-tested on Fortran reference files (self-comparison). Re-run `DIIID_juliaValidation.jl` after the code fixes above, then point the scripts at the new `GPU_*` or `CPU_*` output directory.
