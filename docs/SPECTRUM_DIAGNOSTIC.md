# Spectrum diagnostic (`PROCESS_IN=3`)

Besides the critical-gradient **threshold** scans (`PROCESS_IN=4/5`), TJLFEP supports
the TGLF-EP **spectrum diagnostic** (`PROCESS_IN=3`). It does *not* scan or search for a
critical factor; it just computes the linear growth-rate / real-frequency spectra
γ(ky), ω(ky) at fixed plasma conditions, by running the full TGLF transport model over a
fixed `ky` grid (`nky=30`, `ky=0.15`, `KY_MODEL=0`, `NBASIS=32`) at three drives:

| `mode_in` | file | drive |
|-----------|------|-------|
| 1 | `out.eigenvalue_m1` | background plasma **+** EPs (thermal gradients retained) |
| 2 | `out.eigenvalue_m2` | EP drive only |
| 4 | `out.eigenvalue_m4` | EP-only with TAE/EPM filtered (typically leaves ITG/TEM) |

Comparing the three shows which part of the spectrum is EP-driven (AE/EPM) vs thermal
(ITG/TEM). With `WIDTH_IN_FLAG=true` the fixed `WIDTH_IN` is used; with
`WIDTH_IN_FLAG=false` a `TGLFEP_ky_widthscan` (EP-only) first picks the
max-growth-rate width (and marks `kymark`), then the spectrum is computed at that width.

Because it yields spectra rather than an `SFmin`/critical gradient, this mode is **not**
routed through `ActorTJLFEP` (which would error); call TJLFEP directly. A single radius
runs in a few minutes on a CPU login node; a full radial scan, a GPU run, or any long
job should go to Slurm (NERSC login nodes are shared and kill heavy/long processes).

A ready example input is `examples/DIIID_202017C42_500ms_v3.1/input_spectrum.TGLFEP`
(see also the manual `test/smoke_spectrum.jl`).

### Quick login-node run (one radius, CPU)

```bash
module load julia/1.11.7
export JULIA_DEPOT_PATH="$HOME/.julia:$PSCRATCH/.julia${JULIA_DEPOT_PATH:+:$JULIA_DEPOT_PATH}"
cd TJLFEP
julia --project=. -t 8 --startup-file=no
```

```julia
using TJLFEP   # loaded standalone, TJLFEP stays light (no IMAS/FUSE), CPU-only here

CASE   = "examples/DIIID_202017C42_500ms_v3.1"
gacode = joinpath(CASE, "input.gacode")
tglfep = joinpath(CASE, "input_spectrum.TGLFEP")   # PROCESS_IN=3

# run ONE radius (scan_index=1); writes out.eigenvalue_m{1,2,4} into out_dir
r = run_gacode_scan_task(gacode, tglfep, 1; out_dir="spectrum_out", printout=true)

r.spectra[1]   # (ky, gamma, freq) for mode_in=1; [2]=EP-only, [4]=ITG/TEM
```

`gamma`/`freq` are `nky × nmodes` matrices (row `i` = `ky[i]`). Useful threads scale up
to ~`nky=30` (the `ky` loop is threaded); keep it modest on a login node. For a few
radii, loop `scan_index` over `1:SCAN_N`.

### Slurm (full scan, or GPU)

Use a job when running many radii, on GPU, or unattended:

```bash
#!/bin/bash -l
#SBATCH -N 1 -n 1 -c 128 -C cpu -q debug -t 00:30:00 -A m3739 -J spectrum
module load julia/1.11.7
export JULIA_DEPOT_PATH="$HOME/.julia:$PSCRATCH/.julia${JULIA_DEPOT_PATH:+:$JULIA_DEPOT_PATH}"
cd TJLFEP
julia --project=. -t 32 --startup-file=no -e '
  using TJLFEP
  CASE="examples/DIIID_202017C42_500ms_v3.1"
  for i in 1:SCAN_N    # set to your input’s SCAN_N
      run_gacode_scan_task(joinpath(CASE,"input.gacode"), joinpath(CASE,"input_spectrum.TGLFEP"), i;
                           out_dir="spectrum_out", printout=true)
  end'
```

For GPU, use `-C gpu`, load `cudatoolkit`, and pass `use_gpu=true`. At `NBASIS=32` the
per-`ky` eigenmatrix is `(NS−NS0+1)·15·NBASIS = 1440²` (complex), so the GPU eigensolver
helps each solve, but the `ky` solves *serialize* on one GPU (no cross-`ky` batching),
so the real GPU payoff is across many radii via the same `inner=:threads`/`:mps_team`
multi-GPU layout used by the `PROCESS_IN=5` scans.

### Validation

`test/runtests_regression_spectrum.jl` compares the Julia spectra to a Fortran golden
(`test/fixtures/spectrum/out.eigenvalue_m{1,2,4}_r040`); agreement is ~3e-5 relative.
See `test/fixtures/spectrum/README.md` for how the golden was generated. Note the public
Fortran `TGLFEP_driver` has a bug on the `PROCESS_IN=3` + EXPRO path (`q_scale`
uninitialized → zeroed `q` → `DSYEV` failure), so a one-line-patched reference binary was
used.
