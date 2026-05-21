# Reproducing Fortran–Julia TGLF-EP match (DIII-D 202017C42)

Validated on Perlmutter CPU (premium QoS). See `FORTRAN_JULIA_COMPARISON.md` for physics fixes.

## Repositories

```bash
# Sibling checkouts (Project.toml uses path dep ../TJLF)
git clone git@github.com:ProjectTorreyPines/TJLFEP.jl.git TJLFEP
cd TJLFEP && git checkout fortran_match

git clone git@github.com:ProjectTorreyPines/TJLF.jl.git ../TJLF
cd ../TJLF && git checkout gpu_new   # or branch named in your site setup
```

## Julia environment

```bash
module load julia/1.11.7
export JULIA_DEPOT_PATH=$PSCRATCH/.julia   # optional, site-specific

cd TJLFEP
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
```

File-based runs (no IMAS/FUSE):

```bash
export TJLFEP_FILE_ONLY=1
export GACODE_DUMP=$PWD/src/DIIIDfiles/202017C42_500ms_v3.1/input.gacode
```

## Fortran TGLF-EP

Build `TGLFEP_driver` from GACODE add-on `TGLF-EP` (Perlmutter CPU). Edit paths in batch scripts if needed:

```bash
export GACODE_ROOT=.../gacode
export GACODE_ADD_ROOT=.../gacode_add
export GACODE_PLATFORM=PERLMUTTER_CPU
```

## Reference case

`src/DIIIDfiles/202017C42_500ms_v3.1/` (`dump.profile`, `input.gacode`).

## SCAN_N=20 on GPU (5 nodes, gacode-only)

Recommended: **one job on 5 GPU nodes** — `srun -n 20`, **4 tasks/node**, **1 A100 + 32 CPU threads** per radius.

```bash
cd build
./submit_gacode_scan20_gpu.sh
# or manually:
sbatch batch_run_gacode_scan20_gpu_5nodes.sh
```

Scan and merge run in one job (`srun` then `finalize_gacode_scan` on the batch head node CPU).

Alternative (20 separate array tasks): `sbatch batch_run_gacode_scan20_gpu_array.sh`, then
`OUT_DIR=build/gacode_scan20_<ARRAY_JOB_ID>_tasks sbatch batch_merge_gacode_scan20.sh`.

Inputs: `input.gacode` + `build/debug_nb6/input_scan20.TGLFEP` (`N_BASIS=6`, `SCAN_N=20`).

Outputs per task: `task_<i>.jls`, optional `out.scalefactor_r*`. After merge: `alpha_*_crit.input`, `sfmin_scan.txt`.

## N_BASIS=6, SCAN_N=20, 10 nodes (CPU, legacy file inputs)

| Step | Command |
|------|---------|
| Fortran | `cd build && sbatch batch_debug_nb6_fortran_scan20_10n.sh` |
| Julia | `cd build && sbatch batch_debug_nb6_julia_scan20_10n.sh` |
| Compare | `FORTRAN_DIR=fortran_runs/debug_nb6_scan20_10n_<FJOB> JULIA_DIR=debug_out_nb6_scan20_<JJOB>_dist FILE_DIR=debug_nb6/fileInput_scan20_10n_<JJOB> ./compare_debug_nb6_scan20.sh` |

Inputs: `build/debug_nb6/input_scan20.TGLFEP` (`N_BASIS=6`, `SCAN_N=20`, `IRS=2`).

Expected agreement at `IR_EXP` radii (2026-05-19 jobs 53171364 / 53171385):

- SFmin: max relative error ~0.03%
- α(dn/dr), α(dp/dr): max relative error ~0.5%

Plots: `build/compare_nb6_scan20_plots/` (created by compare script; gitignored).

## Other scripts

| N_BASIS | Single radius | SCAN_N=20 (1 node) | SCAN_N=20 (10 nodes) |
|---------|---------------|--------------------|----------------------|
| 6 | `batch_debug_nb6_{fortran,julia}.sh` | `batch_debug_nb6_*_scan20.sh` | `batch_debug_nb6_*_scan20_10n.sh` |
| 16 | `batch_debug_nb16_*` | — | — |
| 32 | `batch_debug_nb32_*` | — | `batch_prod_nb32_*` |

Text diff of outputs: `src/DIIIDfiles/compare_fortran_julia.jl` with `FORTRAN_DIR` / `JULIA_DIR` set.
