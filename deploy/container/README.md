# TJLFEP container

Self-contained TJLFEP images on `ghcr.io/projecttorreypines/tjlfep` for
**standalone** file-based runs (`input.TGLFEP` / `input.gacode` workflows,
`run_gacode_scan_task` / `runTHD_from_gacode`), with CPU **and** GPU capability
in one image. For the full-FUSE use case, use the
[FUSE container](https://github.com/ProjectTorreyPines/FUSE.jl/tree/master/deploy/perlmutter-container)
instead; for multi-node `SlurmClusterManager` scans, use the bare-metal CFS
install (see Limitations).

## Tags

| Tag | Contents | Built by |
|---|---|---|
| `:<v>`, `:latest` | lean: TJLFEP + TJLF + CUDA, pkgimages only | GitHub Actions (`.github/workflows/container.yml`) |
| `:<v>-imas`, `:latest-imas` | + IMAS/GACODE/TurbulentTransport (`TJLFEPIMASExt`) | GitHub Actions |
| `:<v>-gpu` | lean + GPU-traced multi-target sysimage | Manual, Perlmutter A100 (`gpu-sysimage/`) |

`latest` never points at a `-gpu` tag. All images are linux/amd64, CPU targets
`generic;znver3;znver4;sapphirerapids` (Perlmutter Milan, Defiant-era hosts,
generic fallback). Every image carries the **CUDA 12.9 runtime artifact**
(pinned via `CUDA_Runtime_jll` Preferences, CUDA.jl 5.8.5) and works fully
offline; without an injected GPU driver, `CUDA.functional()` is `false` and
everything falls back to the CPU path — one image serves both.

There is no file-only flag: the lean image simply lacks the IMAS trigger
packages, so the `TJLFEPIMASExt` extension never loads; in the `-imas` image it
loads automatically.

## NERSC Perlmutter (podman-hpc)

One-command setup:

```bash
bash <(curl -fsSL https://raw.githubusercontent.com/ProjectTorreyPines/TJLF-EP.jl/master/deploy/container/install_tjlfep_container_nersc.sh)
```

or by hand:

```bash
podman-hpc pull ghcr.io/projecttorreypines/tjlfep:v2.0.14
podman-hpc tag ghcr.io/projecttorreypines/tjlfep:v2.0.14 localhost/tjlfep:v2.0.14
podman-hpc migrate tjlfep:v2.0.14

# CPU (any node)
podman-hpc run --rm -it localhost/tjlfep:v2.0.14

# GPU (inside a GPU allocation; --gpu injects the host driver)
podman-hpc run --rm --gpu -it localhost/tjlfep:v2.0.14
```

The in-image `tjlfep` launcher starts Julia against the baked project
(`/opt/tjlfep`), auto-loads the GPU sysimage when present (`-gpu` tags), and
prepends a writable `$HOME/.julia_tjlfep_container` depot when the image
depot is read-only.

Batch templates: `test_slurm_cpu.sbatch`, `test_slurm_gpu.sbatch` (the GPU one
also documents CUDA MPS setup for SPMD team runs). Acceptance suite:
`acceptance.sh` (CPU, login-node safe; see header for options).

## ORNL Defiant (Apptainer) — prospective

```bash
apptainer pull tjlfep.sif docker://ghcr.io/projecttorreypines/tjlfep:v2.0.14-gpu
apptainer exec --nv tjlfep.sif tjlfep -e 'import CUDA; @show CUDA.functional()'
```

Caveats (untested until access exists — please report back):

- `--nv` injects the host driver; the image runs its 12.9 runtime via CUDA
  minor-version compatibility on newer drivers (CUDA.jl warns, proceeds). If a
  host driver rejects the 12.9 runtime, re-pin in your writable home depot:
  `tjlfep -e 'import CUDA; CUDA.set_runtime_version!(v"12.x")'` — the launcher's
  home-depot fallback makes the new preference stick.
- SIF images are read-only: the launcher automatically prepends
  `$HOME/.julia_tjlfep_container` as the first (writable) depot.

## Building

```bash
./deploy/container/build.sh                     # lean, local podman-hpc build
TJLFEP_VARIANT=imas ./deploy/container/build.sh # imas variant
TJLFEP_ENVIRONMENT=v2.0.14 ./deploy/container/acceptance.sh
TJLFEP_ENVIRONMENT=v2.0.14 ./deploy/container/publish_ghcr.sh
```

CI (`workflow_dispatch` or release) builds and pushes the lean + imas tags.
The version resolves from FuseRegistry's `Versions.toml` (this repo publishes
no GitHub releases). **First push only:** the GHCR package is created private —
an org owner must flip it public for unauthenticated pulls.

### GPU sysimage bake (`:<v>-gpu`)

CI cannot bake the sysimage (the trace executes real GPU solves). On a
Perlmutter login node:

```bash
TJLFEP_ENVIRONMENT=v2.0.14 ./deploy/container/gpu-sysimage/bake_and_publish.sh
```

pulls + migrates the CI image, runs the bake on an A100 (sbatch, m5377_g
premium, ~1–4 h), layers `sys_tjlfep_gpu.so` in with `Containerfile.gpu`
(COPY layer — reproducible, no `podman-hpc commit`), runs CPU acceptance, and
submits GPU acceptance. Publish after the GPU job passes:

```bash
TJLFEP_ENVIRONMENT=v2.0.14 TAG_SUFFIX=-gpu ./deploy/container/publish_ghcr.sh
```

## Limitations (v1)

- **Single-node only**: threads + local `addprocs` + CUDA MPS. Multi-node
  `SlurmClusterManager` runs would need every srun-spawned worker to start
  in-container — use the bare-metal CFS install for those.
- linux/amd64 only (no arm64/multi-arch manifests).
- The GPU sysimage is baked manually per release; plain tags run the same GPU
  code via JIT (slower first solve, identical throughput after).
