# ROHSA on NU Irgetas

How we run the ROHSA decomposition, the velocity-layer cubes and the
de-rotation on the Irgetas cluster. Login node: use it only to copy files and
submit jobs. VGT runs on the layer cubes use the VGT repo's own Slurm script
(`deps/VGT/hpc/irgetas/vgt.slurm`, see below).

| Script | What it does |
|---|---|
| `submit_rohsa.sh <config>` | Preflight a ROHSA config, size the job from its crop and N, submit `rohsa.slurm` |
| `rohsa.slurm` | Rebuild ROHSA if its sources changed, run `rohsa/cli/preprocess.jl` with one OpenMP thread per allocated CPU |
| `rohsa_preflight.jl <config>` | Read the cube dimensions and estimate peak memory |
| `rohsa_layers.slurm` | Build low/main/high/model layer cubes and their VGT configs from a finished ROHSA run |
| `derotate.slurm` | De-rotate cubes by the disk velocity and write a VGT config for each |

Every script reads the repository root from `ROHSA_REPO` (set by
`submit_rohsa.sh`; pass it with `--export` for the others). Run configs for
experiments live in LMC-cold-gas; method examples are in `configs/rohsa/`.

## Access

Connect to the NU GlobalProtect VPN when off campus, then
`ssh YOUR_USERNAME@<irgetas login node>`. Our copy of the repository is
`/shared/home/YOUR_USERNAME/rohsa-lmc`. It is not a git checkout: sync tracked
files from a laptop checkout with

```bash
rsync -az --files-from=<(git ls-files --recurse-submodules | grep -v '^rohsa/upstream/ROHSA/src/ROHSA$') \
  ./ YOUR_USERNAME@<irgetas login node>:/shared/home/YOUR_USERNAME/rohsa-lmc/
```

The ROHSA binary is excluded on purpose: `rohsa.slurm` rebuilds it on the
cluster whenever a Fortran source is newer than the binary. Large cubes go in
`data/` with `rsync -az --partial`, which resumes after a VPN drop.

`module` only exists in a login shell: use `bash -lc "module load GCC && ..."`
over ssh, and `#!/bin/bash -l` in ad-hoc batch scripts. Julia packages live in
the depot `.julia-irgetas` inside the repository (`JULIA_DEPOT_PATH`); set it up
once with `julia --project=. -e 'using Pkg; Pkg.instantiate()'`.

## ROHSA runs

```bash
bash hpc/irgetas/submit_rohsa.sh configs/rohsa/full_cube.toml
squeue --me
```

The job is named after the config. Logs:
`validation_output/rohsa_irgetas_logs/<config name>-<job id>.out|.err`.
Results: the config's `output_dir`.

ROHSA here is our OpenMP build of the upstream Fortran code
(`rohsa/upstream/ROHSA/src`; changes and parity checks in
`docs/validation/rohsa-openmp-parity-2026-10-17.md`). About half of the runtime
is still single-threaded optimizer bookkeeping, so beyond about 32 cores more
CPUs mostly sit idle.

### Resources

`submit_rohsa.sh` picks these from `crop_size` and `n_gauss`; override with
`CPUS=`, `MEM=` or `TIME=` in front of the command.

| Run | CPUs | Measured peak memory | Measured wall time | Requested |
|---|---|---|---|---|
| 512² crop, N = 6-8 | 16 | 5-6 GB | 20-31 min | 12 GB, 3 h |
| 1024² crop | 16 | not yet measured with OpenMP | | 40 GB, 12 h |
| 2048² crop, N = 6 | 32 | job 31677 | ~6 h | 100 GB, 1 day |
| 2048² crop, N = 8 | 32 | 28 GB before the last level | job 31697 | 140 GB, 1.5 days |
| whole cube (`crop_size = 0`) | 32 | not run | | 180 GB, 3 days |

Check a finished job with
`sacct -j JOB_ID -o JobID,State,Elapsed,AllocCPUS,ReqMem,MaxRSS` and tighten the
table when a peak is known. Long time limits do not slow a run, but the
scheduler starts short jobs sooner. `TotalCPU` stays near zero because ROHSA
runs as a child of Julia; it is not a CPU-efficiency measure here.

The optimizer workspace grows as 25 × n × 8 bytes, where n = 3 × N × pixels
(N = 8 at 2048²: n = 100.7 M, 20 GB). Its offsets are 64-bit, so there is no
size limit other than memory.

### Disk

A run writes its maps, `rohsa_n<N>.dat` (the fitted parameters) and a FITS copy
of the crop. Cubes are not written by default and ROHSA's input copy is deleted
after success. Home quota is 400 GB.

## Velocity layers

After a ROHSA run has finished:

```bash
sbatch --export=ALL,ROHSA_REPO=$PWD,REGION_DIR=<run dir>/<region>,N_GAUSS=<N>,LAYER_DIR=<out dir> \
  hpc/irgetas/rohsa_layers.slurm
```

This writes `<region>_{model,low,main,high}.fits`, the layer fraction maps, and
one VGT config per cube plus `vgt_data.toml` for the original data (the
control). Pass extra `rohsa/cli/layer_cubes.jl` options (for example
`--split-by width --cold-max-fwhm 4.8 --warm-min-fwhm 4.8`) with
`LAYER_ARGS="..."` in `--export`.

Run VGT on each cube with the VGT repo's script, from the VGT checkout
(`deps/VGT`), whose environment variable is `VGT_REPO`:

```bash
sbatch --job-name=vgt_<cube> --export=ALL,VGT_REPO=$PWD,VGT_CONFIG=<LAYER_DIR>/vgt_<cube>.toml \
  deps/VGT/hpc/irgetas/vgt.slurm
```

Each config already has its own `output_dir`.

## De-rotation

```bash
sbatch --export=ALL,ROHSA_REPO=$PWD,REFERENCE=<reference_velocity_channels.fits>,OUT_DIR=<out>,CONFIG_DIR=<dir with vgt_<name>.toml>,CUBES="name=cube.fits ..." \
  hpc/irgetas/derotate.slurm
```

`reference_velocity_channels.fits` comes from the layer step. Options of the
underlying CLI are in the header of `rohsa/cli/derotate_cube.jl`.

## Citations

Work using ROHSA cites Marchal et al. 2019, A&A 626, A101. NU asks users to
acknowledge Nazarbayev University Research Computing when its HPC facilities
contribute materially to published results.
