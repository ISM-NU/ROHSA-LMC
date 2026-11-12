# ROHSA-LMC

Method repo for the ROHSA Gaussian decomposition of H I cubes, as used on the
LMC sightline. It wraps the upstream Fortran code ROHSA (Marchal et al. 2019,
A&A 626, A101; vendored in `rohsa/upstream/ROHSA` with OpenMP patches) in a
Julia package, `RohsaLMC`. The package prepares inputs, runs ROHSA, derives
FITS products from its fitted parameters, builds velocity-layer cubes and
de-rotates cubes by the disk velocity. It does not reimplement the optimizer.
Experiments that use it (scans, nulls, injections, cold-gas maps) live in
[LMC-cold-gas](https://github.com/ISM-NU/LMC-cold-gas), not here.

## Setup

Usually this repository sits in the LMC workspace folder (`~/Projects/LMC`), whose `setup.sh` clones all related repositories next to each other and links the shared data; its AGENTS.md has the rules that span repositories.

```bash
git clone --recursive git@github.com:ISM-NU/ROHSA-LMC.git
cd ROHSA-LMC
ln -s ../lmc-data data          # shared data folder ~/Projects/LMC/lmc-data, ignored by git
git submodule foreach --recursive 'ln -sfn "$toplevel/data" data'   # same data inside submodules
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

`deps/VGT` is a git submodule (ISM-NU/VGT-Denis, branch `split/vgt-method`).
It provides the Julia package `Astronomy` (FITS I/O, noise estimates, VGT).
Update the pin with `git -C deps/VGT pull`, then commit the new pointer.

The Fortran binary is built on first use by the runner (needs `gfortran`). On
Irgetas, `rohsa.slurm` rebuilds it when a source is newer.

## Run

From the repository root:

```bash
julia --project=. test/runtests.jl                       # tests
julia --project=. rohsa/cli/preprocess.jl                # decomposition, default configs/rohsa/full_cube.toml
julia --project=. rohsa/cli/layer_cubes.jl <region_dir> <n_gauss> <out_dir> [options]
julia --project=. rohsa/cli/derotate_cube.jl <reference.fits> <out_dir> <config_dir> name=cube.fits ...
```

Each CLI prints its options in the header comment of the file. The full-cube
config is rejected on a laptop unless `force_full_cube = true`. Use
`configs/rohsa/full_cube_laptop_stride8.toml` for a small local run. Heavy runs
go to the cluster: see `hpc/irgetas/README.md`.

ROHSA workflow details: `rohsa/README.md`.

## Layout

- `src/` - package `RohsaLMC`: `rohsa.jl` (parameters, noise maps, products), `rohsa_layers.jl` (layer cubes), `derotate.jl`.
- `rohsa/` - runner, config loading, memory preflight, report writers, CLIs in `rohsa/cli/`.
- `rohsa/upstream/ROHSA/` - vendored Fortran source with OpenMP patches.
- `configs/rohsa/` - method example configs. Experiment configs live in LMC-cold-gas.
- `hpc/irgetas/` - Slurm scripts and runbook.
- `test/` - tests.
- `docs/` - validation notes, upstream pull requests, the ROHSA paper. Map: `docs/README.md`.
- `deps/VGT` - submodule, the VGT method.

## Related repos

- [VGT-Denis](https://github.com/ISM-NU/VGT-Denis) - the VGT method (`Astronomy` package). Pinned here as `deps/VGT`.
- [LMC-foreground](https://github.com/ISM-NU/LMC-foreground) - Milky Way foreground, Planck and starlight comparisons. Experiments.
- [LMC-cold-gas](https://github.com/ISM-NU/LMC-cold-gas) - ROHSA runs and cold H I work. Pins this repo.
- [LMC-shells](https://github.com/ISM-NU/LMC-shells) - optical supergiant-shell research. Experiments.
