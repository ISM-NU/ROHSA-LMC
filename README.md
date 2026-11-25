# ROHSA-LMC

Method repo for the ROHSA Gaussian decomposition of H I cubes, as used on the
LMC sightline. It wraps the upstream Fortran code ROHSA (Marchal et al. 2019,
A&A 626, A101; vendored in `rohsa/upstream/ROHSA` with OpenMP patches) in a
Julia package, `RohsaLMC`. The package prepares inputs, runs ROHSA, derives
FITS products from its fitted parameters, builds velocity-layer cubes and
de-rotates cubes by the disk velocity. It does not reimplement the optimizer.
Experiments that use it (scans, nulls, injections, cold-gas maps) live in
the group's LMC-cold-gas repository, which stays private until publication.

## Setup

Group members usually work from the LMC workspace (the private ISM-NU/LMC-workspace repository), whose `setup.sh` clones the related repositories side by side and links the shared `lmc-data` folder.

`deps/VGT` is a private repository, so the clone below needs access to the ISM-NU organisation. Without it `--recursive` stops with "Repository not found" and the package cannot be built.

```bash
git clone --recursive https://github.com/ISM-NU/ROHSA-LMC.git
cd ROHSA-LMC
ln -s ../lmc-data data          # shared data folder next to the repo, ignored by git
git submodule foreach --recursive 'ln -sfn "$toplevel/data" data'   # same data inside submodules
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

`deps/VGT` is a git submodule (ISM-NU/VGT.jl).
It provides the Julia package `VGT` (FITS I/O, noise estimates, VGT).
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
- `scripts/` - stand-alone evaluation scripts (`scripts/README.md`).
- `test/` - tests.
- `docs/` - validation notes, upstream pull requests, the ROHSA paper. Map: `docs/README.md`.
- `deps/VGT` - submodule, the VGT method.

## Related repos

- [VGT.jl](https://github.com/ISM-NU/VGT.jl) (private, ISM-NU members only) - the VGT method (`VGT` package). Pinned here as `deps/VGT`.
- LMC-foreground (private until publication) - Milky Way foreground, Planck and starlight comparisons. Experiments.
- LMC-cold-gas (private until publication) - ROHSA runs and cold H I work. Pins this repo.
- [LMC-shells](https://github.com/ISM-NU/LMC-shells) - optical supergiant-shell research. Experiments.

## License

MIT (`LICENSE`) for the code in this repository outside `rohsa/upstream/ROHSA`.
The vendored ROHSA source is Copyright (c) 2019 Antoine Marchal, MIT license
(`rohsa/upstream/ROHSA/LICENSE`); its bundled L-BFGS-B keeps its BSD-3-Clause
notice. Work using ROHSA should cite Marchal et al. 2019, A&A 626, A101.
