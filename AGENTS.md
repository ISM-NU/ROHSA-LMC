# AGENTS.md

## What this repo is

Method repo: the `RohsaLMC` Julia package, the ROHSA runner and CLIs, the
vendored Fortran ROHSA with OpenMP patches, and the Slurm scripts that run them.
Read `README.md` first, then `rohsa/README.md`.

## Scope

- Belongs here: code that prepares, runs or post-processes a ROHSA decomposition, and its tests, validation notes and cluster scripts.
- Does not belong here: experiments (scans, null and injection tests, cold-gas maps, figures, comparisons with stars or Planck). They go to `LMC-cold-gas` (ROHSA experiments) or `LMC-foreground` (foreground work).
- VGT itself is not changed here. It comes from the sibling checkout `../VGT.jl` (ISM-NU/VGT.jl); change it there, not here.

## Language

Use Julia where performance needs it. Simple scripts can be Python.

## Science rules

- ROHSA defaults (`lambda` terms, `n_gauss`, bounds, noise handling, iteration limits, smoothing) are science settings. A change needs a stated reason and a validation, and goes through the docs.
- Do not change what a product means (amplitude, width, centre, layer split) without a documented decision. Experiment repos pin this repo and compare runs across versions.
- Keep parity checks. The output must not depend on the number of threads (`test/rohsa_threads.jl`). A change that alters rounding must pass the ensemble test in "Equivalence of results" of `docs/validation/rohsa-deterministic-openmp-2026-11-13.md`. The 2026-10-17 note is the first single-run check, kept as history.
- `rohsa/upstream/ROHSA` is the upstream base (ab4f299) plus the series in `docs/upstream-prs/rohsa/patches`. If you change the vendored Fortran, add a patch at the end of the series and update the text in `docs/upstream-prs/rohsa`. `scripts/check_rohsa_patches.sh` fails on any difference and runs in CI.
- Do not treat a smoother or prettier decomposition as a better one.

## Data and Slurm

- Data live in `~/Projects/LMC/lmc-data`. `data` in the repo is a symlink made with `ln -s ../lmc-data data`. It is ignored by git. Never commit cubes or run outputs.
- Heavy computation runs on Irgetas via Slurm, not on laptops. Seconds-long scripts and tests may run locally.
- The repo root on the cluster is `ROHSA_REPO`. Runbook: `hpc/irgetas/README.md`. Read it before submitting.
- Do not remove run outputs unless the run finished successfully or the user asks.

## Tests and docs

- Run `julia --project=. test/runtests.jl` before committing code changes.
- Add a test only for a durable behavior: a contract, an invariant, a data format, a real regression.
- Docs are short and have a status (active guidance, validation evidence, historical, archive). Add new docs to `docs/README.md`.
- Papers go under `docs/references/` and are summarized in the docs, not copied.

## VGT dependency

VGT comes from the sibling checkout `../VGT.jl` (`Manifest.toml` has `path = "../VGT.jl"`). There is no pin: a change pushed to VGT.jl main is used here on the next run. Run this repo's tests after pulling VGT.jl.
