# AGENTS.md

## What this repo is

Method repo: the `RohsaLMC` Julia package, the ROHSA runner and CLIs, the
vendored Fortran ROHSA with OpenMP patches, and the Slurm scripts that run them.
Read `README.md` first, then `rohsa/README.md`.

## Scope

- Belongs here: code that prepares, runs or post-processes a ROHSA decomposition, and its tests, validation notes and cluster scripts.
- Does not belong here: experiments (scans, null and injection tests, cold-gas maps, figures, comparisons with stars or Planck). They go to `LMC-cold-gas` (ROHSA experiments) or `LMC-foreground` (foreground work).
- VGT itself is not changed here. It comes from `deps/VGT` (ISM-NU/VGT-Denis). Do not edit anything under `deps/`.

## Language

Use Julia where performance needs it. Simple scripts can be Python.

## Science rules

- ROHSA defaults (`lambda` terms, `n_gauss`, bounds, noise handling, iteration limits, smoothing) are science settings. A change needs a stated reason and a validation, and goes through the docs.
- Do not change what a product means (amplitude, width, centre, layer split) without a documented decision. Experiment repos pin this repo and compare runs across versions.
- Keep parity checks: a speed change must reproduce the previous fit within the tolerance in `docs/validation/rohsa-openmp-parity-2026-10-17.md`.
- Keep the vendored Fortran in sync with `docs/upstream-prs/rohsa`. If you patch `rohsa/upstream/ROHSA`, update the matching patch and text there, and the reverse.
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

## Submodule

After changing `deps/VGT` pointers in a consumer repo, commit the pointer there. Here, update with `git -C deps/VGT pull` and commit the new pointer.
