# ROHSA: deterministic OpenMP and fused data term (2026-11-13)

Status: validation evidence. Supersedes the thread-count statements of
`rohsa-openmp-parity-2026-10-17.md`.

## Change

Vendored `rohsa/upstream/ROHSA/src`, upstream PR branches
`openmp-cost-function` and `openmp-lbfgsb` (v2):

- `f_g_cube_fast` works one spectrum at a time (`data_term`): each Gaussian is
  evaluated once per channel and reused for the model and the three
  derivatives, so there is no residual cube. Parameters are read from `beta` and
  the gradient is written to `g` in place.
- Each spectrum's share of f is stored and added in serial order; the
  regularisation sums are taken per sky column, then over the columns in order.
- L-BFGS-B and BLAS sums over the variables are cut into blocks of 65536
  (`omp_block.h`), each summed in index order, block sums added in order.

Every sum therefore has a fixed order: the output is bit-for-bit the same at any
`OMP_NUM_THREADS` and in the serial build. Against the previous code the result
changes by rounding only (vectorised `exp` values reused in the derivatives,
regularisation grouped by column, blocked sums above 65536 variables). With
`exp` recomputed as before, the new cost function is bitwise identical to the
old one on one thread.

Guard: `test/rohsa_threads.jl` builds the sources and checks that 1 and 3 threads
give identical output with n above one block (fails on the previous code).

## Speed and memory

Same node type (AMD EPYC 9684X), GASKAP-HI LMC tile, N = 16, `noclip`,
smoothing 10, external noise map.

| Case | Previous | New |
|---|---|---|
| 512 × 512 × 159, 400 iterations, 32 threads (same node, back to back) | 724 s, 4.44 GiB | 539 s, 3.94 GiB |
| same, 16 threads | | 552 s, output identical to 32 threads |
| 128 × 128 × 159, 400 iterations, 1 thread (13 runs each, 13 in parallel) | 452–519 s | 325–368 s |
| 2048 × 2048 × 159, 10 iterations, 32 threads (same node) | 323 s | 260 s, peak 64 GiB |

`rohsa/memory.jl` no longer counts the residual and derivative cubes; its
estimate (4.1 GiB for the 512² tile, 65.6 GiB for 2048²) stays above the
measured peaks.

## Equivalence of results

The previous code is not reproducible run to run (thread-dependent summation
order), so the new code is compared with the spread of the previous one under
input perturbation: every voxel moved by ±1 float32 ulp at random, 128² sub-crops
of the tile, all runs at 1 thread.

1. First ensemble (13 runs per build, one crop). Pre-set check (a), "new run on
   the unperturbed input is no farther from the previous unperturbed run than
   the farthest perturbed previous run", **failed** on the centroid p99
   (0.945 vs 0.917 channels; integrated-model p99 passed). Check (b), no
   ensemble statistic differs (Mann–Whitney p ≥ 0.01), passed. Recomputing the
   final cost showed why: the new unperturbed run had converged to a worse local
   minimum (cost 2.6 % above every other run). Both builds produce such runs
   (1 of 13 previous, 2 of 13 new), in the same part of the crop. A single run
   against a 12-run envelope fails about once in 13 by chance, so this check
   could not decide.
2. Second ensemble, criterion written into the analysis script before the runs:
   24 fresh perturbations per build on two crops (the first crop again and a
   second sky area). Tests: permutation test on mean cross-build minus mean
   within-build pairwise distance (integrated model p99, centroid p99), and
   Mann–Whitney on final cost, total flux, flux share by width class
   (σ < 2, 2–4, 4–8, ≥ 8 channels) and Gaussians per pixel. PASS if all 18 tests
   survive Holm correction at family α = 0.05. **Result: PASS**, smallest
   p = 0.06; permutation p ≥ 0.78. Within-build and cross-build distances agree
   (centroid p99 medians 0.35 / 0.37 / 0.36 channels on crop 1). Runs whose
   final cost is more than 1 % above the crop median: 2/24 previous vs 3/24 new
   on crop 1, 8/24 vs 11/24 on crop 2 (95 % intervals overlap widely).

On the 512² tile the new run lies inside the spread of three earlier production
replicates of the previous code (per-pixel integrated-model p99 0.054–0.069
against them; replicate pairs 0.043–0.075).

## Consequences

- Repeating a run no longer gives a replicate. Replicates must come from
  perturbed inputs (for example ±1 ulp) or different settings.
- A bump of the ROHSA-LMC pin changes ROHSA results by rounding; say so in the
  commit of the dependent repo.
- Some fraction of runs on a crop ends in a worse local minimum under either
  build; compare final costs before reading a single run.
