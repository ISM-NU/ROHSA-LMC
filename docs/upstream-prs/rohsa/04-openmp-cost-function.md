# Optional OpenMP parallelisation of the cost function

Branch `openmp-cost-function` · base `master` · **depends on #2 and #3**
(their commits are included; please review the last two commits,
"Optional OpenMP parallelisation of the cost function" and "Skip exp() where it
is exactly zero").

## Summary

`make OPENMP=1` builds ROHSA with OpenMP, so one decomposition can use several
cores of one machine. Every sum is taken in a fixed order, so the result is
bit-for-bit the same at any `OMP_NUM_THREADS` and in the serial build.

Updated 2026-11-13: the first version used OpenMP reductions, so the result
changed with the number of threads. This version removes that and also makes
the serial code faster (one spectrum at a time, no residual cube).

## Why go from one core to several

**Survey cubes have outgrown a single core.** ROHSA fits every pixel of a map
together (the spatial regularisation couples neighbours, and the multigrid
descent refines the whole map level by level), so a large map cannot be split
into independent pieces without changing the result. For the GASKAP-HI cube of
the Large Magellanic Cloud we decompose a 2048 × 2048 × 159 region: 75 million
parameters with N = 6, 101 million with N = 8.

**Measured speed-up (first version).** We timed a 512 × 512 × 159 window of that cube (N = 6,
λ′σ = 30). Each run had its own allocation of exactly that many cores on one
AMD EPYC 9684X node. The nodes were shared with other jobs, so expect about
±10 % per number. Wall time, and speed-up over one thread of the same build:

| threads | this PR | with #5 |
|---:|---:|---:|
| 1 | 5207 s (1.0×) | 5466 s (1.0×) |
| 2 | 3084 s (1.7×) | 2874 s (1.9×) |
| 4 | 1569 s (3.3×) | 1295 s (4.2×) |
| 8 | 1139 s (4.6×) | 809 s (6.8×) |
| 16 | 1199 s (4.3×) | 665 s (8.2×) |
| 32 | 955 s (5.5×) | 419 s (13.0×) |

This PR alone levels off at about 5×: the cost function is then fast, and the
serial L-BFGS-B update takes most of the time. #5 parallelises that part. The
full 2048 × 2048 × 159 region (N = 8) took 2 h 58 min on 32 threads with all
PRs. At the one-thread rate above, it would take roughly 13 times as long, about
1.6 days. That is an extrapolation from the window; we did not make a serial run
of the full region.

**This version against the first one** (both with #5, same node, N = 16):
512 × 512 × 159 at 32 threads, 724 s → 539 s, peak memory 4.4 → 3.9 GiB;
128 × 128 × 159 at one thread, 452–519 s → 325–368 s; 2048 × 2048 × 159 for
10 iterations at 32 threads, 323 s → 260 s.

**Choosing hyper-parameters needs many runs.** N, λ and λ′σ have to be chosen
by comparing decompositions (the paper's recommendation, and what we do on
smaller windows first). At a day per run on one core, a small grid of settings
on one map takes weeks; with threads it fits in a few days on one node.

**Why OpenMP and not the alternatives:**

- *Splitting the map into tiles and running them in parallel* changes the
  result at the tile edges, because the regularisation couples neighbouring
  pixels. Threads compute exactly the same cost function and gradient.
- *MPI* would need a domain decomposition of the same coupled problem and a
  distributed L-BFGS-B; much more code for a single-node gain.
- *GPU* would mean rewriting the cost function and the optimizer.
- OpenMP needs only directives on loops that are already independent (one
  spectrum or one sky column per iteration), keeps one address space, and
  disappears completely from a serial build.

**What it costs:**

- *Reproducibility*: the result does not depend on the number of threads.
  Against #2 + #3 the serial result changes by rounding only: the derivatives
  reuse the `exp()` values of the model (gfortran vectorises that loop, and the
  vector `exp()` can differ from the scalar one in the last bit), and the
  regularisation sum is grouped by sky column. With `exp()` recomputed as
  before, the cost function is bitwise identical to #2 + #3.
- *Memory*: lower than before. The residual cube is gone, and no large array is
  copied per thread (each thread keeps one spectrum's model and Gaussians).
  f is stored per pixel (one map) so that it can be added in serial order.
- *Maintenance*: directives on loops over sky columns, one new routine
  (`data_term`), no new dependencies (OpenMP ships with gfortran).

## What runs in parallel

All loops over independent sky columns:

1. `data_term` for each spectrum: model, residual, the spectrum's share of f
   and its derivatives, written straight into `g`; the shares of f are then
   added in pixel order. `exp()` is not called more than 39σ + 2 channels from
  a line centre, where its argument is below -760 and the result is exactly 0
  (the vector `exp()` falls back to slow scalar code there; output unchanged);
2. per Gaussian, the extraction of the three parameter maps from `beta` and the
   regularisation loop (each pixel writes only its own gradient entries; the
   regularisation parts of f and of the `b_params` gradient are summed per
   column, then over the columns in order);
3. `convolution_2D_mirror`: zero fill, interior copy and the convolution loop
   (rewritten in #2 so that columns are independent).

## How we checked it

- Serial build and `OPENMP=1` at 1, 4 and 8 threads: bit-for-bit identical
  (48 × 48 × 159 window of the LMC cube, N = 8; 512 × 512 × 159 window, N = 16,
  at 16 and 32 threads).
- With `exp()` recomputed in the derivative loop as before, the serial build is
  bit-for-bit identical to #2 + #3, so the rounding differences come only from
  the two sources above.
- Equivalence with #2 + #3 under rounding: we moved every input value by one
  float32 ulp at random, 24 times, on two 128 × 128 × 159 windows, and ran both
  versions. Pairwise differences between the two versions are the same size as
  within each version (permutation test p ≥ 0.78), and no summary statistic
  (final cost, total flux, flux by line width, components per pixel) differs
  (smallest Mann–Whitney p = 0.06, 18 tests). Some runs of both versions end in
  a worse local minimum, at similar rates.

The first version's timings and thread comparisons are in
`check/results-timing-2026-10-19.txt`.

## How to use

```
cd src && make clean && make OPENMP=1
OMP_NUM_THREADS=16 ./ROHSA parameters.txt
```

`OMP_NUM_THREADS` should be set explicitly; unset, OpenMP uses every core of
the machine. A README section is added.

