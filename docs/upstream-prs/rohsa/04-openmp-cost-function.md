# Optional OpenMP parallelisation of the cost function

Depends on #5 and #6, whose commits are included. Review the last two commits: "Optional OpenMP parallelisation of
the cost function" and "Skip exp() where it is exactly zero".

## What this does

`make OPENMP=1` builds ROHSA with OpenMP, so one decomposition can use all cores of one machine. Every sum runs in a
fixed order, so the output is bit-for-bit the same for any `OMP_NUM_THREADS` and in the serial build. Without
`OPENMP=1` the compiler ignores the directives and the build is serial.

## Changed since the first version

The first version used OpenMP reductions, so the output changed with the number of threads. This version sums in a
fixed order instead. The serial code is also faster now: it fits one spectrum at a time and no longer stores a
residual cube.

## Why

ROHSA fits every pixel of a map together, because the smoothing term couples neighbouring pixels. A large map
cannot be split into tiles without changing the result at the tile edges. Our GASKAP-HI region of the Large
Magellanic Cloud is 2048 × 2048 × 159 voxels, which is 101 million parameters with N = 8. With this PR and #9 it
takes 2 h 58 min on 32 threads. At the one-thread rate it would take about 1.6 days. That figure is extrapolated
from a 512 × 512 window, not measured. Choosing N and the λ values needs a grid of such runs, so one core is too
slow.

OpenMP fits the problem because the loops are already independent per spectrum or per sky column. MPI would need a
distributed L-BFGS-B. A GPU port would mean rewriting the cost function and the optimiser.

## What runs in parallel

All loops over independent sky columns run in parallel:

- `data_term`, a new routine, computes the model, the residual, the share of f and the derivatives of one
  spectrum. It writes the derivatives straight into `g`. The shares of f are then added in pixel order.
- For each Gaussian, the extraction of the three parameter maps from `beta` and the regularisation loop. Each pixel
  writes only its own gradient entries. The regularisation parts of f are summed per column, then over the columns
  in order.
- `convolution_2D_mirror`, which #5 rewrote so that columns are independent.

`exp()` is not called more than 39σ + 2 channels from a line centre. There its argument is below −760 and the
result is exactly 0, so the output does not change.

## Differences from #5 and #6

The serial output differs from #5 + #6 at the rounding level, for two reasons:

- The derivatives reuse the `exp()` values of the model. gfortran vectorises that loop, and the vector `exp()` can
  differ from the scalar one in the last bit.
- The regularisation sum is grouped by sky column.

With `exp()` recomputed as before, the cost function is bitwise identical to #5 + #6.

Memory use is lower than before. The residual cube is gone, and each thread keeps only the model of one spectrum.

## How we checked it

- The serial build and `OPENMP=1` at 1, 4 and 8 threads give bit-for-bit the same output on a 48 × 48 × 159 window
  of the LMC cube with N = 8. The same holds at 16 and 32 threads on a 512 × 512 × 159 window with N = 16.
- We tested whether the rounding differences from #5 + #6 matter. We moved every input value by one float32 ulp at
  random, 24 times, on two 128 × 128 × 159 windows, and ran both versions. The differences between the two
  versions are the same size as the differences within each version (permutation test, p ≥ 0.78). No summary
  statistic differs over 18 tests (final cost, total flux, flux by line width, and components per pixel). The
  smallest Mann–Whitney p is 0.06.
- With #9 on one node (N = 16), 512 × 512 × 159 at 32 threads went from 724 s with the first version to 539 s.
  Peak memory went from 4.4 to 3.9 GiB. The timings of the first version are in
  `check/results-timing-2026-10-19.txt`.

## How to use

```
cd src && make clean && make OPENMP=1
OMP_NUM_THREADS=16 ./ROHSA parameters.txt
```

Set `OMP_NUM_THREADS`. If it is unset, OpenMP uses every core of the machine. The README has a new section on
this.
