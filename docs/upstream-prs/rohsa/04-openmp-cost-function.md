# Optional OpenMP parallelisation of the cost function

Branch `openmp-cost-function` · base `master` · **depends on #2 and #3**
(their commits are included; please review the last commit,
"Optional OpenMP parallelisation of the cost function").

## Summary

`make OPENMP=1` builds ROHSA with OpenMP, so one decomposition can use several
cores of one machine. The default `make` is unchanged: without `-fopenmp` the
directives are comments and the code is the serial code of #2 and #3.

## Why go from one core to several

**Survey cubes have outgrown a single core.** ROHSA fits every pixel of a map
together (the spatial regularisation couples neighbours, and the multigrid
descent refines the whole map level by level), so a large map cannot be split
into independent pieces without changing the result. For the GASKAP-HI cube of
the Large Magellanic Cloud we decompose a 2048 × 2048 × 159 region: 75 million
parameters with N = 6, 101 million with N = 8.

TODO-SCALING-NUMBERS

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

- *No change unless asked for*: the default build is serial and identical to
  #2 + #3.
- *Reproducibility*: with `OMP_NUM_THREADS=1` an OpenMP build gives bit-for-bit
  the serial result. With more threads the objective is summed in a different
  order, so results differ by rounding only: TODO-MULTITHREAD-DIFF. For
  scale: changing the input values in their last bit only (≤ 5.6e-17) already
  moves the same quantity by median 5.8e-5 — the optimizer amplifies
  rounding-level differences, so a different summation order is the same
  kind of difference, not a different solution.
- *Memory*: no large array is copied per thread. The private variables are
  scalars and one spectrum-sized model array; peak memory is the same as the
  serial build.
- *Maintenance*: 4 directive pairs in two routines, no new dependencies
  (OpenMP ships with gfortran).

## What runs in parallel

All loops over independent sky columns:

1. the model/residual/objective loop of `f_g_cube_fast` (one spectrum per
   pixel; the objective is a reduction);
2. the gradient and regularisation loop for each Gaussian (each pixel writes
   only its own derivatives; the per-Gaussian sums from #3 are reductions,
   which is why #3 introduced them);
3. the zeroing of `deriv`;
4. `convolution_2D_mirror` (rewritten in #2 so that columns are independent).

## How we checked it

On the synthetic 256 × 256 × 100 cube of `check/make_synthetic.py` (N = 4, one thread):

- serial build of this PR: bit-for-bit identical to #2 + #3 (serial);
- `OPENMP=1` build at `OMP_NUM_THREADS=1`: bit-for-bit identical to the serial
  build.

TODO-MULTITHREAD-CHECKS

## How to use

```
cd src && make clean && make OPENMP=1
OMP_NUM_THREADS=16 ./ROHSA parameters.txt
```

`OMP_NUM_THREADS` should be set explicitly; unset, OpenMP uses every core of
the machine. A README section is added.

