# OpenMP in L-BFGS-B and BLAS

Branch `openmp-lbfgsb` · base `master` · **depends on #4** (please review the
last commit, "OpenMP in L-BFGS-B: cauchy, subsm, formk, cmprlb, projgr, BLAS").

## Summary

With the cost function parallel (#4), about half of the remaining run time at
12 threads was spent in the serial parts of L-BFGS-B, which loop over all n
parameters several times per iteration. In an `OPENMP=1` build these loops now
run in parallel:

| Routine | Change |
|---|---|
| `cauchy` | The main loop is split into two passes. The first classifies each variable and accumulates `f1` and `p` (parallel); the second builds the breakpoint and free-variable lists in index order (serial), so `iorder`, `t`, `nbreak`, `nfree`, `ibkmin` and `bnded` are exactly those of the original loop. |
| `subsm`, `cmprlb`, `formk` | The 2m dot products with the stored correction pairs are accumulated in **one** pass over the free variables (array reductions, each sum still in index order) instead of 2m passes; the per-variable updates are reordered so each variable is updated by one thread with its terms in the original order. |
| `projgr` | Max reduction. |
| `ddot`, `dcopy`, `dscal` | Long unit-stride vectors (n > 100 000) in parallel when more than one thread is available. All new lines use the `!$` sentinel, so a serial build compiles the original BLAS code, and with one thread the original unrolled `ddot` is used. |

The breakpoint heap (`hpsolb`), `freev` and the small 2m × 2m algebra stay
serial.

## Reproducibility

In a serial build, and in an OpenMP build with `OMP_NUM_THREADS=1`, every sum is
formed in the original order: **bit-for-bit identical** to #4. With more
threads the dot products are summed in a different order (rounding only).

## Checks

- On the synthetic 256 × 256 × 100 cube of `check/make_synthetic.py` (N = 4, one thread): the serial build and the `OPENMP=1` build at
  `OMP_NUM_THREADS=1` are both bit-for-bit identical to #4's serial build.
- TODO-PR5-SPEED (run time at several thread counts vs #4).

