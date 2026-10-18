# Faster serial cost function

Branch `perf-cost-function-serial` · base `master` · independent of #2.
Three commits, so the bit-for-bit part can be reviewed separately from the two
that change rounding.

## 1. Evaluate each Gaussian once, drop redundant work (bit-for-bit identical)

In `f_g_cube_fast`:

- the Gaussian `exp()` and the weighted residual `residual/std²` are computed
  once per channel and shared by the three derivatives (they were evaluated
  three times), and the `std_map` test moves outside the channel loop;
- the residual is written straight into the `residual` cube instead of a
  temporary array allocated and freed for every spectrum;
- `params`, `residual`, `g` and the convolution/image work arrays were zeroed
  although they are overwritten before use. Only `deriv` and the last `n_gauss`
  entries of `g` (the only part `ravel_3D` does not write) are zeroed now.

Same expressions in the same order: identical output.

## 2. Sum the regularisation term per Gaussian (rounding only)

The smoothness and width terms of each Gaussian are summed in a local variable
and added to `f` once per Gaussian, instead of term by term. This is what lets
the OpenMP version in #4 give exactly the serial result on one thread. The
objective changes only by the grouping of a floating-point sum; the gradient is
unchanged.

## 3. Integer exponents (rounding only)

`x**2._xp` and `x**3._xp` are real powers; gfortran evaluates many of them with
`pow()`, which was about a fifth of the run time in a profile. `x**2` and `x**3`
compile to multiplications. Mathematically identical; results differ by the
rounding of `pow()` versus `x*x`.

## Checks

On the synthetic 256 × 256 × 100 cube of `check/make_synthetic.py` (N = 4, one thread):

- Commit 1: bit-for-bit identical to `master`.
- Whole PR vs `master`: per-pixel integrated model differs by median 3.2e-5,
  99th percentile 6.6e-4 (relative). For scale: changing the input values
  only in their last bit (≤ 5.6e-17) moves the same quantity by median
  5.8e-5, so this is rounding-level, not a change of the solution.
- Serial run time on the synthetic cube: 518 s against 888 s for `master`,
  **1.71× faster** (mean of three runs, spread ±3 s, same setup as #2). With
  #1, #2, #3 and #6 combined: 491 s (1.81×).

