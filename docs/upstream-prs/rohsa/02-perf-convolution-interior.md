# Speed up convolution_2D_mirror

Branch `perf-convolution-interior` · base `master` · independent.

## Summary

`convolution_2D_mirror` is called six times per Gaussian in every evaluation of
the cost function. It convolved the whole padded (ny+4) × (nx+4) array, with a
bounds test on each of the 9 kernel terms, allocated a second padded array on
every call, and then kept only the interior ny × nx points.

It now computes only the interior points, directly into `conv`. Their
neighbours are always inside the padded array, so the bounds test and the
second array are not needed. The kernel terms are added in the same order, so
**the result is bit-for-bit identical**.

This is also what makes the convolution parallel-friendly (each output column
is independent), used by #4.

## Checks

- On the synthetic 256 × 256 × 100 cube of `check/make_synthetic.py` (N = 4, one thread): decomposition bit-for-bit identical to `master`.
- TODO-PR2-SPEED (serial run time vs `master`).

