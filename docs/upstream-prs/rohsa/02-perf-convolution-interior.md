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
- Serial run time on the synthetic cube: 863 s against 888 s for `master`
  (2.8 % faster; mean of three runs, spread ±6 s). Both builds ran at the same
  time on four cores of one node, rotating cores between rounds, so they share
  the same load. The gain is small on its own. The main purpose is to make the
  convolution columns independent for #4.

