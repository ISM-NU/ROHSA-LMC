# OpenMP in L-BFGS-B and BLAS

Branch `openmp-lbfgsb` · base `master` · **depends on #4** (please review the
last commit, "OpenMP in L-BFGS-B and BLAS (results independent of the thread count)").

## Summary

With the cost function parallel (#4), about half of the remaining run time at
12 threads was spent in the serial parts of L-BFGS-B, which loop over all n
parameters several times per iteration. In an `OPENMP=1` build these loops now
run in parallel: `cauchy`, `cmprlb`, `formk`, `subsm`, `projgr`, the line-search
update in `lnsrlb`, and `ddot`, `dcopy`, `daxpy` and `dscal` in `blas.f`.

Updated 2026-11-13: the first version used OpenMP reductions, so the result
changed with the number of threads. This version does not.

## Reproducibility

Sums over the variables are cut into blocks of `nblk` = 65536 elements
(`omp_block.h`). Each block is summed in index order and the block sums are
added in block order, so the result is the same for any number of threads and
with or without OpenMP.

Up to 65536 variables there is one block and every sum is the same as in #4.
For more variables the order of the additions changes, which changes results
at the rounding level.

Loops that update one element per variable run in parallel as they are. In
`cmprlb` and `subsm` the loop over the variables is now the outer loop, which
keeps the order of the updates of each element. The breakpoint list in
`cauchy`, the step bound in `lnsrlb` and the backtracking in `subsm` stay serial
because their result depends on the order.

## Checks

- n = 82 956 (above one block): the serial build and `OPENMP=1` at 1, 4 and 8
  threads give bit-for-bit the same output. With n below one block the output
  is bit-for-bit that of #4.
- Real data, both PRs together, against the first version (same node, N = 16):
  512 × 512 × 159 at 32 threads, 724 s → 539 s, and 16 threads give the same
  output as 32; 2048 × 2048 × 159 (n = 201 million) for 10 iterations at
  32 threads, 323 s → 260 s.
- Equivalence with the previous summation order under ±1-ulp input
  perturbations: see #4.
