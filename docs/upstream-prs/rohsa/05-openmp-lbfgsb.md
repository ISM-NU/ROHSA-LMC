# OpenMP in L-BFGS-B and BLAS

Depends on #8. Review the last two commits: "OpenMP in L-BFGS-B and BLAS (results independent of the thread count)"
and "L-BFGS-B: breakpoint list and step bound in parallel (same result)".

## What this does

With #8 the cost function runs in parallel. At 12 threads about half of the remaining time then goes to the serial
loops of L-BFGS-B. In an `OPENMP=1` build these loops now run in parallel: `cauchy`, `cmprlb`, `formk`, `subsm`,
`projgr`, the line-search update in `lnsrlb`, and `ddot`, `dcopy`, `daxpy` and `dscal` in `blas.f`. The output does
not depend on the number of threads.

## Changed since the first version

The first version used OpenMP reductions, so the output changed with the number of threads. This version does not.

## How the output stays the same

Sums over the variables are cut into blocks of `nblk` = 65536 elements, set in `omp_block.h`. Each block is summed
in index order, and the block sums are added in block order. The result is therefore the same for any number of
threads, with or without OpenMP. Up to 65536 variables there is one block, and every sum matches #8. With more
variables the order of the additions changes, which changes results at the rounding level.

Loops that update one element per variable run in parallel as they are. In `cmprlb` and `subsm` the loop over the
variables is now the outer loop, so each element is updated in the same order as before. `cauchy` builds its
breakpoint list in two passes over the blocks. The first pass counts, and the second fills at the counted offsets,
so the list comes out in the same order as before. The step bound in `lnsrlb` is a minimum and a zero test, and
neither depends on the order. The breakpoint heap in `hpsolb` and the backtracking in `subsm` stay serial.

## How we checked it

- With n = 82 956, which is more than one block, the serial build and `OPENMP=1` at 1, 4 and 8 threads give
  bit-for-bit the same output. With n below one block the output is bit-for-bit that of #8.
- The last commit and the `exp()` change in #8 do not change the output. Together they save a further 12 to 15 %.
  In paired runs on one node, 100 iterations on 512 × 512 × 159 at 32 threads went from 118–122 s to 103 s.
- With #8 on real data (N = 16, one node), 512 × 512 × 159 at 32 threads went from 724 s with the first version to
  539 s, and 16 threads give the same output as 32. On 2048 × 2048 × 159 (201 million variables), 10 iterations at
  32 threads went from 323 s to 260 s.
- #8 describes the test of the new summation order under random one-ulp changes of the input.
