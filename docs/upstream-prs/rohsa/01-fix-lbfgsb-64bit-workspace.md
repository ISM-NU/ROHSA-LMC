# Fix integer overflow in the L-BFGS-B workspace for large problems

Branch `fix-lbfgsb-64bit-workspace` · base `master` · independent of the other PRs.

## Problem

ROHSA crashes with a segmentation fault at the start of the last multigrid
level for large maps:

```
Program received signal SIGSEGV: Segmentation fault - invalid memory reference.
#1  in cauchy_   at lbfgsb.f:1473
#2  in mainlb_   at lbfgsb.f:682
#3  in setulb_   at lbfgsb.f:280
#4  in __mod_minimize_MOD_minimize   at mod_minimize.f90:143
```

`setulb` keeps the offsets of its sub-arrays in `wa` (`ws`, `wy`, …, `d`, `t`,
`xp`, `wa`) in default 32-bit integers. The offset of `d` is about
`n·(2m + 2)`. With `n = 3·N·nx·ny` parameters this passes 2³¹ − 1 =
2 147 483 647 for large maps: N = 8 on 2048 × 2048 gives n = 100.7 million and,
with m = 10, an offset of 2.2 billion. The offset wraps to a negative number and
`cauchy` writes outside `wa`. The array size passed to `allocate` in
`mod_minimize` overflows the same way. (N = 6 on the same map, 1.66 billion,
still fits, which is why the limit is easy to miss.)

## Change

- The offsets in `setulb` are `integer*8`, computed from `n` and `m` on every
  call (instead of being stored in `isave(1:16)` on the first call).
- The `wa` dummy argument of `setulb` is assumed-size (`wa(*)`), so its
  declared size is not computed in 32-bit arithmetic either.
- `mod_minimize` allocates `wa` with a 64-bit size expression.

Nothing else in L-BFGS-B depends on these offsets; `isave(22:44)`, which callers
read, is unchanged.

## Checks

- On the synthetic 256 × 256 × 100 cube of `check/make_synthetic.py` (N = 4, one thread): output bit-for-bit identical to `master`
  (`check/results-bitwise-2026-10-17.txt`).
- The run that crashed (N = 8, 2048 × 2048 × 159, m = 10) passes the start of
  the last level with the fix: TODO-PR1-BIGRUN.

