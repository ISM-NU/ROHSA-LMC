# ROHSA OpenMP patch: parity check (2026-10-17)

Status: validation evidence.

Change: `rohsa/upstream/ROHSA/src/mod_optimize.f90` (`f_g_cube_fast`) parallelizes
the per-pixel residual and gradient loops with OpenMP (`-fopenmp` in the
Makefile) and evaluates the Gaussian `exp` once per channel instead of three
times. The expressions are otherwise unchanged. `hpc/irgetas/rohsa.slurm`
requests 16 CPUs and sets `OMP_NUM_THREADS` from Slurm; without it OpenMP would
start one thread per node core inside a 1-CPU allocation.

Test: 30 Dor 512x512x159 window (`configs/rohsa/dor30_512_n6_defaultnoise.toml` inputs,
noise=false, N=6, lambda=100), original binary vs patched binary, same
parameter file.

| Metric | Original vs patched |
|---|---|
| Wall time (laptop) | 84.7 min (1 core) vs 24.5 min (12 threads, ~580% CPU) |
| Per-pixel integrated model | median rel. diff 0.26%, p99 1.9% |
| Per-pixel model centroid | median 0.045 km/s, p99 0.24 km/s |
| Total model emission | +0.16% |
| Column share sigma <2 / 2-6.4 / >6.4 km/s | 0.02/23.0/76.9% vs 0.02/23.3/76.7% |

Results are not bitwise identical: the OpenMP reduction changes the summation
order of the objective, and L-BFGS-B then follows a slightly different path.
Individual Gaussians can swap or shift (the decomposition is not unique), but
the reconstructed model and derived statistics agree well within the
differences between hyper-parameter settings. Speed-up is limited by the
serial parts (L-BFGS-B vector updates, convolutions); profile those before
expecting more from extra threads.

## Follow-up: convolution and integer powers (same day)

Profiling (gprof single-thread; stack sampling of the main thread at 12
threads) showed the remaining serial time in `convolution_2D_mirror`
(full padded-array pass with a bounds test per term, two allocations per call,
single-threaded), in `pow()` from real-valued exponents (`x**2._xp`) in the
gradient, and in the L-BFGS-B internals (`cauchy`, `formk`, `subsm`, `ddot`).

- `convolution_2D_mirror` now computes only the interior points, adds the
  terms in the same order, and runs columns in parallel. Output is bitwise
  identical to the previous binary (64² crop, 1 thread).
- `mod_optimize.f90` uses integer exponents (`**2`, `**3`), which compile to
  multiplications instead of `pow()` calls. Not bitwise identical (rounding
  of `pow` vs `x*x`); per-pixel integrated model vs previous build: median
  9e-5, p99 3.5e-3 relative.

256² crop around 30 Dor, N=6, external noise, 12 threads (laptop):

| Build | Wall time |
|---|---|
| OpenMP gradient only | 371 s |
| + interior-only parallel convolution | 339 s |
| + integer powers | 256 s |
| + parallel L-BFGS-B `cauchy` main loop, parallel long `ddot`/`dcopy` | 240 s |
| + no redundant zero-fills in `f_g_cube_fast` | 214 s |
| + parallel `subsm`, `formk`, `cmprlb`, `projgr`, long `dscal` | 169 s |

The `cauchy` split (parallel classification and sums, serial breakpoint list
in index order) and the zero-fill removal are bitwise identical at 1 thread
(64² crop). At 12 threads the per-pixel integrated model differs from the
previous build by median 1e-4, p99 5e-3 relative (summation order).

## Round 3: rest of L-BFGS-B, binary input, 64-bit workspace

- `subsm`, `cmprlb` and `formk` loop over the free variables once and
  accumulate all 2·m dot products of a pass at the same time (OpenMP array
  reductions); the per-variable updates swap loop order so each variable is
  updated by one thread in the original order. `projgr` (max reduction) and
  `dscal` above 100000 elements run in parallel. Bitwise identical to the
  previous build at 1 thread (256² crop, full run); at 12 threads, per-pixel
  integrated model vs previous build median 1.6e-4, p99 5.5e-3 relative, the
  same size as run-to-run summation-order differences. CPU use at 12 threads
  893% -> 1107%.
- ROHSA reads a binary input (`<crop>.bin`: int32 nv, ny, nx, then Float64 in
  `cube(v,y,x)` order) when the file name ends in `.bin`; the Julia wrapper
  writes it whenever the cube is in memory. Same values in text and binary give
  bitwise identical output (256² crop, 1 thread). Parsing the 20 GB text input
  of a 2048² crop takes ROHSA about 4 min (0.36 µs per line, measured on a
  256² cube) and writing it ~7 min. (An earlier estimate of ~15 min came from
  a misattributed profiler sample.)
- L-BFGS-B workspace offsets are 64-bit (`setulb`, `mod_minimize`): with
  2·m·n above 2³¹ (N = 8 on 2048², m = 10) the default-integer offsets wrapped
  and `cauchy` segfaulted (Irgetas jobs 31682, 31683). Bitwise identical
  (30 Dor 512², 1 thread).

Remaining main-thread time: the breakpoint heap in `cauchy` (`hpsolb`,
inherently sequential), `freev` (ordered index lists) and the small dense
algebra on 2m×2m matrices.
