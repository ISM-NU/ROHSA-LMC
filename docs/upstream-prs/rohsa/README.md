# ROHSA pull requests

Status: **submitted on 2026-10-19** from the fork `unsound-mint/ROHSA`:
#1 https://github.com/antoinemarchal/ROHSA/pull/4,
#2 https://github.com/antoinemarchal/ROHSA/pull/5,
#3 https://github.com/antoinemarchal/ROHSA/pull/6,
#4 https://github.com/antoinemarchal/ROHSA/pull/8,
#5 https://github.com/antoinemarchal/ROHSA/pull/9,
#6 https://github.com/antoinemarchal/ROHSA/pull/7
(upstream numbers differ from ours; the PR texts were rewritten to the upstream
numbers). Upstream: https://github.com/antoinemarchal/ROHSA
(`master` at `e36bc49`, version 1.2.1; last commit February 2021, one
maintainer). Our vendored copy `rohsa/upstream/ROHSA` was identical to that
commit before our changes.

The work is split into six PRs so that each one is small, has one purpose, and
says exactly how it changes the results:

| # | Branch | Base | What | Effect on results |
|---|---|---|---|---|
| 1 | `fix-lbfgsb-64bit-workspace` | master | Fix the 32-bit overflow in the L-BFGS-B workspace (segfault on large maps) | bit-for-bit identical where it ran before |
| 2 | `perf-convolution-interior` | master | Convolve only the interior points | bit-for-bit identical |
| 3 | `perf-cost-function-serial` | master | Gaussian once per channel, no redundant zeroing (bit-for-bit); per-Gaussian regularisation sum and integer exponents (rounding) | 3 commits, see PR |
| 4 | `openmp-cost-function` | master, needs #2 #3 | Optional OpenMP build (`make OPENMP=1`) for the cost function, one spectrum at a time; **justification for multi-core** | same result at any thread count; rounding vs #2 + #3 |
| 5 | `openmp-lbfgsb` | master, needs #4 | OpenMP in L-BFGS-B and BLAS, sums in fixed blocks | same result at any thread count; identical to #4 up to 65536 variables, rounding above |
| 6 | `binary-input` | master | Optional binary input (`.bin`) + `ROHSApy.cube2bin` | text input unchanged; same values give identical output |

#1, #2, #3 and #6 are independent and can be merged in any order. #4 contains
the commits of #2 and #3 and #5 those of #4 (stacked PRs); each says which
commit to review. All six merge together without conflicts (branch
`all-prs-combined`), and that combination is what our pipeline runs.

## Evidence

`check/` holds everything needed to reproduce the claims without our data:

- `make_synthetic.py`: a 256 × 256 × 100 synthetic cube (three spatially
  coherent Gaussians plus noise), written with ROHSApy's `cube2dat` and the new
  `cube2bin`.
- `params_template.txt`, `run_one.sh`, `compare.py`: run a build on it and
  compare outputs (bit-for-bit, otherwise per-pixel integrated-model
  differences).
- `results-*.txt`: our runs.

Timings on real data use windows of the GASKAP-HI LMC cube on 192-core AMD
EPYC (Zen 4) nodes of the NU Irgetas cluster. The 2026-11-13 versions of #4 and
#5 (sums in a fixed order) were checked as described in
`../../validation/rohsa-deterministic-openmp-2026-11-13.md`; the patches in
`patches/openmp-*` are those versions, force-pushed to the PR branches the
same day. Two further commits that leave the output bit-for-bit
unchanged (`exp()` skipped where it is exactly zero; breakpoint list and step
bound of L-BFGS-B in parallel) were added the same evening; they are in the
patches but not yet on the PR branches.

## Submitting (after approval)

1. Fork `antoinemarchal/ROHSA` on GitHub (account of the person submitting).
2. `git am` the patches in `patches/<branch>/` onto `master` in a clone of the
   fork, one branch per PR (or push the branches from our prepared clone).
3. Open #1, #2, #3, #6 against `master`. Open #4 after #2 and #3, and #5 after
   #4, each with a "depends on" line (GitHub cannot base a PR on another fork's
   PR branch).
4. Upstream has been quiet since 2021, so also e-mail the author (address in
   the upstream README) with the PR links and a two-line summary.

PR texts: `01-…md` to `06-…md` (title = first heading).
