# ROHSA pull requests

Status: **submitted on 2026-10-19** from the fork `unsound-mint/ROHSA`:
#1 https://github.com/antoinemarchal/ROHSA/pull/4,
#2 https://github.com/antoinemarchal/ROHSA/pull/5,
#3 https://github.com/antoinemarchal/ROHSA/pull/6,
#4 https://github.com/antoinemarchal/ROHSA/pull/8,
#5 https://github.com/antoinemarchal/ROHSA/pull/9,
#6 https://github.com/antoinemarchal/ROHSA/pull/7,
#7 (submitted 2026-11-25) https://github.com/antoinemarchal/ROHSA/pull/10
(upstream numbers differ from ours; the PR texts were rewritten to the upstream
numbers). Upstream: https://github.com/antoinemarchal/ROHSA
(`master` at `e36bc49`, version 1.2.1; last commit February 2021, one
maintainer). Our vendored copy `rohsa/upstream/ROHSA` was identical to that
commit before our changes.

The work is split into seven PRs so that each one is small, has one purpose, and
says exactly how it changes the results:

| # | Branch | Patches | Base | What | Effect on results |
|---|---|---|---|---|---|
| 1 | `fix-lbfgsb-64bit-workspace` | 0001 | master | Fix the 32-bit overflow in the L-BFGS-B workspace (segfault on large maps) | bit-for-bit identical where it ran before |
| 2 | `perf-convolution-interior` | 0005 | master | Convolve only the interior points | bit-for-bit identical |
| 3 | `perf-cost-function-serial` | 0002-0004 | master | Gaussian once per channel, no redundant zeroing (bit-for-bit); per-Gaussian regularisation sum and integer exponents (rounding) | 3 commits, see PR |
| 4 | `openmp-cost-function` | 0002-0007 | master, needs #2 #3 | Optional OpenMP build (`make OPENMP=1`) for the cost function, one spectrum at a time; **justification for multi-core** | same result at any thread count; rounding vs #2 + #3 |
| 5 | `openmp-lbfgsb` | 0002-0009 | master, needs #4 | OpenMP in L-BFGS-B and BLAS, sums in fixed blocks | same result at any thread count; identical to #4 up to 65536 variables, rounding above |
| 6 | `binary-input` | 0010 | master | Optional binary input (`.bin`) + `ROHSApy.cube2bin` | text input unchanged; same values give identical output |
| 7 | `zero-padding` | 0011 | master | Set the padded part of `cube` and `std_cube` to 0 before `reshape_up` (upstream #10, 2026-11-25) | power-of-two maps unchanged; padded maps no longer depend on uninitialised memory |

`patches/` is one series, numbered in the order that applies to upstream
`e36bc49`. Patches 0001 to 0011 are the seven PRs (the Patches column; a stacked
PR contains the patches of the PRs it needs). 0011 is upstream
antoinemarchal/ROHSA#10 verbatim. 0012 (exact smoothing-term gradient at map
edges, stop-reason log) is not part of any PR yet. 0013 holds this project's
build defaults, plus the one comment line that points the vendored copy of 0011
at its upstream PR; it is not for upstream. 0014 is also local only: a separate
evaluation limit for the full-resolution level (`maxiter_last`) and parameter
checkpoints of that level at chosen evaluation counts (`checkpoint_at`), both off
by default, so a budget curve is read off one long run (test:
`test/rohsa_checkpoint.jl`). 0015 (`filename_range`) is a per-pixel live channel
range: channels outside it carry no data term and pixels with an empty range
none at all, for mosaics whose blocks cover different channels. It is off by
default and could become an upstream PR, but none is drafted (test:
`test/rohsa_live_range.jl`). Applying all fifteen gives
`rohsa/upstream/ROHSA` exactly; `scripts/check_rohsa_patches.sh` checks that
and runs in CI.

#1, #2, #3, #6 and #7 are independent and can be merged in any order. #4 contains
the commits of #2 and #3 and #5 those of #4 (stacked PRs); each says which
commit to review. #1 to #6 merge together without conflicts (branch
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
patches 0006 to 0009 are those versions, force-pushed to the PR branches the
same day. Two further commits that leave the output bit-for-bit
unchanged (`exp()` skipped where it is exactly zero, 0007; breakpoint list and
step bound of L-BFGS-B in parallel, 0009) were added the same evening; they are
in the patches but not yet on the PR branches.

## Submitting (after approval)

1. Fork `antoinemarchal/ROHSA` on GitHub (account of the person submitting).
2. `git am` the patches of each PR (Patches column) onto `master` in a clone of
   the fork, one branch per PR (or push the branches from our prepared clone).
3. Open #1, #2, #3, #6 against `master`. Open #4 after #2 and #3, and #5 after
   #4, each with a "depends on" line (GitHub cannot base a PR on another fork's
   PR branch).
4. Upstream has been quiet since 2021, so also e-mail the author (address in
   the upstream README) with the PR links and a two-line summary.

PR texts: `01-…md` to `06-…md` (title = first heading).
