Status: validation evidence (synthetic data). Opt-in; the default output is unchanged.

# Live channel range per pixel, pixels without data, floored per-pixel noise (2026-11-22)

Why: survey cubes made as mosaics of observing blocks can give each pixel its own contiguous range of observed
channels, with fill values (exact 0) outside it and NaN beyond the survey edge. The wrapper used to write NaN as 0 and
fit every channel, so fill channels and pixels beyond the edge were fitted as measured zeros, and one noise value
weighted blocks with different noise. ROHSA's cost had no per-channel weight.

## What changed

- Fortran (`rohsa/upstream/ROHSA/src`, opt-in through the namelist key `filename_range`):
  - `mod_array.f90`: the switch and the two integer maps `range_lo`, `range_hi` (no cube-sized weight array).
  - `mod_inout.f90`: `read_range` (header `ny nx`, rows `y x first last`, 0-based; `first > last` = no data).
  - `mod_optimize.f90` `data_term`: residuals outside `kf..kl` are set to 0 before f and the gradient are formed.
    On the full range this is a no-op.
  - `mod_functions.f90`: `mean_array_live` and `mean_map_live` for the coarse levels; in `update`, a pixel without
    data gets the level's largest amplitude bound instead of `maxval` of its zero spectrum (which pinned its
    amplitudes at 0, so the smoothing term pulled its neighbours towards 0).
  - `mod_rohsa.f90`: reads the file, builds each level's ranges, stops if combined with `lym`.
- Coarse-level rule: union (hull) of the children's ranges, each channel averaged over the children observed in it;
  coarse noise = mean over children with data and noise > 0. Intersection was rejected: it drops a block's extra
  channels wherever a coarse pixel straddles a block boundary and, at the coarsest levels, cuts the band down to the
  narrowest block. Dividing by all children (the upstream `mean_array`) dilutes partly covered coarse pixels.
- Wrapper: `missing_data = "skip"` writes the ranges from the raw crop and gives pixels without data noise 0;
  `noise_floor_quantile` floors the per-pixel line-free map. Both off by default. The line-free estimator already
  skipped exact zeros and NaN; a test now covers a zeroed channel block.

## Checks

1. Byte identity, one-off (2026-11-26, laptop, gfortran 16.2, 4 OpenMP threads), when the change was ported onto the
   RohsaConfig code as patch 0015: the build before it (5c68855) and the new build without a range file, N = 3, on a
   16 x 16 x 60 and a padded 13 x 20 x 60 synthetic cube, with the exact smoothing gradient on and off, a per-pixel
   noise map, ROHSA's own noise and a run to convergence (526 evaluations). out.dat and every `lbfgsb_stop` line are
   identical in all seven cases.
2. Byte identity, permanent (`test/rohsa_live_range.jl`): no range file vs an all-live range file on a 32 x 32 map,
   for a flat and a varying noise map and for ROHSA's own noise estimate. On a map that is not a power of two the two
   differ by design: with a range file the padding has an empty range (no data) instead of being fitted as zeros, so
   the coarse levels average over the real pixels only.
3. Recovery (`test/rohsa_live_range.jl`): 32 x 32 x 80 cube, two Gaussians (3 K at channel 30 + 0.15 x, σ 3; 1.5 K at
   50 - 0.1 y, σ 2), noise 0.1 K, N = 3, λ = 10, 200 evaluations. Channels 1-28 zeroed in the right half (they cut the
   lower wing of the first line) and a 6 x 6 corner without data. Error of the model against the injected noise-free
   model over all channels in the zeroed half (RMS, K):

   | fit | RMS error | column error |
   |---|---|---|
   | cube without the zeroed block | 0.0275 | -0.2 % |
   | zeroed cube, zeros fitted as data (old behaviour) | 0.0600 | -3.2 % |
   | zeroed cube with ranges | 0.0296 | -0.4 % |

   At the corner without data the old fit has amplitudes pinned at 0; with ranges the model column there is 1.006 of
   the injected one, from the smoothing term alone.

Not checked here: real-data runs (the experiment repo records them); spatially varying noise between channel blocks
inside one pixel's range (needs a per-channel noise, not a range).
