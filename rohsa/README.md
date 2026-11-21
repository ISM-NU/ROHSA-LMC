# ROHSA Preprocessing

This folder contains the project-side ROHSA preprocessing workflow.

For the NU Irgetas supercomputer, start with the beginner runbook at
`hpc/irgetas/README.md`. It adds Slurm and memory guardrails without changing
the ROHSA estimator.

- `upstream/ROHSA/` is the nested upstream ROHSA source tree. Keep it separate
  from project adapter code.
- `../configs/rohsa/` holds the editable run configurations (default `full_cube.toml`).
- `configs.jl` loads the TOML config and serializes the Fortran parameter file.
- `inputs.jl` writes full-cube or explicit-crop FITS, ROHSA `.dat`, and
  RMS-map inputs.
- `runner.jl` writes the Fortran parameter file, finds or builds the nested
  upstream binary, and invokes it.
- `products_stage.jl` consumes ROHSA `.dat` output and derives FITS products.
- `selection.jl` selects a decomposition from reporting metrics when multiple
  `N` values are tested.
- `visuals.jl` writes optional PNG panels and the HTML visual index.
- `report.jl` writes CSV tables, `manifest.toml`, and the primary markdown
  report.
- `cli/preprocess.jl` is the user-facing entrypoint.

The actual Gaussian decomposition is performed by the upstream Fortran
executable at `upstream/ROHSA/src/ROHSA`. Project Julia code prepares inputs,
records trace files, and derives products from ROHSA's fitted parameters; it
does not reimplement the ROHSA optimizer.

Run the main preprocessing wrapper from the repository root:

```bash
julia --project=. rohsa/cli/preprocess.jl
```

The default config is `configs/rohsa/full_cube.toml`. It describes the exact
full-resolution input intent, but the wrapper now runs a memory preflight before
writing large `.dat` files or invoking ROHSA. On laptops, the full-resolution
LMC cube is expected to exceed memory by a large margin and will be rejected
unless `force_full_cube = true` is set in the config.

For a laptop-safe full-footprint probe, use the explicit reduced-resolution
config:

```bash
julia --project=. rohsa/cli/preprocess.jl --config configs/rohsa/full_cube_laptop_stride8.toml
```

`full_cube_laptop_stride8.toml` keeps the full sky footprint but sets
`spatial_stride = 8`, disables immediate product reconstruction, and disables
visual generation. This is not equivalent to the paper's full-resolution
whole-cube fit; it is a controlled way to exercise ROHSA on the same footprint
within a small RAM budget.

To inspect memory feasibility without writing inputs or running ROHSA, set
`preflight_only = true` in a config and run the same command. The preflight
reports the padded ROHSA grid, estimated peak memory, and a suggested
`spatial_stride` when the requested run is unsafe.

If the upstream binary is already built, it is used directly. If it is absent,
the wrapper builds `upstream/ROHSA/src` with `make`. The nested upstream tree is
the only automatic ROHSA source.

The wrapper derives ROHSA products only. It does not run VGT and does not
compare ROHSA products against default VGT outputs.

Edit `configs/rohsa/full_cube.toml` to choose the input cube, output directory, crop,
visual-output setting, FITS output extension, and every ROHSA parameter. The
checked-in values mirror `upstream/ROHSA/ROHSApy.ipynb`; there is no
project-specific hidden preset.

Set `fits_extension = ".fits.fz"` to write ROHSA FITS products through the same
path-based compression helper used by VGT outputs. This requires `fpack` on
`PATH`; leave the default `.fits` for uncompressed products.

Set `spectral_kernel` in `[input]` (`"hanning"` or `"gauss<σ in channels>"`, default
`"none"`) to smooth every spectrum before ROHSA. Smoothing runs on the raw crop,
before negatives are clipped, and the noise passed to ROHSA is the raw estimate
times sqrt(Σk²) (the channel-difference estimator cannot measure correlated
noise). Fitted widths include the kernel in quadrature; the kernel's equivalent
σ is written to the crop header (`SPKSIGCH`, channels) and the run summary, and
analyses must deconvolve with it before applying width cuts.

Missing data and noise, for mosaics whose pointings cover different channel
ranges (all opt-in in `[input]`; the defaults reproduce earlier runs):

- `missing_data = "skip"`: the wrapper writes each pixel's live channel range
  (first to last channel that is neither 0 nor NaN in the raw crop) to
  `live_channel_range.dat`, and ROHSA leaves the channels outside it out of its
  data term and gradient (`filename_range`, mod_array.f90). Pixels without any
  data get noise 0 and an empty range: no data term, parameters from the
  smoothing term only (do not read them as measurements). The coarse
  multiresolution levels average each channel over the children observed in
  it and take the union of their ranges (mod_functions.f90 `mean_array_live`).
  The default `"zero"` writes NaN as 0 and fits those zeros as data.
- `channel_range_file` (with `"skip"`): reliable ranges made beforehand by
  `rohsa/cli/channel_ranges.jl`. In a mosaic, end-of-band channels that only one
  of two overlapping blocks covers can be averaged with the empty block, so
  signal and noise fade towards the seam while staying nonzero. The tool
  measures that weight w from the noise (adjacent-channel differences per
  16 px tile, relative to the most common noise level) and keeps, per pixel,
  the outermost live channels with w >= 0.8 (at most a 20 % amplitude error in
  kept voxels; no 1/w rescaling). Use the same file for a null cube made from
  the data. The line-free noise estimate then uses only these channels.
- `noise_estimate = "line_free"` with `noise_per_pixel = true`: per-pixel noise
  from emission-free voxels (exact zeros and NaN are never used; with
  `"skip"`, only channels in the pixel's range). `noise_floor_quantile` floors
  the map at that quantile over the pixels with data. Attenuated channels bias
  the estimate low, so without reliable ranges a high floor (e.g. 0.25) is
  needed; with them the map has no low tail and a guard (e.g. 0.01) keeps real
  noise variations. A per-run quantile gives a data cube and its null different
  floors, so keep it low when they are compared.

Evidence: `docs/validation/rohsa-live-channel-range-2026-11-22.md`.

To use a different config file:

```bash
julia --project=. rohsa/cli/preprocess.jl --config path/to/rohsa.toml
```

When several `n_gauss` values are supplied in the config, the wrapper runs
every candidate, keeps every candidate product and trace, and auto-selects one
result. The full candidate list is written to `rohsa_fit_metrics.csv` and
`rohsa_product_paths.csv`; the chosen products are written separately to
`rohsa_selection_summary.csv`, `selected_rohsa_products.csv`, and
`selected_rohsa_config.toml`.

## Auto-Config

ROHSA auto-config is optional and stays inside preprocessing. Enable it with
`[auto_config] enabled = true` in `configs/rohsa/full_cube.toml`. The normal `[rohsa]`
table remains the base config; any auto-config axis that is omitted uses the
single value from `[rohsa]`.

Resolution algorithm:

1. Build a deterministic candidate grid from listed axes such as `n_gauss`,
   `lambda_amp`, `lambda_mu`, `lambda_sig`, and `lambda_var_sig`.
2. If the grid exceeds `max_candidates`, keep a deterministic coarse spread of
   candidates instead of random sampling.
3. Run every retained candidate with an isolated candidate ID, so repeated
   `n_gauss` values with different regularization do not overwrite outputs.
4. Score successful candidates as
   `reconstruction_score + physical_penalty + stability_penalty`.
5. Select the lowest score per region and write the complete ranking to
   `rohsa_auto_config_candidates.csv`.

`reconstruction_score` is the existing flux/moment residual score.
`physical_penalty` discourages negligible components, sigma-bound pileup,
excess broad-component fraction, and duplicate components. `stability_penalty`
discourages rough amplitude, centroid, and linewidth maps. This selector is a
preprocessing diagnostic, not downstream VGT or starlight validation.
