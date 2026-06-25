# ROHSA Preprocessing

This folder contains the project-side ROHSA preprocessing workflow.

- `upstream/ROHSA/` is the nested upstream ROHSA source tree. Keep it separate
  from project adapter code.
- `config.toml` is the editable run configuration.
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

The default config is `rohsa/config.toml`. It describes the exact
full-resolution input intent, but the wrapper now runs a memory preflight before
writing large `.dat` files or invoking ROHSA. On laptops, the full-resolution
LMC cube is expected to exceed memory by a large margin and will be rejected
unless `force_full_cube = true` is set in the config.

For a laptop-safe full-footprint probe, use the explicit reduced-resolution
config:

```bash
julia --project=. rohsa/cli/preprocess.jl --config rohsa/config.laptop.toml
```

`config.laptop.toml` keeps the full sky footprint but sets
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

Edit `rohsa/config.toml` to choose the input cube, output directory, crop,
visual-output setting, and every ROHSA parameter. The checked-in values mirror
`upstream/ROHSA/ROHSApy.ipynb`; there is no project-specific hidden preset.

To use a different config file:

```bash
julia --project=. rohsa/cli/preprocess.jl --config path/to/rohsa.toml
```

When several `n_gauss` values are supplied in `config.toml`, the wrapper runs
every candidate, keeps every candidate product and trace, and auto-selects one
result. The full candidate list is written to `rohsa_fit_metrics.csv` and
`rohsa_product_paths.csv`; the chosen products are written separately to
`rohsa_selection_summary.csv`, `selected_rohsa_products.csv`, and
`selected_rohsa_config.toml`.
