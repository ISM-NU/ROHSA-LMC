# ROHSA Full-Cube Laptop Design, 2026-06-26

Status: active implementation design, not validated science output.

## Goal

Run ROHSA on the full LMC H I footprint at full spatial resolution on a laptop
class machine without freezing the operating system.

The scientific motivation is valid: ROHSA's spatial regularization is designed
for a whole PPV cube, not disconnected crops. The engineering problem is that
the upstream implementation assumes the whole numeric optimization state fits
in RAM.

## Current Constraint

For the default HNV cube,
`data/exp_raw/Tb_combined_mw_full_hnv_resOG_reprojected_on_gass.fits`, the image
shape is `3624 x 3624 x 106`. ROHSA pads the spatial grid to powers of two for
the multiresolution descent, so this becomes a `4096 x 4096 x 106` working grid.

The corrected preflight estimate for the production `f_g_cube_fast` path at
`N = 6`, `m = 10` is about `128 GiB` peak. An earlier estimate of about 400
GiB incorrectly included the four-dimensional derivative tensor allocated by
the older, unused `f_g_cube` routine. This remains far above laptop memory but
fits a normal Irgetas ZEN4 node.
This is not caused by FITS storage format. The large allocations happen after
the data are decompressed and loaded as numeric arrays.

Relevant upstream allocation points:

- `rohsa/upstream/ROHSA/src/mod_inout.f90` reads the full `.dat` input into
  `cube(nv, ny, nx)`.
- `rohsa/upstream/ROHSA/src/mod_functions.f90` flattens the full parameter cube
  into `beta`, `lb`, and `ub` arrays for L-BFGS-B.
- `rohsa/upstream/ROHSA/src/mod_optimize.f90` production routine
  `f_g_cube_fast` allocates full residual and 3N gradient-map work arrays. The
  larger `dF_over_dB(3*n_gauss, dim_v, dim_y, dim_x)` belongs to the retained
  but unused `f_g_cube` routine.

The exact upstream full-resolution path still does not fit a 14 GiB laptop.
The production allocation is dominated by padded data/residual cubes and the
L-BFGS-B vectors and workspace, which scale with the full padded sky grid.

## Does This Require Changing Fortran?

Yes, for a credible full-resolution laptop run.

Wrapper-only changes can reduce disk use, avoid Julia-side full-cube copies, and
stop unsafe runs before they crash the machine. They cannot remove the upstream
Fortran optimizer's full-grid arrays.

There are two viable implementation choices:

1. Patch the upstream Fortran source in place.
2. Add a separate project-owned out-of-core optimizer that reproduces ROHSA's
   objective closely enough for validation, while retaining the upstream binary
   as the reference implementation on small crops.

The first path is closer to upstream but harder to make maintainable. The second
path is easier to test in Julia, but it must be treated as a new implementation
until crop parity is demonstrated.

## Why `.fz` Compression Is Not Enough

FITS tile compression, such as `.fits.fz`, can reduce disk footprint and may
make region reads cheaper. It does not reduce the optimizer's live memory once
pixels are read into numeric arrays.

Useful roles for compression:

- Store source cubes and intermediate caches more compactly.
- Read row, plane, or tile chunks into the wrapper without materializing the
  whole FITS cube.
- Reduce disk pressure for generated model and residual products.

Not useful for:

- Avoiding `cube(nv, ny, nx)` inside current ROHSA.
- Avoiding full parameter, bound, gradient, or derivative arrays inside current
  L-BFGS-B.
- Making the exact upstream full-resolution run fit in 14 GiB RAM.

## Proposed Architecture

The laptop implementation should be a staged out-of-core ROHSA-compatible path.
It should keep the upstream Fortran binary as the reference path and introduce a
new experimental runner only after bounded parity tests pass.

### Stage 1: Memory-Safe Input Cache

Create a disk-backed cache for the full cube:

- Read FITS or `.fits.fz` in slabs using FITSIO region reads.
- Sanitize nonfinite and nonpositive values consistently with the current ROHSA
  wrapper.
- Write a binary `Float32` or quantized cache with explicit metadata:
  `nx`, `ny`, `nv`, velocity axis, WCS, noise estimate, stride, and source path.
- Prefer memory-mapped arrays or chunked JLD2/HDF5/Zarr-style storage over the
  current text `.dat` format for the experimental runner.

This stage is wrapper-side and does not require Fortran changes.

### Stage 2: Streaming Objective Kernels

Replace the largest ROHSA work arrays with kernels that accumulate objective
and gradient terms tile by tile.

Required change:

- Do not allocate `dF_over_dB(3*n_gauss, dim_v, dim_y, dim_x)`.
- Do not allocate a full residual cube.
- Compute each spectrum's model, residual, and parameter derivatives locally,
  then accumulate into the local parameter-gradient tile.
- Compute spatial regularization terms from parameter maps using halo reads.

This stage requires changing Fortran or writing a new optimizer implementation.

### Stage 3: Block Coordinate Optimization

The current upstream final update flattens every parameter into one global
L-BFGS-B vector. That vector and its work arrays are too large for the laptop.

Instead:

- Store parameter maps on disk.
- Optimize spatial tiles with a halo wide enough for the regularization stencil.
- Hold only one tile's spectra, parameters, bounds, gradients, and L-BFGS state
  in memory.
- Sweep over all tiles repeatedly until global objective improvement stalls.
- Keep global width anchors `b_params` as shared state updated after each sweep
  or after a full pass over component maps.

This is not bit-for-bit upstream ROHSA. It is an out-of-core approximation to
the same regularized objective. It must be validated against upstream on crops.

### Stage 4: Product Derivation Without Full Component Cubes

After fitting, derive products from parameter maps in slabs:

- model cube
- residual cube
- component cubes when explicitly requested
- moment0, NHI, centroid, sigma, Doppler-temperature maps
- grouped broad/intermediate/narrow products

Do not reconstruct all `N` component cubes simultaneously. The existing product
derivation can be kept for small crops, but full-cube products need a streaming
path.

## Validation Plan

Validation must happen before using the laptop full-cube result scientifically.

1. Run upstream ROHSA on small crops already known to fit in memory.
2. Run the out-of-core path on the same crops with identical hyperparameters.
3. Compare:
   - total model/data flux
   - voxel RMSE and MAE
   - moment0 residual RMSE and bias
   - residual skewness
   - component amplitude, centroid, and sigma maps
   - grouped NHI maps
4. Set acceptance thresholds before full-cube execution.
5. Only then run the full-resolution laptop path.
6. Treat the first full-cube output as experimental until VGT/starlight
   validation improves relative to non-ROHSA baselines.

## Implementation Tasks

1. Add a binary cube-cache writer and reader.
2. Add a streaming product writer for ROHSA parameter maps.
3. Prototype objective/gradient evaluation without full derivative tensors.
4. Run crop parity against upstream ROHSA for one tile-sized crop.
5. Add block-coordinate sweeps with halo regularization.
6. Add objective trace and convergence diagnostics per sweep.
7. Run multi-crop parity tests.
8. Run one full-footprint laptop job with conservative checkpointing.

## Risks

- Block-coordinate optimization may converge to a different local solution than
  global L-BFGS-B.
- Halo handling can create tile boundary artifacts.
- `Float32` caches reduce memory and disk pressure but may shift fitted
  parameters; this needs crop parity evidence.
- Runtime may still be long, even if memory is controlled.
- The first successful full-cube run is not automatically a better VGT input.

## Decision

Use `.fits.fz` or other compression only as a storage and streaming aid.

For true full-resolution laptop execution, implement an experimental out-of-core
ROHSA-compatible optimizer and validate it against upstream ROHSA on crops. Do
not claim scientific equivalence until parity and downstream validation are
documented.
