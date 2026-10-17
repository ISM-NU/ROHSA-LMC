# Optional binary input format

Branch `binary-input` · base `master` · independent.

## Problem

ROHSA reads its input cube from a text file with one line per voxel. For a
2048 × 2048 × 159 cube that is 667 million lines and 20 GB: writing it takes
several minutes (7 min from Julia on our cluster) and parsing it in
`read_cube` about 4 minutes on one core, before the fit starts.

## Change

- `read_cube`: a file name ending in `.bin` is read as a binary stream: three
  4-byte integers `nv, ny, nx`, then the `nv·ny·nx` 8-byte values in Fortran
  order `cube(v,y,x)`. Any other name is read as text, exactly as before.
- ROHSApy: `cube2bin(filename)` writes that format from the same `(v, y, x)`
  NumPy cube that `cube2dat` takes.
- README: a short "Binary input" section.

The text format stays the default; nothing changes for existing scripts.

## Checks

On the synthetic 256 × 256 × 100 cube of `check/make_synthetic.py` (N = 4):

- `.dat` input: decomposition bit-for-bit identical to `master`.
- `.bin` file holding the values parsed from the `.dat` file: bit-for-bit
  identical decomposition to the `.dat` run.
- `.bin` written directly by `cube2bin`: the values differ from the
  16-decimal text by at most 5.6e-17, and the decomposition by rounding
  (per-pixel integrated model: median 5.8e-5 relative). This is the size of
  difference that a last-bit change of the input causes.
- `read_cube`: 2.37 s for the 191 MB text file, 0.03 s for the 52 MB binary
  file.

