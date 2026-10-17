"""Synthetic test cube for the ROHSA pull requests.

Three spatially coherent Gaussian components plus noise, written with
ROHSApy's cube2dat (text) and cube2bin (binary), a noise map, and a
parameter file. Usage: python make_synthetic.py <ROHSApy dir> <out dir>
"""
import os, sys
import numpy as np
from scipy.ndimage import gaussian_filter
os.environ.setdefault("MPLBACKEND", "Agg")
sys.path.insert(0, sys.argv[1])
from ROHSApy import ROHSA

out = sys.argv[2]; os.makedirs(out, exist_ok=True)
rng = np.random.default_rng(2026)
nv, ny, nx = 100, 256, 256

def field(mean, spread, scale=12):
    f = gaussian_filter(rng.standard_normal((ny, nx)), scale)
    return mean + spread * f / f.std()

v = np.arange(1, nv + 1, dtype=float)[:, None, None]
cube = np.zeros((nv, ny, nx))
for a, mu, sig in [((2.0, 0.5), (40, 3), (6, 1)), ((1.5, 0.4), (55, 4), (2.5, 0.4)), ((1.0, 0.3), (70, 3), (10, 2))]:
    A, M, Sg = field(*a), field(*mu), np.maximum(field(*sig), 1.0)
    cube += A * np.exp(-(v - M) ** 2 / (2 * Sg ** 2))
cube += 0.05 * rng.standard_normal(cube.shape)

core = ROHSA(cube)
core.cube2dat(filename=os.path.join(out, "synthetic.dat"))
core.cube2bin(filename=os.path.join(out, "synthetic.bin"))
core.rms_map(np.full((ny, nx), 0.05), filename=os.path.join(out, "synthetic_rms.dat"))
