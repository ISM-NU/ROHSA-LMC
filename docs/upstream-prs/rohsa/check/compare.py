"""Compare ROHSA outputs: bitwise, else per-pixel integrated-model differences."""
import sys, filecmp, numpy as np
def load(p):
    a = np.loadtxt(p, comments='#'); return a
def integ(a):
    i = a[:, 0].astype(int); j = a[:, 1].astype(int); I = a[:, 2] * a[:, 4] * np.sqrt(2 * np.pi)
    m = np.zeros((i.max() + 1, j.max() + 1)); np.add.at(m, (i, j), I); return m
for pair in sys.argv[1:]:
    x, y = pair.split(':'); px, py = f'runs/{x}.out', f'runs/{y}.out'
    if filecmp.cmp(px, py, shallow=False):
        print(f'{x:14s} vs {y:14s} BIT-FOR-BIT IDENTICAL'); continue
    A, B = integ(load(px)), integ(load(py)); r = np.abs(B - A) / np.maximum(np.abs(A), 1e-12)
    print(f'{x:14s} vs {y:14s} differ: per-pixel integrated model rel. diff median {np.median(r):.1e}, p99 {np.percentile(r, 99):.1e}, max {r.max():.1e}')
