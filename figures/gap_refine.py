"""Refinement of the first-level classification (data/gap_*.npz, written by gap_data.py).

Grid points that gap_data.py leaves with label 4 are tested again with longer fan programs
(p <= 79) and with fan-lag certificates using lags up to 200 on 12000 frequencies. Points that pass
the certificate test receive label 2; points with a feasible fan program receive 200 + p.
Output: data/gap_<kappa><tag>_refined.npz with the same fields.
"""
import sys
import numpy as np
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'part2'))
import hbgeom as H

DATA = Path(__file__).resolve().parents[1] / 'data'
W = np.linspace(1e-4, np.pi - 1e-4, 12000)


def refine(name, kappa=100.0):
    d = dict(np.load(DATA / f'{name}.npz'))
    lab = d['lab'].copy()
    for j, i in np.argwhere(lab == 4):
        beta = d['beta'][i]
        aL = d['ad'][i] + d['s'][j]
        q = H.Qcoef(kappa, aL, beta)
        K = H.gap_K(q)
        new = 4
        for p in range(2, 80):
            if H.periodic_feasible(q, p * K + 1)[0]:
                new = 200 + p
                break
        if new == 4:
            lags = [1 + m * K for m in range(0, 200 // K + 1)] + [-m * K for m in range(1, 200 // K + 1)]
            if H.certificate_margin(q, lags, W) < -1e-10:
                new = 2
        lab[j, i] = new
    d['lab'] = lab
    np.savez_compressed(DATA / f'{name}_refined.npz', **d)
    u, c = np.unique(lab, return_counts=True)
    print(name, dict(zip(u.tolist(), c.tolist())), flush=True)


if __name__ == '__main__':
    for name in sys.argv[1:] or ['gap_100', 'gap_100_hi']:
        refine(name)
