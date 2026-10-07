"""Explicit classification of the (alpha L, beta) plane by the closed-form criteria of Part II.

Labels: 0 alpha L >= 2(1+beta); 1 quadratic regime Q; 2 Omega_1 between the damping curves;
3 rest of Omega_1; 10+K single-harmonic tongue j=1 (smallest K); 5 tongues j>=2 only; 4 gap.
Output: data/phase_<kappa>.npz
"""
import sys, numpy as np
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'part2'))
import hbgeom as H

OUT = Path(__file__).resolve().parents[1] / 'data'
OUT.mkdir(exist_ok=True)


def classify(kappa, nx=1000, ny=700, aLmax=4.0):
    aL = np.linspace(aLmax / nx, aLmax, nx)
    be = np.linspace(0, 0.995, ny)
    Y, Bt = np.meshgrid(aL, be)
    lab = np.full(Y.shape, 4, np.int16)
    div = Y >= 2 * (1 + Bt)
    om = H.in_omega1(kappa, Y, Bt)
    Qr = H.in_Qregion(kappa, Y, Bt)
    ud = H.underdamped(kappa, Y, Bt)
    lab[om & ~Qr & ~ud] = 3
    lab[om & ud] = 2
    lab[om & Qr] = 1
    rest = ~om & ~div
    x = Y / kappa
    # j = 1 tongues: E(cos 2pi/K) <= 0 (Corollary C(b)); the label is the smallest such K
    bestK = np.zeros(Y.shape, np.int16)
    for K in range(3, 400):
        e = H.E(x, Y, Bt, np.cos(2 * np.pi / K))
        sel = rest & (e <= 0) & (bestK == 0)
        bestK[sel] = K
    lab[rest & (bestK > 0)] = 10 + np.minimum(bestK[rest & (bestK > 0)], 300)
    # tongues j >= 2
    q = H.Qcoef(kappa, Y, Bt)
    other = rest & (bestK == 0)
    for (j, k) in H.coprime_pairs(80):
        if j == 1:
            continue
        w = 2 * np.pi * j / k
        t = other & (H.Qval(q, np.cos(w)) >= np.sin(w) / np.tan(np.pi / k))
        lab[t] = 5
    lab[div] = 0
    return aL, be, lab


if __name__ == '__main__':
    for kappa in [float(a) for a in sys.argv[1:]] or [2, 10, 100, 1000]:
        aL, be, lab = classify(kappa)
        np.savez_compressed(OUT / f'phase_{kappa:g}.npz', aL=aL, beta=be, lab=lab)
        u, c = np.unique(lab, return_counts=True)
        print(kappa, dict(zip(u.tolist(), c.tolist())), flush=True)
