"""The first level unrolled along the boundary of Omega_1.

For fixed kappa and each beta, alpha_d(beta) is the boundary point of Omega_1 (smallest alpha L
outside Omega_1). Points are sampled as alpha L = alpha_d(beta) + s, 0 <= s <= smax, and classified:
  K-tongue (any K)                              -> 100 + K (the figures relabel by the smallest K)
  gap, two-lag inequality (Theorem fan-conv)    -> 1
  gap, fan program feasible (Theorem fan-div), p -> 200 + p  (period pK+1)
  gap, other single-harmonic tongue j/k         -> 3
  gap, certificate with fan lags (LP)           -> 2
  gap, none of the above                        -> 4
  alpha L >= 2(1+beta)                          -> 0
Output: data/gap_<kappa>.npz with beta, alpha_d, rotation number of the tangency, s, labels.
"""
import sys, numpy as np
from pathlib import Path
from scipy.optimize import brentq
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'part2'))
import hbgeom as H

TAG = ''
OUT = Path(__file__).resolve().parents[1] / 'data'

CG = np.cos(np.linspace(1e-3, np.pi - 1e-3, 4000))


def minE(kappa, aL, beta):
    x, y, b = aL / kappa, aL, beta
    return H.min_on_interval_quad(4 * b, 2 * (y - 1 - 2 * b - b * b + b * x), 2 - 2 * x + x * y + 2 * b * b - 2 * b * y)


def boundary(kappa, beta):
    top = 2 * (1 + beta) - 1e-9
    if minE(kappa, top, beta) > 0:
        return None
    lo, hi = 1e-6, top
    if minE(kappa, lo, beta) <= 0:
        return None
    return brentq(lambda a: float(minE(kappa, a, beta)), lo, hi, xtol=1e-13)


def two_lag(q, K):
    """Two-lag inequality: exists r>=0 with Xi[(1-c)+r(1-T_K)] < (1-c^2)(1-r U_{K-1}) on (-1,1)."""
    c = CG
    w = np.arccos(c)
    Xi = H.Qval(q, c)
    TK = np.cos(K * w)
    UK = np.sin(K * w) / np.sin(w)
    a = (Xi * (1 - c) - (1 - c * c)) / (1 - c)
    b = (Xi * (1 - TK) + (1 - c * c) * UK) / (1 - c)
    neg = b < 0
    pos = b > 0
    zero = ~neg & ~pos
    if zero.any() and (a[zero] >= 0).any():
        return False
    lo = max(0.0, np.max(a[neg] / -b[neg])) if neg.any() else 0.0
    hi = np.min(-a[pos] / b[pos]) if pos.any() else np.inf
    return lo < hi


def classify_point(kappa, aL, beta, K, pmax=30):
    q = H.Qcoef(kappa, aL, beta)
    if two_lag(q, K):
        return 1
    for p in range(2, pmax + 1):
        ok, m, wts = H.periodic_feasible(q, p * K + 1)
        if ok:
            return 200 + p
    for (j, k) in H.coprime_pairs(60):
        w = 2 * np.pi * j / k
        if H.Qval(q, np.cos(w)) >= np.sin(w) / np.tan(np.pi / k):
            return 3
    fan = [1 + m * K for m in range(0, 60 // K + 1)] + [-m * K for m in range(1, 60 // K + 1)]
    if H.certificate_margin(q, fan) < -1e-9:
        return 2
    return 4


def run(kappa, nb=420, ns=140, smax=None, bmin=None, bmax=0.97):
    betas = np.linspace(bmin, bmax, nb)
    rows = []
    for beta in betas:
        ad = boundary(kappa, beta)
        rows.append(np.nan if ad is None else ad)
    ad = np.array(rows)
    good = np.isfinite(ad)
    betas, ad = betas[good], ad[good]
    cstar = H.tangency_c(kappa, ad, betas)
    rot = np.arccos(np.clip(cstar, -1, 1)) / (2 * np.pi)
    s = np.linspace(0, smax, ns + 1)[1:]
    lab = np.zeros((len(s), len(betas)), np.int16)
    for i, (beta, a0) in enumerate(zip(betas, ad)):
        for jj, ds in enumerate(s):
            aL = a0 + ds
            if aL >= 2 * (1 + beta):
                lab[jj, i] = 0
                continue
            x = aL / kappa
            best, bestK = np.inf, 0
            for K in range(3, 200):
                e = H.E(x, aL, beta, np.cos(2 * np.pi / K))
                if e <= 0 and e < best:
                    best, bestK = e, K
            if bestK:
                lab[jj, i] = 100 + bestK
                continue
            q = H.Qcoef(kappa, aL, beta)
            K = H.gap_K(q)
            if K is None or K < 2:
                lab[jj, i] = 4
                continue
            lab[jj, i] = classify_point(kappa, aL, beta, K)
        if i % 20 == 0:
            print(kappa, i, len(betas), flush=True)
    np.savez_compressed(OUT / f'gap_{kappa:g}{TAG}.npz', beta=betas, ad=ad, rot=rot, s=s, lab=lab)


if __name__ == '__main__':
    kappa = float(sys.argv[1])
    bmin = float(sys.argv[2]); bmax = float(sys.argv[3]); smax = float(sys.argv[4])
    nb = int(sys.argv[5]); ns = int(sys.argv[6])
    TAG = sys.argv[7] if len(sys.argv) > 7 else ''
    run(kappa, nb=nb, ns=ns, smax=smax, bmin=bmin, bmax=bmax)
