"""Geometry of the heavy-ball method for general parameters (Part II of the paper).

Heavy ball  x_{k+1} = x_k - alpha grad f(x_k) + beta (x_k - x_{k-1})  on F_{mu,L}.
We use L = 1, so that the step parameter of the plots is aL = alpha*L and
x = alpha*mu = aL/kappa, y = alpha*L = aL.

Normalized recurrence (Lemma normalization of the paper):
    e_{k+1} = A e_k + B sigma_k - beta e_{k-1},
    A = 1 + beta - alpha (L+mu)/2,   B = alpha (L-mu)/2.
At the threshold t = 0 the angle function satisfies cot Theta(w) = Xi(cos w)/sin w (Lemma angle);
the quadratic Xi is called Q in this module,
    Q(c) = (B^2 - (1-beta)^2 - A^2 + 2A(1+beta)c - 4 beta c^2) / (2B(1-beta)).
"""
import numpy as np
from math import gcd
from scipy.optimize import linprog

# ----------------------------------------------------------------- parameters
def AB(kappa, aL, beta):
    x = aL / kappa
    A = 1 + beta - (aL + x) / 2
    B = (aL - x) / 2
    return A, B


def Qcoef(kappa, aL, beta):
    """Coefficients (q0, q1, q2) of Q; vectorized in aL and beta."""
    A, B = AB(kappa, aL, beta)
    den = 2 * B * (1 - beta)
    q0 = (B * B - (1 - beta) ** 2 - A * A) / den
    q1 = 2 * A * (1 + beta) / den
    q2 = -4 * beta / den
    return q0, q1, q2


def Qval(q, c):
    return q[0] + q[1] * c + q[2] * c * c


def theta(q, w):
    """Angle function Theta(w) in (0, pi)."""
    return np.arctan2(np.sin(w), Qval(q, np.cos(w)))


def E(x, y, b, c):
    """The lag-one quadratic L(c) of equation (E); Theta > w/2 on (0,pi] iff L > 0 on [-1,1) (Lemma cubic)."""
    return 4 * b * c ** 2 + 2 * (y - 1 - 2 * b - b ** 2 + b * x) * c + (2 - 2 * x + x * y + 2 * b ** 2 - 2 * b * y)


def min_on_interval_quad(a2, a1, a0, lo=-1.0, hi=1.0):
    """Minimum of a2 c^2 + a1 c + a0 over [lo, hi]; vectorized."""
    a2, a1, a0 = np.broadcast_arrays(np.asarray(a2, float), np.asarray(a1, float), np.asarray(a0, float))
    m = np.minimum(a2 * lo * lo + a1 * lo + a0, a2 * hi * hi + a1 * hi + a0)
    with np.errstate(divide='ignore', invalid='ignore'):
        cv = np.where(a2 > 0, -a1 / (2 * a2), np.nan)
    inside = (a2 > 0) & (cv > lo) & (cv < hi)
    mv = np.where(inside, a0 - a1 * a1 / (4 * np.where(a2 > 0, a2, 1)), np.inf)
    return np.minimum(m, mv)


def in_omega1(kappa, aL, beta):
    """Lag-one region Omega_1 (Corollary C(a)); explicit."""
    x, y, b = aL / kappa, aL, beta
    a2 = 4 * b
    a1 = 2 * (y - 1 - 2 * b - b ** 2 + b * x)
    a0 = 2 - 2 * x + x * y + 2 * b ** 2 - 2 * b * y
    return (min_on_interval_quad(a2, a1, a0) > 0) & (aL < 2 * (1 + beta))


def tangency_c(kappa, aL, beta):
    """Vertex c_* of E (tangency point on the boundary of Omega_1)."""
    x, y, b = aL / kappa, aL, beta
    return (1 + 2 * b + b * b - b * x - y) / (4 * b)


def root_moduli(kappa, aL, beta):
    x, y, b = aL / kappa, aL, beta
    out = []
    for p in (1 + b - x, 1 + b - y):
        disc = p * p - 4 * b
        r = np.where(disc <= 0, np.sqrt(b), (np.abs(p) + np.sqrt(np.maximum(disc, 0))) / 2)
        out.append(r)
    return out  # rho_a, rho_b


def underdamped(kappa, aL, beta):
    x, y, b = aL / kappa, aL, beta
    sb = np.sqrt(b)
    return (x >= (1 - sb) ** 2) & (y <= (1 + sb) ** 2)


def Kmu(r, b, y, c):
    return (4 * b * r ** 3 * c ** 2 + 2 * r * (-b ** 2 - b * r ** 3 + b * r ** 2 - 2 * b * r + r ** 3 * y - r ** 3) * c
            + b ** 2 * r ** 2 + b ** 2 - b * r ** 4 + b * r ** 3 - b * r ** 2 * y - b * r ** 2 - b * r * y + b * r
            + r ** 5 + r ** 4 * y - r ** 3 * y + r ** 3)


def KL_at_minus1(r, b, x):
    return (-3 * b ** 2 * r ** 2 - 2 * b ** 2 * r + b ** 2 + b * r ** 4 - 7 * b * r ** 3 + 3 * b * r ** 2 * x
            - 9 * b * r ** 2 + b * r * x - b * r + 3 * r ** 5 - r ** 4 * x + 2 * r ** 4 + r ** 3 * x - r ** 3)


def in_Qregion(kappa, aL, beta):
    """The explicit quadratic regime Q (Theorem B(b), Corollary C(c)); vectorized."""
    x, y, b = aL / kappa, aL, beta
    ra, rb = root_moduli(kappa, aL, beta)
    ud = underdamped(kappa, aL, beta)
    ok = np.zeros(np.broadcast(aL, beta).shape, bool)
    # case rho_a >= rho_b
    r = ra
    a2 = 4 * b * r ** 3
    a1 = 2 * r * (-b ** 2 - b * r ** 3 + b * r ** 2 - 2 * b * r + r ** 3 * y - r ** 3)
    a0 = Kmu(r, b, y, 0.0)
    caseI = (ra >= rb) & ~ud & (min_on_interval_quad(a2, a1, a0) >= 0)
    caseIII = (rb > ra) & ~ud & (KL_at_minus1(rb, b, x) >= 0)
    ok = caseI | caseIII
    return ok & in_omega1(kappa, aL, beta)


def coprime_pairs(kmax):
    return [(j, k) for k in range(2, kmax + 1) for j in range(1, (k + 1) // 2) if gcd(j, k) == 1 and 2 * j < k]


def tongue(kappa, aL, beta, j, k):
    """Single-harmonic tongue Theta(2 pi j/k) <= pi/k (Theorem A(ii), Corollary C(b))."""
    q = Qcoef(kappa, aL, beta)
    w = 2 * np.pi * j / k
    return (Qval(q, np.cos(w)) >= np.sin(w) / np.tan(np.pi / k)) & (aL < 2 * (1 + beta))


# ----------------------------------------------------------------- linear programs at t = 0
_W = np.linspace(1e-4, np.pi - 1e-4, 6000)
_W = np.unique(np.r_[_W, [2 * np.pi * j / k for k in range(2, 80) for j in range(1, k // 2 + 1) if gcd(j, k) == 1]])
_C, _S = np.cos(_W), np.sin(_W)


def lag_rows(q, lags, W=None):
    """Rows of the threshold problem at t = 0 divided by sin(Theta)/sin(w) (equation Lambdapoly).

    For a lag n>0 (lower) the row is Q(c)(1-cos nw) - sin w sin nw; for n<0 (upper)
    Q(c)(1-cos |n|w) + sin w sin |n|w; the minimizer row is Q(c)."""
    if W is None:
        W, C, S = _W, _C, _S
    else:
        C, S = np.cos(W), np.sin(W)
    Qc = Qval(q, C)
    cols = []
    for n in lags:
        m = abs(n)
        cols.append(Qc * (1 - np.cos(m * W)) - np.sign(n) * S * np.sin(m * W))
    cols.append(Qc)
    return np.array(cols).T


def certificate_margin(q, lags, W=None):
    """min over certificates (weights on lags + minimizer, sum 1) of max_w combination/(1-c+1e-6).
    Negative value: a strict certificate exists on the grid."""
    M = lag_rows(q, lags, W)
    C = _C if W is None else np.cos(W)
    nrm = (1 - C) + 1e-6
    nv = M.shape[1]
    r = linprog(np.r_[np.zeros(nv), 1], A_ub=np.hstack([M / nrm[:, None], -np.ones((len(nrm), 1))]),
                b_ub=np.zeros(len(nrm)), A_eq=[np.r_[np.ones(nv), 0]], b_eq=[1],
                bounds=[(0, None)] * nv + [(None, None)], method='highs')
    return r.x[-1] if r.status == 0 else np.nan


def periodic_feasible(q, period, harmonics=None):
    """Feasibility of the periodic program (Proposition periodic, Theorem fan-div(a)): weights on harmonics j of
    period `period`, all lags 1..period-1, minimizer row. Returns (feasible, margin, weights)."""
    if harmonics is None:
        harmonics = list(range(1, (period - 1) // 2 + 1))
    j = np.array(harmonics)
    w = 2 * np.pi * j / period
    th = theta(q, w)
    n = np.arange(1, period)[:, None]
    Lmb = np.cos(th)[None, :] - np.cos(n * w[None, :] - th[None, :])
    rows = np.vstack([Lmb, np.cos(th)[None, :]])
    p = len(j)
    r = linprog(np.r_[np.zeros(p), -1], A_ub=np.hstack([-rows, np.ones((rows.shape[0], 1))]),
                b_ub=np.zeros(rows.shape[0]), A_eq=[np.r_[np.ones(p), 0]], b_eq=[1],
                bounds=[(0, None)] * p + [(None, 1)], method='highs')
    if r.status != 0:
        return False, np.nan, None
    return r.x[-1] >= -1e-12, r.x[-1], r.x[:p]


def gap_K(q, npts=20000):
    """Return K if the dip {Theta <= w/2} lies in the rotation interval (1/(K+1), 1/K)."""
    w = np.linspace(1e-3, np.pi - 1e-3, npts)
    c = np.cos(w)
    d = Qval(q, c) - (1 + c)
    if not (d > 0).any():
        return None
    rot = w[d > 0] / (2 * np.pi)
    K = int(np.floor(1 / rot.max()))
    return K if (rot.min() > 1 / (K + 1) and rot.max() < 1 / K) else -1
