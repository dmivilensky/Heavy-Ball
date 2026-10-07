"""Lattice and fan statements of Sections levels and level1 (Appendix levels).

Exact integer checks:
  Theorem farey(a): for every reduced j/k < 1/2 with k <= 80, the lattice lines n xi - l eta = 1
    through (j, k) with |n| < k are exactly n = b and n = -d (Farey parents a/b < j/k < c/d), and
    they pass through (c, d) and (a, b).
  Theorem farey(b): for one atom at 2 pi j/k with phase Theta, the violated lags n in 1..k-1 are those
    with j n mod k in {1, ..., ceil(k Theta/pi) - 1} (compared with the sign of Lambda_n, many Theta).
  Lemma fan: Farey neighbors, parents, the line U and the fan lines, and the congruences, K, p <= 40.
Exact symbolic check:
  equation Lambdapoly: Lambda_n = (sin Theta / sin w){Xi (1 - T_|n|(c)) - sgn(n)(1 - c^2) U_{|n|-1}(c)}.
Floating-point checks (random data, consistency of the closed forms):
  Theorem lattice(b): sign rule of Lambda_n in lattice coordinates.
  Theorem fan-div(b): the four values Lambda_{1}, Lambda_{-K} at w_p, w_{p-1}; the existence of a
    two-atom measure for the lags 1, -K is equivalent to Delta_p >= 0 (bisection on the ratio).
  Theorem fan-div(c): Delta_2 / sin(pi/k) equals the displayed mediant expression.
  Lemma pin: D(w_j) = E(w_j) = 0 and the derivative formulas, by finite differences on Theta of a tuning.
"""
import math
import random
from fractions import Fraction
from math import gcd
import numpy as np
import sympy as sp

FAILED = []


def report(name, ok):
    print(('OK    ' if ok else 'FAIL  ') + name)
    if not ok:
        FAILED.append(name)


def farey_parents(j, k):
    # a/b < j/k < c/d with b j - a k = 1 and c k - d j = 1, 0 < b, d < k
    b = pow(j, -1, k)
    a = (b * j - 1) // k
    d = k - b
    c = j - a
    return a, b, c, d


def check_farey():
    ok_a, ok_b = True, True
    for k in range(3, 81):
        for j in range(1, (k + 1) // 2):
            if gcd(j, k) != 1 or 2 * j >= k:
                continue
            a, b, c, d = farey_parents(j, k)
            ok_a &= (b * j - a * k == 1) and (c * k - d * j == 1) and (b * c - a * d == 1) and (a + c == j) and (b + d == k)
            sols = sorted(n for n in range(-k + 1, k) if n != 0 and (n * j - 1) % k == 0)
            ok_a &= sols == sorted([b, -d])
            ok_a &= (b * c - a * d == 1) and (-d * a + c * b == 1)
            for th in np.linspace(0.01, math.pi - 0.01, 37):
                w = 2 * math.pi * j / k
                viol = {n for n in range(1, k) if math.cos(th) - math.cos(n * w - th) < -1e-12}
                R = math.ceil(k * th / math.pi) - 1
                pred = {n for n in range(1, k) if 1 <= (j * n) % k <= R}
                # exclude exact boundary cases k*th/pi integer
                if abs(k * th / math.pi - round(k * th / math.pi)) > 1e-9:
                    ok_b &= viol == pred
    report('Theorem farey(a): lines through (j,k) with |n|<k are the parent lags b and -d (k <= 80)', ok_a)
    report('Theorem farey(b): violated lags of one atom at 2 pi j/k (k <= 80, 37 phases each)', ok_b)


def check_fan():
    ok = True
    for K in range(2, 41):
        for p in range(1, 41):
            k = p * K + 1
            kp = (p - 1) * K + 1
            ok &= Fraction(p, k) > Fraction(p - 1, kp) and Fraction(p, k) < Fraction(1, K)
            ok &= kp * p - (p - 1) * k == 1 and k - p * K == 1 and kp - (p - 1) * K == 1
            if p >= 2:
                ok &= farey_parents(p, k) == (p - 1, kp, 1, K)
            ok &= -K * p + k == 1
            ok &= ((-p * K) - 1) % k == 0 and (((p - 1) * K + 1) + K) % k == 0
        for m in range(0, 41):
            ok &= (m * K + 1) * 1 - m * K == 1
            ok &= (m * K + 1) * (m + 1) - m * ((m + 1) * K + 1) == 1
    report('Lemma fan: neighbors, parents, U, fan lines and congruences (K, p, m <= 40)', ok)


def check_lambdapoly():
    th, w = sp.symbols('theta omega')
    c = sp.cos(w)
    ok = True
    for n in list(range(1, 9)) + list(range(-8, 0)):
        lam = sp.cos(th) - sp.cos(n * w - th)
        m = abs(n)
        Xi = sp.cos(th) / sp.sin(th) * sp.sin(w)
        rhs = sp.sin(th) / sp.sin(w) * (Xi * (1 - sp.chebyshevt(m, c)) - sp.sign(n) * (1 - c ** 2) * sp.chebyshevu(m - 1, c))
        ok &= sp.simplify(sp.expand_trig(lam - rhs)) == 0
    report('equation Lambdapoly for |n| <= 8 (exact)', ok)


def check_lattice_sign():
    random.seed(2)
    ok = True
    for _ in range(4000):
        xi = random.uniform(0, 3); eta = random.uniform(1.05, 20); n = random.choice([i for i in range(-25, 26) if i])
        lam = 2 * math.sin(math.pi * n * xi / eta) * math.sin(math.pi * (n * xi - 1) / eta)
        neg = any(0 < n * xi - l * eta < 1 for l in range(-200, 201))
        if abs(lam) > 1e-9:
            ok &= (lam < 0) == neg
    report('Theorem lattice(b): Lambda_n < 0 iff 0 < n xi - l eta < 1 (4000 random points)', ok)


def lam(n, w, th):
    return math.cos(th) - math.cos(n * w - th)


def check_fan_div():
    random.seed(3)
    ok_vals, ok_delta, ok_med = True, True, True
    for _ in range(3000):
        K = random.randint(2, 12); p = random.randint(2, K + 1); k = p * K + 1
        wp, wq = 2 * math.pi * p / k, 2 * math.pi * (p - 1) / k
        thp = random.uniform(math.pi / k, math.pi * p / k)
        thq = random.uniform(math.pi * (p - 1) / k, math.pi * (K + 1) / k)
        vals = [lam(1, wp, thp), lam(-K, wp, thp), lam(1, wq, thq), lam(-K, wq, thq)]
        pred = [2 * math.sin(math.pi * p / k) * math.sin(math.pi * p / k - thp),
                2 * math.sin(math.pi / k) * math.sin(math.pi / k - thp),
                2 * math.sin(math.pi * (p - 1) / k) * math.sin(math.pi * (p - 1) / k - thq),
                2 * math.sin(math.pi * (K + 1) / k) * math.sin(math.pi * (K + 1) / k - thq)]
        ok_vals &= max(abs(u - v) for u, v in zip(vals, pred)) < 1e-12
        ok_vals &= vals[0] > 0 > vals[1] and vals[3] > 0 > vals[2]
        D = (math.sin(math.pi * p / k) * math.sin(math.pi * p / k - thp) * math.sin(math.pi * (K + 1) / k)
             * math.sin(math.pi * (K + 1) / k - thq)
             - math.sin(math.pi * (p - 1) / k) * math.sin(thq - math.pi * (p - 1) / k) * math.sin(math.pi / k)
             * math.sin(thp - math.pi / k))
        lo = -vals[1] / vals[3]
        hi = vals[0] / -vals[2]
        feasible = lo <= hi
        if abs(D) > 1e-10:
            ok_delta &= feasible == (D >= 0)
        if p == 2:
            v1, v2 = thq - math.pi / k, thp - math.pi / k
            med = (2 * math.cos(math.pi / k) * math.cos(math.pi / (2 * k)) * math.sin(math.pi / k - v2)
                   * math.cos(v1 + math.pi / (2 * k)) - math.sin(math.pi / k) * math.sin(v1) * math.sin(v2))
            ok_med &= abs(D / math.sin(math.pi / k) - med) < 1e-12
    report('Theorem fan-div(b): values of Lambda_1, Lambda_{-K} at the fan frequencies and their signs', ok_vals)
    report('Theorem fan-div(b): a two-atom measure satisfies lags 1 and -K iff Delta_p >= 0', ok_delta)
    report('Theorem fan-div(c): Delta_2 / sin(pi/k) equals the mediant expression', ok_med)


def check_pin():
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import hbgeom as H
    q = H.Qcoef(100.0, 70471 / 31800, 6169 / 9900)
    Th = lambda w: float(H.theta(q, w))
    ok = True
    for K, p in ((4, 3), (4, 2), (3, 2), (5, 3)):
        k = p * K + 1
        for j in range(1, (k + 1) // 2):
            wj = 2 * math.pi * j / k
            D = lambda w: lam(-p * K, w, Th(w)) - lam(1, w, Th(w))
            E = lambda w: lam((p - 1) * K + 1, w, Th(w)) - lam(-K, w, Th(w))
            e = 1e-6
            dD = (D(wj + e) - D(wj - e)) / (2 * e)
            dE = (E(wj + e) - E(wj - e)) / (2 * e)
            th = Th(wj)
            ok &= abs(D(wj)) < 1e-12 and abs(E(wj)) < 1e-12
            ok &= abs(dD + k * math.sin(wj - th)) < 1e-5 * k and abs(dE + k * math.sin(K * wj + th)) < 1e-5 * k
    report('Lemma pin: D, E vanish at 2 pi j/k with the stated derivatives (finite differences)', ok)


def main():
    check_farey()
    check_fan()
    check_lambdapoly()
    check_lattice_sign()
    check_fan_div()
    check_pin()
    if FAILED:
        raise SystemExit('FAILED: ' + ', '.join(FAILED))
    print('ALL OK')


if __name__ == '__main__':
    main()
