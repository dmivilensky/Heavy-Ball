"""Theorem convex: the formula for Theta'' and the identity 2 Xi M^2 - N(1 - c^2 + Xi^2) = W/2,
where W is the polynomial written in sans-serif type in the paper.

With Xi(c) = q0 + q1 c + q2 c^2 and h = cot Theta = Xi(cos w)/sin w:
  (1) h' = -M/sin^2 w and h'' = N/sin^3 w with the polynomials M, N of the proof (exact);
  (2) 2 Xi M^2 - N(1 - c^2 + Xi^2) = W/2 as a polynomial identity in (q0, q1, q2, c) (exact);
  (3) Theta'' = W/(2 sin^5 w (1 + h^2)^2), compared with direct differentiation of arccot h at
      random points (floating point, a consistency check of the assembled formula).
"""
import random
import sympy as sp

w = sp.symbols('omega')
q0, q1, q2, c = sp.symbols('q0 q1 q2 c')
Xi = lambda cc: q0 + q1 * cc + q2 * cc ** 2
A1, A2 = -Xi(1), -Xi(-1)
X, Y, Q = A1 * (1 + c) ** 2, A2 * (1 - c) ** 2, Xi(c)
W = (Q * (X - Y) ** 2 + Q ** 2 * (X * (3 * c + 2) + Y * (2 - 3 * c)) + 2 * Q ** 3 * (1 + c ** 2) + 2 * Q * (1 - c ** 2) ** 2
     + A1 * (1 + c) ** 3 * (1 - c) * (2 - c) + A2 * (1 - c) ** 3 * (1 + c) * (2 + c))
M = (q1 + 2 * q2 * c) * (1 - c ** 2) + Q * c
N = q2 * c ** 4 + (q0 - q2) * c ** 2 + 2 * q1 * c + q0 + 2 * q2


def main():
    h = Xi(sp.cos(w)) / sp.sin(w)
    s = sp.sin(w)
    ok1 = sp.simplify(sp.diff(h, w) + M.subs(c, sp.cos(w)) / s ** 2) == 0
    ok2 = sp.simplify(sp.diff(h, w, 2) - N.subs(c, sp.cos(w)) / s ** 3) == 0
    print('OK    ' if ok1 and ok2 else 'FAIL  ', "h' = -M/sin^2 w and h'' = N/sin^3 w")
    ok3 = sp.expand(2 * Q * M ** 2 - N * (1 - c ** 2 + Q ** 2) - W / 2) == 0
    print('OK    ' if ok3 else 'FAIL  ', '2 Xi M^2 - N(1 - c^2 + Xi^2) = W/2 (polynomial identity)')
    Tpp = sp.diff(sp.acot(h), w, 2)
    random.seed(0)
    worst = 0.0
    for _ in range(25):
        v = {q0: random.uniform(-1, 1), q1: random.uniform(-1, 1), q2: random.uniform(-2, 0), w: random.uniform(0.1, 3.0)}
        lhs = float(Tpp.subs(v))
        hv = float(h.subs(v))
        rhs = float(W.subs({q0: v[q0], q1: v[q1], q2: v[q2], c: sp.cos(v[w])})) / (2 * float(sp.sin(v[w])) ** 5 * (1 + hv ** 2) ** 2)
        worst = max(worst, abs(lhs - rhs) / (1 + abs(lhs)))
    ok4 = worst < 1e-9
    print('OK    ' if ok4 else 'FAIL  ', "Theta'' = W/(2 sin^5 w (1 + cot^2 Theta)^2) at 25 random points (rel. diff %.1e)" % worst)
    if not (ok1 and ok2 and ok3 and ok4):
        raise SystemExit('FAILED')
    print('ALL OK')


if __name__ == '__main__':
    main()
