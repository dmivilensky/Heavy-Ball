"""Regularity of the boundary arcs (Theorem 3, regularity), exact.

For a curve P(G, u) = 0 and an arc u = S(G), the rate is analytic in log q where P_u != 0 and
d/dt[t - V(t)] != 0. The script computes Res_u(P, P_u) and Res_u(P, N) with
N = P_u^2 (1 + u^2) - P_G^2 (G^2 - 1) and counts their real roots in the G-ranges of the spiral,
P2 and P1 (Sturm sequences)."""
import sympy as sp, pickle, sys
sys.path.insert(0, '.')
from cs_form import cs, ss
d = pickle.load(open('P12.pkl', 'rb')); P1, P2 = d['P1'], d['P2']
Psp = 2*cs*ss - 2*cs + 1          # spiral
def check(P, name, lo, hi):
    P = sp.factor_list(P)[1]; P = [f for f, k in P if f.has(ss)][0]
    Ps, Pc = sp.diff(P, ss), sp.diff(P, cs)
    r1 = sp.Poly(sp.resultant(P, Ps, ss), cs)
    Nn = sp.expand(Ps**2*(1 + ss**2) - Pc**2*(cs**2 - 1))
    r2 = sp.Poly(sp.resultant(P, Nn, ss), cs)
    n1 = r1.count_roots(lo, hi) if r1.degree() > 0 else 0
    n2 = r2.count_roots(lo, hi) if r2.degree() > 0 else 0
    print("%s: #real roots in [%s, %s] of Res(P,P_s): %d ; of Res(P,N): %d" % (name, lo, hi, n1, n2))
    assert n1 == n2 == 0, (name, n1, n2)
    return r1, r2
check(Psp, "spiral", 1, sp.oo)
# two-atom arc: |t| in [t_a, t_1]: use the safe c-range [ch(0.2), c_1], c_1 root of 32c^3-36c^2+4c-1
c1 = sp.real_roots(sp.Poly(32*cs**3 - 36*cs**2 + 4*cs - 1, cs))[0]
lo2 = sp.Rational(10201, 10000)          # = 1.0201 < ch(0.2) = 1.02007...? check below
print("ch(0.2) =", sp.N(sp.cosh(sp.Rational(1, 5))), " ch t_a (t_a~0.2149) ~", sp.N(sp.cosh(sp.Rational(2149, 10000))), " c_1 =", sp.N(c1))
check(P2, "two-atom P2", sp.Rational(102, 100), sp.Rational(104, 100))
# one-wave arc: |t| >= t_W, ch t_W = (1+sqrt3)/2 ~ 1.366
check(P1, "one-wave P1", sp.Rational(136, 100), sp.oo)
