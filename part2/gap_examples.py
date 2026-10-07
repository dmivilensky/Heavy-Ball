"""Propositions "A two-harmonic thirteen-cycle" (prop:13) and "Convergence in a gap" (prop:conv-gap):
both behaviors inside the gap I_4 at kappa = 100.

Proposition prop:conv-gap is proved in exact rational arithmetic:
  * the dip {L(c) < 0} of the lag-one quadratic is a nonempty interval (c_-, c_+)
    with 0 < c_- < c_+ < cos(2 pi/5) = (sqrt5 - 1)/4, i.e. rotation numbers in (1/5, 1/4);
  * P(c) = (1-c^2)(1 - r U_3(c)) - Xi(c)[(1-c) + r(1 - T_4(c)) + gamma] > 0 on [-1, 1]
    for r = 13/100, gamma = 1/1000, decided by a Sturm count of the real roots of P
    on [-1, 1] and the sign of P at one point.
Proposition prop:13 is proved in interval arithmetic (mpmath.iv, 200 bits):
  * the same gap membership, in exact arithmetic;
  * the spectrum nu = (149/797) delta_{2.2pi/13} + (648/797) delta_{3.2pi/13} at t = 0 has all raw
    lag slacks >= 0 (lags 1..12, both signs; lags divisible by 13 vanish) and gbar > 0 with C = 0,
    so that w = |g - 2is| nu satisfies the hypothesis of Theorem fan-div(a);
  * Theta > pi/6 on (0, pi]: Xi(c) <= 0 or Xi(c) < cot(pi/6) sqrt(1 - c^2) on 4000 boxes covering [-1, 1].
"""
from fractions import Fraction as Fr
import sympy as sp
from mpmath import iv, mpf

iv.prec = 200
c = sp.symbols('c')


def data(kappa, aL, beta):
    x = aL / kappa; y = aL
    A = 1 + beta - (aL + x) / 2; B = (aL - x) / 2
    den = 2 * B * (1 - beta)
    Xi = ((B * B - (1 - beta) ** 2 - A * A) + 2 * A * (1 + beta) * c - 4 * beta * c ** 2) / den
    Lc = 4 * beta * c ** 2 + 2 * (y - 1 - 2 * beta - beta ** 2 + beta * x) * c + (2 - 2 * x + x * y + 2 * beta ** 2 - 2 * beta * y)
    return A, B, sp.nsimplify(sp.expand(Xi)), sp.expand(Lc), x, y


def in_gap_I4(Lc, y, x, beta):
    """Exact check that the dip of L lies strictly inside (0, cos 2pi/5)."""
    roots = sp.Poly(Lc, c).nroots(n=50)
    rr = [sp.re(z) for z in roots if abs(sp.im(z)) < 1e-40]
    assert len(rr) == 2
    # exact: L(0) > 0, L(cos 2pi/5) > 0, L at the vertex < 0, vertex in (0, cos 2pi/5)
    # cos(2pi/5) = (sqrt5 - 1)/4 is the positive root of 4c^2 + 2c - 1 and lies in (lo5, hi5)
    lo5, hi5 = sp.Rational(309, 1000), sp.Rational(3091, 10000)
    assert 4 * lo5 ** 2 + 2 * lo5 - 1 < 0 < 4 * hi5 ** 2 + 2 * hi5 - 1
    vert = sp.Rational(-sp.Poly(Lc, c).all_coeffs()[1], 2 * sp.Poly(Lc, c).all_coeffs()[0])
    ok = (Lc.subs(c, 0) > 0) and (Lc.subs(c, lo5) > 0) and (sp.Poly(Lc, c).count_roots(lo5, hi5) == 0) \
        and (Lc.subs(c, vert) < 0) and (0 < vert) and (vert < lo5)
    # (y - x)(1 - beta)(1 + c - Xi) = L, so the dip is {Theta < omega/2}
    return bool(ok), [float(v) for v in rr]


def prop_conv():
    kappa, aL, beta = Fr(100), Fr(111, 50), Fr(31, 50)
    A, B, Xi, Lc, x, y = data(sp.Rational(100), sp.Rational(111, 50), sp.Rational(31, 50))
    assert sp.expand(Lc - (y - x) * (1 - sp.Rational(31, 50)) * (1 + c - Xi)) == 0
    ok, rr = in_gap_I4(Lc, y, x, beta)
    print('Proposition prop:conv-gap: tuning kappa=100, alpha L=111/50, beta=31/50')
    print('  Xi(c) =', Xi)
    print('  dip of L in (0, cos 2pi/5):', ok, ' roots of L:', rr)
    r, g = sp.Rational(13, 100), sp.Rational(1, 1000)
    T4 = sp.chebyshevt(4, c); U3 = sp.chebyshevu(3, c)
    P = sp.expand((1 - c ** 2) * (1 - r * U3) - Xi * ((1 - c) + r * (1 - T4) + g))
    nroots = sp.Poly(P, c).count_roots(-1, 1)
    val = P.subs(c, 0)
    print('  P(c) =', sp.Poly(P, c).as_expr())
    print('  real roots of P on [-1,1] (Sturm):', nroots, '; P(0) =', val)
    assert ok and nroots == 0 and val > 0
    print('  two-lag inequality holds with r = 13/100, gamma = 1/1000: True')


def prop_13():
    kappa, aL, beta = Fr(100), Fr(70471, 31800), Fr(6169, 9900)
    A, B, Xi, Lc, x, y = data(sp.Rational(100), sp.Rational(70471, 31800), sp.Rational(6169, 9900))
    ok, rr = in_gap_I4(Lc, y, x, beta)
    print('Proposition prop:13: tuning kappa=100, alpha L=70471/31800, beta=6169/9900')
    print('  dip of L in (0, cos 2pi/5):', ok, ' roots of L:', rr)
    assert ok
    Ar = Fr(str(A)); Br = Fr(str(B))
    I = lambda q: iv.mpf(q.numerator) / q.denominator
    Ai, Bi, bi = I(Ar), I(Br), I(beta)
    pi = iv.pi
    L = 13
    weights = {2: Fr(149, 797), 3: Fr(648, 797)}

    def gs(j):
        w = 2 * pi * j / L
        cw = iv.cos(w); sw = iv.sin(w)
        re = ((1 + bi) * cw - Ai) / Bi
        im = (1 - bi) * sw / Bi
        return 1 - (re ** 2 + im ** 2), im
    mins = None
    for n in range(1, L):
        for sg in (1, -1):
            val = iv.mpf(0)
            for j, p in weights.items():
                g_, s_ = gs(j); w = 2 * pi * j / L
                val += I(p) * (g_ * (1 - iv.cos(n * w)) - sg * 2 * s_ * iv.sin(n * w))
            if mins is None or val.a < mins.a:
                mins = val
    gbar = sum((I(p) * gs(j)[0] for j, p in weights.items()), iv.mpf(0))
    print('  min lag slack, lags 1..12 both signs, lies in', mins)
    print('  gbar lies in', gbar)
    assert mins.a > 0 and gbar.a > 0
    # the numbers stated in Proposition 9.4
    assert mins.a >= iv.mpf('3.8594e-4').b and iv.mpf('0.79908').b < gbar.a and gbar.b < iv.mpf('0.79909').a
    q = sp.Poly(Xi, c).all_coeffs()[::-1]
    q0, q1, q2 = [I(Fr(int(sp.fraction(v)[0]), int(sp.fraction(v)[1]))) for v in q]
    theta0 = pi / 6; ct = iv.cos(theta0) / iv.sin(theta0); N = 4000; okth = True
    for i in range(N):
        C = iv.mpf([mpf(-1) + mpf(2 * i) / N, mpf(-1) + mpf(2 * (i + 1)) / N])
        Q = q0 + q1 * C + q2 * C * C
        lo = 1 - C * C
        rhs = ct * iv.sqrt(iv.mpf([max(mpf(0), lo.a), lo.b]))
        if not (Q.b <= 0 or Q.b < rhs.a):
            okth = False
            break
    print('  Theta > pi/6 on (0, pi]:', okth)
    assert okth


if __name__ == '__main__':
    prop_conv()
    prop_13()
    print('ALL OK')
