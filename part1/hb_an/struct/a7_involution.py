"""Lemma 36(i) and the normalized functionals of Lemma 37, exact.

Checks: the map (t, q, C) -> (-t, q e^{-2t}, -C) at fixed (v, omega) preserves a_m and the C-terms;
the spiral polynomial satisfies F(1/rho, q/rho^2) rho^6 = F(rho, q); the normalized functionals
L_m, U_m are even in t at fixed v; their values at c = cosh t and at c = +-1 (the values
-sinh^2 v r_m (cosh mt -+ (+-1)^m) used in the proof of Lemma 37)."""
import sympy as sp
from common import *
t, v, w, C, m = sp.symbols('t v w C m', real=True)
# (1) single-wave constraint quantities as functions of (t, v, w, C) for symbolic m
gg = sp.sin(w)**2 - sp.sinh(v)**2
am = gg*(sp.cosh(m*t) - sp.cos(m*w)); bm = 2*sp.sinh(v)*sp.sin(w)*sp.sin(m*w)
xm = C*sp.sinh(m*t) - bm
iota = {t: -t, C: -C}
print("(1) a_m(iota) - a_m      =", sp.simplify(am.subs(iota, simultaneous=True) - am))
print("    x_m(iota) - x_m      =", sp.simplify(xm.subs(iota, simultaneous=True) - xm))
print("    star |C|<=g: g unchanged, |C| unchanged (trivial)")
# HB orbit: s = (lam + q^2/lam)/(2q) = ch(v + i w), v = t - ln q: check with lam = e^{t+iw}
qq = sp.symbols('qq', positive=True)
lam = sp.exp(t + sp.I*w)
s_orbit = (lam + qq**2/lam)/(2*qq)
print("(1') s_orbit - ch(v+iw) at v=t-ln q:", sp.simplify((s_orbit - sp.cosh(t - sp.log(qq) + sp.I*w)).rewrite(sp.exp)))
# (2) F invariance: F(1/rho, q/rho^2) = rho^{-6} F(rho, q)
rho = sp.symbols('rho', positive=True)
print("(2) F(1/rho,q/rho^2)*rho^6 - F(rho,q) =", sp.expand(F(1/rho, q/rho**2)*rho**6 - F(rho, q)))
# H vs F: H_q(rho) and F(rho,q) define the same branch: check  H = F ? (both quartics)
print("    H_q(rho) - F(rho,q) =", sp.expand(H(rho, q) - F(rho, q)))
# (3) F in (X,Y) coordinates: F(X, X/Y) * Y^2/X^2 = S(X,Y)
print("(3) F(X,X/Y)*Y^2/X^2 - S =", sp.expand(F(X, X/Y)*Y**2/X**2 - S))
# and S = 2Y(1+X^2) (sh v - 1 + X/(1+X^2)),  X/(1+X^2) = 1/(2 ch t)
print("    S - 2Y(1+X^2)(shv - 1 + 1/(2 cht)) =", sp.simplify(S - 2*Y*(1 + X**2)*(shv - 1 + 1/(2*chX(1)))))
# (4) evenness of the C-free normalized functionals in t at fixed v: L_m(1/X) = L_m(X), U_m likewise
for mm in range(1, 8):
    for fn, nm in ((Lfun, 'L'), (Ufun, 'U')):
        e = fn(mm)
        assert sp.simplify(e.subs(X, 1/X) - e) == 0, (nm, mm)
print("(4) L_m, U_m even in t (X->1/X) at fixed Y for m=1..7: OK  [general m: a_m, r_m even, b_m t-free]")
# (5) the 3-term certificate symbol Kt = L1 + W U2 + (1-W) U3: numerator R, and R(c = ch t) = 0
Kt = sp.together(Lfun(1) + W*Ufun(2) + (1 - W)*Ufun(3))
num, den = sp.fraction(sp.factor(Kt))
R = sp.expand(num)
print("(5) denominator of Kt:", sp.factor(den))
print("    R(X,Y,c=ch t,W) =", sp.simplify(R.subs(c, chX(1))))
print("    R(1/X)*X^k - R : k with identity:", [k for k in range(0, 20) if sp.expand(R.subs(X, 1/X)*X**k - R) == 0])
# general fact: every L_m, U_m vanishes at c = ch t up to the sign: L_m(ch t) = +2 sh v sh^2 t, U_m(ch t) = -2 sh v sh^2 t
for mm in range(1, 7):
    assert sp.simplify(Lfun(mm).subs(c, chX(1)) - 2*shv*shX(1)**2) == 0
    assert sp.simplify(Ufun(mm).subs(c, chX(1)) + 2*shv*shX(1)**2) == 0
print("    L_m(c=ch t) = 2 sh v sh^2 t = -U_m(c=ch t), m=1..6: OK")
# (6) endpoint values: L_m, U_m at c = +-1 equal -sh^2 v r_m (ch mt -+ ...) < 0
for mm in range(1, 7):
    for cc in (1, -1):
        val = sp.simplify(Lfun(mm).subs(c, cc)); val2 = sp.simplify(Ufun(mm).subs(c, cc))
        target = -shv**2*r(mm)*(chX(mm) - cc**mm)
        assert sp.simplify(val - target) == 0 and sp.simplify(val2 - target) == 0
print("(6) L_m(+-1) = U_m(+-1) = -sh^2 v * r_m (ch mt - (+-1)^m), m=1..6: OK")
