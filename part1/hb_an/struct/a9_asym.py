"""Expansions of Theorem 3 at kappa -> infinity and at kappa_c, exact.

(1) kappa -> infinity: rho_* is the positive root of F(rho, 1) = rho^4 - 2rho^3 + 2rho^2 - 2rho - 1,
    F_rho(rho_*, 1) != 0, and the first two coefficients of the expansion in kappa^{-1/2}.
(2) kappa -> kappa_c: V(t) = log phi + O(t^2) is even, sqrt(kappa_c) = phi^3, and
    rho_HB - 1 = (kappa - kappa_c)/(4 kappa_c) + O((kappa - kappa_c)^2) (equation critical)."""
import sympy as sp, pickle
t, u, e, k, z = sp.symbols('t u e k z')
C = sp.symbols('C')                      # C = c_* placeholder
mp_ = 4*C**4 - 8*C**2 + 4*C - 1
def red(expr):
    """reduce a rational function of C modulo the minimal polynomial (degree 4) -> cubic in C"""
    n, d = sp.fraction(sp.together(expr))
    n = sp.rem(sp.expand(n), mp_, C); d = sp.rem(sp.expand(d), mp_, C)
    # invert d modulo mp_
    s_, t_, g_ = sp.gcdex(d, mp_, C)
    assert sp.simplify(g_ - 1) == 0 or sp.degree(g_, C) == 0
    inv = sp.expand(s_/g_)
    return sp.rem(sp.expand(n*inv), mp_, C)
# clean derivation: implicit differentiation of F(rho,q)=0 at q=1 (F2), rho_* root of m(r)=r^4-2r^3+2r^2-2r-1
r, qv = sp.symbols('r qv')
Fp = r**2*(r**2 + 1) - 2*qv*r*(r**2 - r + 1) - qv**2*(r**2 + 1)
mr = r**4 - 2*r**3 + 2*r**2 - 2*r - 1
def redr(expr):
    n, d = sp.fraction(sp.together(expr))
    n = sp.rem(sp.expand(n), mr, r); d = sp.rem(sp.expand(d), mr, r)
    s_, t_, g_ = sp.gcdex(d, mr, r); inv = sp.expand(s_/g_)
    return sp.rem(sp.expand(n*inv), mr, r)
Fr, Fq = sp.diff(Fp, r), sp.diff(Fp, qv)
Frr, Frq, Fqq = sp.diff(Fp, r, 2), sp.diff(Fp, r, qv), sp.diff(Fp, qv, 2)
at1 = lambda ex: ex.subs(qv, 1)
print("F(r,1) - m(r) =", sp.expand(at1(Fp) - mr))
print("F_r(rho_*,1) != 0: resultant(F_r, m) =", sp.resultant(at1(Fr), mr, r))
rp = redr(-at1(Fq)/at1(Fr))
rpp = redr(-(at1(Fqq) + 2*at1(Frq)*rp + at1(Frr)*rp**2)/at1(Fr))
c1r = redr(2*rp); c2r = redr(2*rpp + 2*rp)          # q(e) = (1-e)/(1+e): q_e = -2, q_ee = 4
print("rho'(1) =", rp); print("c1 in Q(rho_*):", c1r); print("c2 in Q(rho_*):", c2r)
# conversion to c_* = ch t_* = (rho + 1/rho)/2: express via resultants and compare minimal polynomials / values
rho = red(C + (2*C - 1)/(2*C))                       # rho_* = ch + sh, sh t_* = (2c-1)/(2c)
print("rho_* in Q(c_*):", rho, ";  m(rho_*) mod =", red(rho**4 - 2*rho**3 + 2*rho**2 - 2*rho - 1))
c1 = red(c1r.subs(r, rho)); c2 = red(c2r.subs(r, rho))
print("c1 =", c1); print("c2 =", c2)
cstar = [r_ for r_ in sp.Poly(mp_, C).real_roots() if r_ > 1][0]
print("c_* =", sp.N(cstar, 25))
print("c1 ~", sp.N(c1.subs(C, cstar), 25), "   c2 ~", sp.N(c2.subs(C, cstar), 25))
# minimal polynomials
for name, val in (("c1", c1), ("c2", c2)):
    P = sp.resultant(sp.Symbol('Z') - val, mp_, C)
    print("minimal polynomial of", name, ":", sp.factor(P))
# cross-check against old formula (2 rho^3 + 2 rho + 2)/(4 rho^3 - 6 rho^2 + 4 rho - 2) * 2
old = 2*(2*rho**3 + 2*rho + 2)/(4*rho**3 - 6*rho**2 + 4*rho - 2)
print("c1 - 2*drho/dq|_{q=1} ([HB] formula) =", red(old - c1))
# numerical sanity check (not part of the proof): rho_kappa at kappa = 1e8
import mpmath as mpm
mpm.mp.dps = 40
for kap in (mpm.mpf(10)**6, mpm.mpf(10)**8):
    qq = (mpm.sqrt(kap) - 1)/(mpm.sqrt(kap) + 1)
    rr = mpm.findroot(lambda x: x**2*(x**2 + 1) - 2*qq*x*(x**2 - x + 1) - qq**2*(x**2 + 1), 1.7)
    ee = 1/mpm.sqrt(kap); rs = mpm.mpf(sp.N(rho.subs(C, cstar), 40))
    print("  kappa=%s: (rho - rho_* + c1 e)/e^2 =" % mpm.nstr(kap, 3), mpm.nstr((rr - rs + mpm.mpf(sp.N(c1.subs(C, cstar), 40))*ee)/ee**2, 12))
pickle.dump({'c1': c1, 'c2': c2}, open('a9.pkl', 'wb'))
# ---------------- (2) kappa -> kappa_c ----------------
print("\n(2) near kappa_c")
tt, dd = sp.symbols('tt dd')
Vt = sp.asinh(1 - 1/(2*sp.cosh(tt)))
vser = sp.series(Vt, tt, 0, 5).removeO()
print("V(t) =", vser, "  [even in t: V(-t)=V(t) exactly, since it depends on ch t]")
v2 = sp.nsimplify(vser.coeff(tt, 2))                 # = sqrt5/10
# u - u_c = t - (V(t) - V(0)) = t - v2 t^2 + O(t^4)  =>  t = w + v2 w^2 + O(w^3),  w = ln q - ln q_c
phi = (1 + sp.sqrt(5))/2
# ln q as function of kappa: d/dk ln q = k^{-1/2}(k-1)^{-1};  sqrt(kc) = phi^3, kc - 1 = 4 phi^3
A = 1/(phi**3*4*phi**3)
B = -sp.Rational(1, 2)/(phi**9*4*phi**3) - 1/(phi**3*(4*phi**3)**2)
print("check sqrt(kc) = phi^3:", sp.simplify(phi**6 - (9 + 4*sp.sqrt(5))), "; kc - 1 - 4phi^3:", sp.simplify(8 + 4*sp.sqrt(5) - 4*phi**3))
# w = A dd + B dd^2/2 + ... ; rho = e^t = 1 + w + (v2 + 1/2) w^2 + ...
a1 = sp.radsimp(sp.simplify(A)); a2 = sp.radsimp(sp.simplify(B/2 + (v2 + sp.Rational(1, 2))*A**2))
kc = 9 + 4*sp.sqrt(5)
print("rho - 1 = a1 (k-kc) + a2 (k-kc)^2 + O((k-kc)^3):  a1 =", a1, "=", sp.N(a1, 15), " = 1/(4kc)?", sp.simplify(a1 - 1/(4*kc)) == 0)
print("   a2 =", a2, "=", sp.N(a2, 15))
mpm.mp.dps = 30
for dk in (mpm.mpf('1e-3'), mpm.mpf('-1e-3')):
    kk = mpm.mpf(sp.N(kc, 30)) + dk; qq = (mpm.sqrt(kk) - 1)/(mpm.sqrt(kk) + 1)
    rr = mpm.findroot(lambda x: x**2*(x**2 + 1) - 2*qq*x*(x**2 - x + 1) - qq**2*(x**2 + 1), 1.0)
    print("  dk=%s: (rho-1-a1 dk)/dk^2 =" % mpm.nstr(dk, 2), mpm.nstr((rr - 1 - mpm.mpf(sp.N(a1, 30))*dk)/dk**2, 8))
