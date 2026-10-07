"""Lemma 56 (slacks of the initial pairs) and equation startupidentities, exact.

Checks, for the explicit data of equation explicitdata (z = rho e^{i omega} in the paper):
  (I1) g = 1 - |m|^2 = 2 rho (rho^2 - rho + 1)/(rho^2 + 1)^2 = 1/2 - (rho - 1)^4/(2(rho^2 + 1)^2) on F = 0,
       and m + Delta = (z + q^2)/(2q);
  (I2) R_{0k} -+ T_{0k} = a_-+ X^2 + Re(b_-+ Z) + c_-+ + (1 - (-q^2)^{ceil(k/2)})^2 (equation startupslack),
       for an arbitrary complex Z = z^k with X = |Z|, together with the form through E_{k-+1};
  (I3) the closed forms A(rho) + B(rho) q of all coefficients used in Table 5."""
import sympy as sp
from exactlib import *
N, Nb = sp.symbols('N Nb')
lam, mu = L, M
s = (lam + q**2/lam)/(2*q); sb = swap(s)
g = 2*r*(r**2 - r + 1)/(r**2 + 1)**2
c = (r + 1/r)/2
t_th = (r**2 - 1)/(r**2 + 1)           # th t
C = g*t_th
# (I1)
print("(I1) 1-|s|^2 == 2rho(rho^2-rho+1)/(rho^2+1)^2 mod F:", is_zero_modF(Re(1 - s*sb) - g))
print("(I1') g == 1/2-(rho-1)^4/(2(rho^2+1)^2):", sp.simplify(g - (sp.Rational(1, 2) - (r - 1)**4/(2*(r**2 + 1)**2))) == 0)
print("     also g == (rho^2-q^2)/(q(1+rho^2)) mod F:", is_zero_modF(g - (r**2 - q**2)/(q*(1 + r**2))))
# data
Res_ = Re(s)
cpsi = (C + 2*Res_)/4
Dl = (q/2)*(1 - 1/lam); Dlb = swap(Dl)
sig0 = s + Dl; sig0b = swap(sig0)
print("     sigma_0 == (lam+q^2)/(2q):", sp.simplify(sig0 - (lam + q**2)/(2*q)) == 0)
alp = Dlb*(s + 1)/2; alm = Dlb*(s - 1)/2
print("     alpha_+ == (1-1/mu)(lam+q)^2/(8 lam):", sp.simplify(alp - (1 - 1/mu)*(lam + q)**2/(8*lam)) == 0,
      "  alpha_- == (1-1/mu)(lam-q)^2/(8 lam):", sp.simplify(alm - (1 - 1/mu)*(lam - q)**2/(8*lam)) == 0)
def ReN(z): return (z + swap(z, extra=((N, Nb),)))/2
X2 = N*Nb
T = cpsi*X2 - cpsi - ReN(sp.expand((sig0b + sb*Nb)*(N - 1)))/2
R = ((N - 1)*(Nb - 1) - (s*N - sig0)*(sb*Nb - sig0b))/4
em = mu/r; ep = lam/r           # e^{-iw}, e^{iw}
am = g/(4*c*r); ap = g*r/(4*c)
bm = g*c*em + alp; bp = g*c*ep + alm
gm = g*r/(4*c) - ReN(alp) - Dl*Dlb/4
gp = g/(4*c*r) - ReN(alm) - Dl*Dlb/4
quad_m = am*X2 + ReN(bm*N) + gm
quad_p = ap*X2 + ReN(bp*N) + gp
# F3 form: (gX/2c) E_{k-+1}, X E_n = X ch(n t) + 2c^2 X cos(n w), X ch((k-+1)t) = (X^2 rho^{-+1} + rho^{+-1})/2,
# X cos((k-+1) w) = Re(N e^{-+iw})
F3m = (g/(2*c))*((X2/r + r)/2 + 2*c**2*ReN(N*em)) + ReN(alp*(N - 1)) - Dl*Dlb/4
F3p = (g/(2*c))*((X2*r + 1/r)/2 + 2*c**2*ReN(N*ep)) + ReN(alm*(N - 1)) - Dl*Dlb/4
def zero_identity(expr):
    """expr polynomial in N, Nb with coefficients rational in lam, mu, rho, q: check every coefficient vanishes on
    lam+mu=e1, lam mu = rho^2, F=0.  Coefficients of N^i Nb^j are not symmetric individually; we check
    c_{ij}(lam,mu) by reducing with mu = rho^2/lam and lam^2 - e1 lam + rho^2 = 0, then mod F."""
    e = sp.together(sp.expand(expr)); n, d = sp.fraction(e)
    P = sp.Poly(sp.expand(n), N, Nb)
    for mon, co in P.terms():
        co = sp.expand(co.subs(mu, e1 - lam))
        Pl = sp.Poly(co, lam); Rm = sp.rem(Pl, sp.Poly(lam**2 - e1*lam + e2, lam))
        for cc in [Rm.coeff_monomial(1), Rm.coeff_monomial(lam)]:
            if not is_zero_modF(cc): return False
    return True
print("(I2) R-T == a_-|N|^2 + Re(beta_- N) + gamma_-  :", zero_identity((R - T) - quad_m))
print("(I2) R+T == a_+|N|^2 + Re(beta_+ N) + gamma_+  :", zero_identity((R + T) - quad_p))
print("(I2) R-T == (gX/2c)E_{k-1} + Re(alpha_+(N-1)) - |Delta|^2/4 :", zero_identity((R - T) - F3m))
print("(I2) R+T == (gX/2c)E_{k+1} + Re(alpha_-(N-1)) - |Delta|^2/4 :", zero_identity((R + T) - F3p))
# wave pair (0,k) check of F3 itself: replace sigma_0 by s (Delta=0): R^w -+ T^w = (gX/2c) E_{k-+1}
Tw = cpsi*X2 - cpsi - ReN(sp.expand((sb + sb*Nb)*(N - 1)))/2
Rw = ((N - 1)*(Nb - 1) - (s*N - s)*(sb*Nb - sb))/4
print("(F3) R^w-T^w == (gX/2c)E_{k-1}:", zero_identity((Rw - Tw) - (g/(2*c))*((X2/r + r)/2 + 2*c**2*ReN(N*em))))
print("(F3) R^w+T^w == (gX/2c)E_{k+1}:", zero_identity((Rw + Tw) - (g/(2*c))*((X2*r + 1/r)/2 + 2*c**2*ReN(N*ep))))
# (I3) closed forms
print("(I3) closed forms  (value = A(rho) + B(rho) q on F=0):")
items = [('a_-', am), ('a_+', ap), ('gamma_-', symred(gm)), ('gamma_+', symred(gp)),
         ('|beta_-|^2', Re(bm*swap(bm))), ('|beta_+|^2', Re(bp*swap(bp))),
         ('B_0^- = Re beta_-', Re(bm)), ('B_1^- = Re(beta_- e^{iw})', Re(bm*ep)), ('B_2^- = Re(beta_- e^{2iw})', Re(bm*ep**2)),
         ('B_0^+ = Re beta_+', Re(bp)), ('B_1^+ = Re(beta_+ e^{iw})', Re(bp*ep)), ('B_2^+ = Re(beta_+ e^{2iw})', Re(bp*ep**2)),
         ('|Delta|^2', Re(Dl*Dlb)), ('Re Delta', Re(Dl)), ('|sigma_0|^2', Re(sig0*sig0b)), ('C', C),
         ('T_{0*} = -c_psi + Re(sigma_0)/2', -cpsi + Re(sig0)/2)]
for nm, e in items:
    A, B = modF(e); print("   %-28s A = %s\n   %-28s B = %s" % (nm, A, '', B))
print("     T_{0*} == Re(Delta)/2 - C/4 :", is_zero_modF(-cpsi + Re(sig0)/2 - (Re(Dl)/2 - C/4)))
print("     Re Delta == (q/2)(rho^2+2)/(rho^2+1):", is_zero_modF(Re(Dl) - (q/2)*(r**2 + 2)/(r**2 + 1)))
# (I4) kappa = kappa_c exactly: rho = 1, q = qc = (sqrt5-1)/2, lam = e^{2 pi i/3}; spiral parts of R -+ T for k mod 3
print("(I4) kappa_c: exact values in Q(sqrt5)")
qc = (sp.sqrt(5) - 1)/2
l3 = sp.Rational(-1, 2) + sp.sqrt(3)*sp.I/2
sc = sp.expand((l3 + qc**2/l3)/(2*qc)); Dc = sp.expand((qc/2)*(1 - 1/l3)); s0c = sc + Dc; cpc = sp.re(sc)/2
nz = lambda e: sp.nsimplify(sp.radsimp(sp.simplify(sp.expand(e))))
print("   s =", nz(sc), ", |s|^2 =", nz(sc*sp.conjugate(sc)), ", |Delta|^2 =", nz(Dc*sp.conjugate(Dc)), ", |sigma_0|^2 =", nz(s0c*sp.conjugate(s0c)))
for k in (1, 2, 3):
    Nk = sp.expand(l3**k)
    Tk = -sp.re(sp.expand(sp.conjugate(s0c + sc*Nk)*(Nk - 1)))/2
    Rk = (sp.expand((Nk - 1)*sp.conjugate(Nk - 1)) - sp.expand((sc*Nk - s0c)*sp.conjugate(sc*Nk - s0c)))/4
    print("   k = %d mod 3:  R-T (spiral) = %s,  R+T (spiral) = %s" % (k % 3, nz(Rk - Tk), nz(Rk + Tk)))
print("   T_{0*} =", nz(sp.re(s0c)/2 - cpc), ";  (1-qc^4)^2 =", nz((1 - qc**4)**2), "; (289-126 sqrt5)>0:", 289**2 > 126**2*5)
