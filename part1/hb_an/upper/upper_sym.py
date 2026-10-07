"""Lemma 11 (dual factorization) and the spiral certificate (equation spiralcertificate), exact.

Notation: X = e^t, Y = e^v, G = cosh t, u = sinh v, c = cos(omega), k = 1 - c^2, g = k - u^2,
h = 2uk, so that b_m = h U_{m-1}(c); the symbol with lower lag 1 and upper lags J is
K = (a_1 - b_1) + sum_m W_m r_m (a_m + b_m), r_m = sinh t / sinh(mt), sum_m W_m = 1.
Checks:
  (S1) K = (G - c)(g alpha_J - h beta_J) for J = {2,3,5} with general weights; r_m = 1/U_{m-1}(G);
       for J = {2,3} the polynomials M, N of Lemma 46(b);
  (S2) the identity behind equation spiralpoly on the spiral branch, and g = u/G;
  (S3) the weight W(G) = 2(1 + 2G - 2G^2)/(4G + 1) that makes the zero at c = -1/(2G) double,
       and the derivative of Proposition 41 (rank two of the active Jacobian);
  (S4) the factorization of equation spiralcertificate, the residual quadratic
       2Gc^2 - (2G - 1)c + 4G(G - 1), its vertex value and discriminant 32G^3 - 36G^2 + 4G - 1;
  (S5) the cubic 32G^3 - 36G^2 + 4G - 1 has one real root G_1 (Sturm), and G_W = (1 + sqrt 3)/2;
  (S6) minimal polynomials of rho, q and kappa at kappa_W, kappa_1, kappa_2 (Table 1)."""
import sympy as sp
X, Y, c, G, u, q = sp.symbols('X Y c G u q', positive=True)
W2, W3, W5, W = sp.symbols('W2 W3 W5 W')
def T(m, x): return sp.chebyshevt(m, x)
def U(m, x): return sp.chebyshevu(m, x)
def ok(e):
    assert sp.simplify(sp.together(e)) == 0, e
    return True
print("=== (S1) factorisation K = (G-c)(g*alpha - h*beta), general weights, lags {2,3,5}")
Wd = {2: W2, 3: W3, 5: 1 - W2 - W3}
chm = lambda m: (X**m + X**-m)/2; shm = lambda m: (X**m - X**-m)/2
shv = (Y - 1/Y)/2; k = 1 - c**2; gX = k - shv**2
aX = lambda m: gX*(chm(m) - T(m, c)); bX = lambda m: 2*shv*k*U(m - 1, c)
K_orig = aX(1) - bX(1) + sum(Wd[m]*shm(1)/shm(m)*(aX(m) + bX(m)) for m in Wd)
def alpha_beta(Wd, G, c):
    al = 1 + sum(Wd[m]*sp.cancel((T(m, G) - T(m, c))/(G - c))/U(m - 1, G) for m in Wd)
    be = sum(Wd[m]*sp.cancel((U(m - 1, G) - U(m - 1, c))/(G - c))/U(m - 1, G) for m in Wd)
    return al, be
al, be = alpha_beta(Wd, G, c)
gG = k - u**2; hG = 2*u*k
K_fact = (G - c)*(gG*al - hG*be)
sub = {G: (X + 1/X)/2, u: shv}
print("  identity holds:", ok(K_orig - K_fact.subs(sub)))
print("  r_m = 1/U_{m-1}(ch t):", all(ok(shm(1)/shm(m) - 1/U(m - 1, (X + 1/X)/2)) for m in (2, 3, 5)))
# lags {2,3}: alpha = M/(G(4G^2-1)), beta = N/(G(4G^2-1))
al23, be23 = alpha_beta({2: W, 3: 1 - W}, G, c)
N = W*(4*G**2 - 1) + 4*(1 - W)*G*(G + c)
M = 8*G**3 + 4*G**2*c - G + (1 - W)*(4*G*c**2 - 3*G) - W*(G + c)
print("  lags {2,3}: alpha = M/(G(4G^2-1)), beta = N/(G(4G^2-1)):", ok(al23 - M/(G*(4*G**2 - 1))), ok(be23 - N/(G*(4*G**2 - 1))))

print("=== (S2) master equation <=> u = 1 - 1/(2G)")
H = X**4 - 2*q*X**3 + (1 + 2*q - q**2)*X**2 - 2*q*X - q**2
lhs = ((X/q - q/X)/2 - (1 - 1/(X + 1/X)))
print("  (Y-1/Y)/2 - 1 + 1/(2G) = H_q(X)/(2qX(X^2+1)) :", ok(lhs - H/(2*q*X*(X**2 + 1))))
ub = 1 - 1/(2*G); ck = -1/(2*G)
print("  on the branch g = u/G at c = c_k:", ok((1 - ck**2 - ub**2) - ub/G))

print("=== (S3) the tight weight W(G) and the rank condition D2 != D3")
K23 = sp.together((G - c)*((1 - c**2 - u**2)*al23 - 2*u*(1 - c**2)*be23))
Kb = sp.together(K23.subs(u, ub))
print("  K(t_k, c_k) = 0 :", ok(Kb.subs(c, ck)))
dK = sp.together(sp.diff(Kb, c).subs(c, ck))      # d/dw = -sin w d/dc
Wsol = sp.solve(sp.numer(dK), W); assert len(Wsol) == 1
WG = 2*(1 + 2*G - 2*G**2)/(4*G + 1)
print("  W(G) = 2(1+2G-2G^2)/(4G+1) :", ok(Wsol[0] - WG), ";  1-W = (4G^2-1)/(4G+1):", ok(1 - WG - (4*G**2 - 1)/(4*G + 1)))
# D_j as in [HB] (H3)(i): D1 = d_w F1, D2 = d_w (r2 F2), D3 = d_w (r3 F3) at (t_k, w_k); d_w = -sin w d_c
coefW = sp.together(sp.diff(K23, W))           # = r2 F2 - r3 F3 (as functions of c)
D23 = sp.factor(sp.together(sp.diff(coefW.subs(u, ub), c).subs(c, ck)))
print("  d/dc (r2F2 - r3F3) at the branch point =", D23)
print("  [so D2 - D3 = -sin(w_k) * that, nonzero for G >= 1]")

print("=== (S4) K(t_k, c) = -(G-c)(2Gc+1)^2 Q(c) / (2 G^3 (4G+1)),  Q = 2Gc^2-(2G-1)c+4G(G-1)")
Q = 2*G*c**2 - (2*G - 1)*c + 4*G*(G - 1)
print("  identity:", ok(Kb.subs(W, WG) + (G - c)*(2*G*c + 1)**2*Q/(2*G**3*(4*G + 1))))
p1 = 32*G**3 - 36*G**2 + 4*G - 1
print("  disc_c Q = -p1(G):", ok(sp.discriminant(Q, c) + p1), "; Q(c_v) = p1/(8G), c_v=(2G-1)/(4G):", ok(Q.subs(c, (2*G - 1)/(4*G)) - p1/(8*G)))
print("  Q(0) = 4G(G-1), Q(1) = (2G-1)^2, Q(-1) =", sp.factor(Q.subs(c, -1)))

print("=== (S5) roots: p1 and 2G^2-2G-1")
P1 = sp.Poly(p1, G)
print("  real roots of p1 (exact count via Sturm):", sp.Poly(p1, G).count_roots(), " isolating:", P1.intervals(eps=sp.Rational(1, 10**12)))
print("  p1(1) =", p1.subs(G, 1), ";  G_W = (1+sqrt3)/2 root of 2G^2-2G-1;  G_W^2 < 2 :", sp.simplify(((1 + sp.sqrt(3))/2)**2 - 2) < 0)

print("=== (S6) minimal polynomials of rho, q, kappa at the thresholds")
kap = sp.symbols('kappa')
Hq = sp.Poly(H.subs(X, sp.Symbol('r')), sp.Symbol('r'))
r = sp.Symbol('r')
def minpolys(polyG, name, numerics):
    # rho + 1/rho = 2G
    pr = sp.factor(sp.numer(sp.together(polyG.subs(G, (r + 1/r)/2))))
    pq = sp.factor(sp.resultant(sp.Poly(pr, r).as_expr(), H.subs(X, r), r))
    pk = sp.factor(sp.resultant(sp.Poly(pq, q).as_expr(), kap*(1 - q)**2 - (1 + q)**2, q))
    print("  [%s] rho:" % name, pr); print("        q  :", pq); print("        kappa:", pk)
    print("        number of irreducible factors (rho, q, kappa):", [len(sp.factor_list(e)[1]) for e in (pr, pq, pk)])
    for val in numerics:
        for fac, _ in sp.factor_list(pk)[1]:
            if abs(sp.N(fac.subs(kap, val), 30)) < 1e-12*max(1, abs(sp.N(sp.diff(fac, kap).subs(kap, val)))):
                print("        kappa = %s is a root of %s ; derivative there %s" % (sp.N(val, 20), fac, sp.N(sp.diff(fac, kap).subs(kap, val), 6)))
import mpmath as mp
mp.mp.dps = 40
def kap_of(Gv, sgn):
    t = sgn*mp.acosh(Gv); uu = 1 - 1/(2*Gv); qq = mp.e**(t - mp.asinh(uu)); return ((1 + qq)/(1 - qq))**2
GWv = (1 + mp.sqrt(3))/2; G1v = mp.findroot(lambda x: 32*x**3 - 36*x**2 + 4*x - 1, 1.03)
kW, k1, k2 = kap_of(GWv, -1), kap_of(G1v, -1), kap_of(G1v, 1)
print("  numerics: kappa_W = %s, kappa_1 = %s, kappa_2 = %s" % (mp.nstr(kW, 25), mp.nstr(k1, 25), mp.nstr(k2, 25)))
minpolys(2*G**2 - 2*G - 1, 'kappa_W', [sp.Float(str(kW), 40)])
minpolys(p1, 'kappa_1,2', [sp.Float(str(k1), 40), sp.Float(str(k2), 40)])

print("=== (S7) auxiliary exact facts")
print("  p1(26/25) =", p1.subs(G, sp.Rational(26, 25)), " (>0);  p1'(G) = 96G^2-72G+4 >= 28 for G >= 1")
cs_ = sp.Poly(4*G**4 - 8*G**2 + 4*G - 1, G)
iv = [I for I in cs_.intervals(eps=sp.Rational(1, 10**9)) if I[0][0] > 1]
print("  c_* (root > 1 of 4c^4-8c^2+4c-1) in", iv, "; G_1 < 26/25 < c_* < 23/20 < G_W:",
      iv[0][0][0] > sp.Rational(26, 25), iv[0][0][1] < sp.Rational(23, 20), sp.Rational(23, 20) < (1 + sp.sqrt(3))/2)
r_ = sp.Symbol('rho', positive=True)
print("  W in rho:", sp.factor(sp.together(WG.subs(G, (r_ + 1/r_)/2))), ";  1-W:", sp.factor(sp.together((1 - WG).subs(G, (r_ + 1/r_)/2))))
print("  kappa >= 3+2sqrt2 <=> q >= sqrt2-1 <=> t_max >= 0 :", sp.simplify(((1 + (sp.sqrt(2) - 1))/(1 - (sp.sqrt(2) - 1)))**2 - (3 + 2*sp.sqrt(2))) == 0)
print("  closed forms sqrt(kappa) = coth((v-t)/2):")
co = lambda x: mp.cosh(x)/mp.sinh(x)
print("    kappa_W:", mp.nstr(co((mp.asinh((3 - mp.sqrt(3))/2) + mp.acosh(GWv))/2)**2, 25))
print("    kappa_1:", mp.nstr(co((mp.asinh(1 - 1/(2*G1v)) + mp.acosh(G1v))/2)**2, 25))
print("    kappa_2:", mp.nstr(co((mp.asinh(1 - 1/(2*G1v)) - mp.acosh(G1v))/2)**2, 25))
phi = (1 + sp.sqrt(5))/2
print("    kappa_c = coth^2(asinh(1/2)/2) = phi^6 = 9+4sqrt5:", sp.simplify(sp.expand(phi**6) - 9 - 4*sp.sqrt(5)) == 0,
      sp.simplify(sp.exp(sp.asinh(sp.Rational(1, 2))) - phi) == 0)
