"""Minimal polynomials of q and kappa at the switching points kappa_W, kappa_1, kappa_2 and at kappa_c
(Table 1), exact. On the branch: cosh t = G, sinh v = 1 - 1/(2G), q = X/Y."""
import sympy as sp
c, X, Y, q, k = sp.symbols('c X Y q k')
def minpoly_q(pc, cnum):
    s = 1 - 1/(2*c)
    e1 = X**2 - 2*c*X + 1; e2 = sp.numer(sp.together((Y**2 - 2*s*Y - 1)))
    e3 = q*Y - X
    r1 = sp.resultant(e3, e2, Y)                  # in q, X, c
    r2 = sp.resultant(r1, e1, X)                  # in q, c
    r3 = sp.resultant(r2, pc, c)                  # in q
    fl = sp.factor_list(r3)[1]
    return fl
def to_kappa(pq):
    # kappa = ((1+q)/(1-q))^2  <=>  q = (sqrt k - 1)/(sqrt k + 1); use K = sqrt(kappa)
    K = sp.symbols('K')
    e = sp.numer(sp.together(pq.subs(q, (K - 1)/(K + 1))))
    e = sp.resultant(e, K**2 - k, K)
    return sp.factor_list(e)[1]
for name, pc in (("kappa_W", 2*c**2 - 2*c - 1), ("kappa_1,2", 32*c**3 - 36*c**2 + 4*c - 1), ("kappa_c", c - 1)):
    fl = minpoly_q(pc, None)
    print(name)
    for f, mult in fl:
        if sp.degree(f, q) > 0:
            rts = [sp.N(r_, 20) for r_ in sp.Poly(f, q).real_roots() if 0 < r_ < 1]
            if rts:
                print("   q-factor (deg %d):" % sp.degree(f, q), f if sp.degree(f, q) <= 12 else str(f)[:200], " roots in (0,1):", rts)
                for g_, m_ in to_kappa(f):
                    rk = [sp.N(r_, 20) for r_ in sp.Poly(g_, k).real_roots() if r_ > 1]
                    if rk and sp.degree(g_, k) > 0: print("      kappa-factor (deg %d): %s ; roots>1: %s" % (sp.degree(g_, k), g_ if sp.degree(g_, k) <= 12 else str(g_)[:200], rk))
