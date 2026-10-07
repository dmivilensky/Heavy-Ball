"""Expansions of Theorem 3 as q -> 0, exact.

Family P1: sinh v = 1 - z/2 + z^2/8 - z^3/4 + ... in z = 1/cosh t, from the regular system
(equations p1system, p1blowup); spiral: sinh v = 1 - z/2. Expansions of rho/q in q and the difference
rho_HB - rho_sp = (10 + 7 sqrt 2)/4 q^3 + O(q^4)."""
import sympy as sp, pickle, sys
sys.path.insert(0, '.')
from cs_form import cs, ss
P1 = pickle.load(open('P12.pkl', 'rb'))['P1']
z, sig = sp.symbols('z sigma')
N = 6
# one-wave: s = 1 + sum_{j>=1} d_j z^j
d = sp.symbols('d1:%d' % (N + 1))
sz = 1 + sum(d[j]*z**(j + 1) for j in range(N))
Pz = sp.expand(P1.subs(cs, 1/z)*z**10)
ser = sp.expand(sp.series(Pz.subs(ss, sz), z, 0, N + 2).removeO())
co = [ser.coeff(z, j) for j in range(0, N + 2)]
print("orders:", [sp.factor(cc) for cc in co[:3]])
sols = {}
# solve order by order
for j in range(0, N + 1):
    cj = sp.factor(co[j].subs(sols))
    if cj == 0: continue
    free = [x for x in d if x in cj.free_symbols and x not in sols]
    if not free: print("inconsistent at order", j, cj); break
    sj = sp.solve(cj, free[0])
    print("order z^%d: solve for %s:" % (j, free[0]), sj)
    sols[free[0]] = [x for x in sj if x.is_real][0] if len(sj) > 1 else sj[0]
s1w = sp.expand(sz.subs(sols)); print("one-wave: sh v =", s1w)
print("spiral  : sh v = 1 - z/2")
# rho/q = e^v = s + sqrt(1+s^2) ; z = 1/ch t = 2X/(1+X^2), X = rho ; then solve rho = q * Y(rho)
X, q, R = sp.symbols('X q R')
def rho_series(sexpr, order=4):
    s_X = sp.series(sexpr.subs(z, 2*X/(1 + X**2)), X, 0, order + 1).removeO()
    Yx = sp.series(s_X + sp.sqrt(1 + s_X**2), X, 0, order + 1).removeO()
    # rho = q Y(rho): fixed point iteration in powers of q
    r_ = sp.Integer(0)
    for _ in range(order + 1):
        r_ = sp.series(q*Yx.subs(X, r_), q, 0, order + 1).removeO()
    return sp.expand(sp.nsimplify(sp.radsimp(sp.expand(r_))))
rs = rho_series(1 - z/2); r1 = rho_series(s1w)
print("spiral   rho/q =", sp.expand(sp.radsimp(rs/q)))
print("one-wave rho/q =", sp.expand(sp.radsimp(r1/q)))
print("difference (one-wave - spiral) =", sp.factor(sp.radsimp(sp.expand(r1 - rs))))
# Theorem 2.3: rho_HB = (1+sqrt2) q - (2+3sqrt2/2) q^2 + O(q^3), rho_HB - rho_sp = (10+7sqrt2)/4 q^3 + O(q^4)
r2 = sp.sqrt(2); Pq = sp.Poly(sp.expand(r1), q); Dq = sp.Poly(sp.expand(r1 - rs), q)
ok = (sp.simplify(Pq.coeff_monomial(q) - (1 + r2)) == 0 and sp.simplify(Pq.coeff_monomial(q**2) + 2 + 3*r2/2) == 0
      and all(sp.simplify(Dq.coeff_monomial(q**j)) == 0 for j in (0, 1, 2))
      and sp.simplify(Dq.coeff_monomial(q**3) - (10 + 7*r2)/4) == 0)
print("coefficients stated in Theorem 2.3 (q -> 0):", ok)
assert ok