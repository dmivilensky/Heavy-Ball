"""The case t = 0 of Proposition 42 and the three-cycle, exact.

At t = 0: a_m = g(1 - T_m(c)), b_m = 2u(1 - c^2) U_{m-1}(c), g = 1 - c^2 - u^2.
  (i)  a_1 - b_1 = (1 - c)[(1 - 2u) - (c + u)^2], hence h_1(c) = -(1 - c)(c + 1/2)^2 at u = 1/2;
  (ii) the three-cycle omega = 2 pi/3, u = 1/2, C = 0 is feasible: |b_m| <= a_m for all m
       (a_m, b_m depend on m mod 3 only), and g = 1/2."""
import sympy as sp
c, s = sp.symbols('c s')
g = 1 - c**2 - s**2
T = lambda m: sp.chebyshevt(m, c); U = lambda m: sp.chebyshevu(m, c)
a = lambda m: g*(1 - T(m)); b = lambda m: 2*s*(1 - c**2)*U(m - 1)
print("(i) a_1 - b_1 - (1-c)((1-2s)-(c+s)^2) =", sp.expand(a(1) - b(1) - (1 - c)*((1 - 2*s) - (c + s)**2)))
w = 2*sp.pi/3
for m in range(1, 4):
    am = (sp.sin(w)**2 - sp.Rational(1, 4))*(1 - sp.cos(m*w)); bm = 2*sp.Rational(1, 2)*sp.sin(w)*sp.sin(m*w)
    print("(ii) m=%d: a_m=%s, b_m=%s, a_m-|b_m|=%s" % (m, sp.nsimplify(am), sp.nsimplify(bm), sp.nsimplify(am - abs(bm))))
print("     (a_m, b_m depend on m mod 3 only)")
