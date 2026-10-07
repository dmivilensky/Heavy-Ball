"""Family P1: the double-root condition disc_c R1|_{W=0} = 0 and its main factor, expressed as a
polynomial P1(cosh t, sinh v) (equation p1system)."""
import sympy as sp, pickle
from common import *
from cs_form import to_cs, cs, ss
d = pickle.load(open('thr.pkl', 'rb'))
R1 = sp.Poly(sp.expand(d['R1']), c)
R10 = sp.Poly(R1.as_expr().subs(W, 0), c)
print("R1|W=0 coefficients:", [sp.factor(x) for x in R10.all_coeffs()])
D1 = sp.factor(sp.discriminant(R10.as_expr(), c))
fl = sp.factor_list(D1)
for f, k in fl[1]: print(sp.degree(f, X), sp.degree(f, Y), k, str(f)[:120])
main = [f for f, k in fl[1] if sp.degree(f, X) > 2][0]
dx = sp.degree(main, X)//2; dy = sp.degree(main, Y)//2
P1 = to_cs(main, dx, dy); print("one-wave curve P1(ch t, sh v) =", P1)
pickle.dump({'D1': main, 'P1': P1, 'R10': R10.as_expr()}, open('D1.pkl', 'wb'))
