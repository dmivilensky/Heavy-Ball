"""The curves of P2 and P1 in (cosh t, sinh v), exact.

Checks: palindromy in X (t -> -t); conversion to (G, u); the restriction of P2 to the spiral
u = 1 - 1/(2G) equals -(32G^3 - 36G^2 + 4G - 1)^2 (tangency at G_1, Lemma 47(a)); the restriction of
P1 contains the factor (2G^2 - 2G - 1)^2 (contact at G_W, Proposition 43)."""
import sympy as sp, pickle
from common import X, Y
from cs_form import to_cs, cs, ss
fl = pickle.load(open('D2.pkl', 'rb'))['Res'][1]
D2 = [f for f, k in fl if sp.degree(f, X) == 8][0]
print("D2(1/X,Y)*X^8 - D2 =", sp.expand(D2.subs(X, 1/X)*X**8 - D2))
P2 = to_cs(D2, 4, 4); P2 = sp.factor(P2/16)
P1 = pickle.load(open('D1.pkl', 'rb'))['P1']; P1 = sp.factor(P1/1024)
print("P2 =", sp.expand(P2)); print("P1 =", sp.expand(P1))
sp_ = 1 - 1/(2*cs)
print("P2 on spiral:", sp.factor(sp.numer(sp.together(P2.subs(ss, sp_)))))
print("P1 on spiral:", sp.factor(sp.numer(sp.together(P1.subs(ss, sp_)))))
p10 = 448*cs**10 - 1152*cs**9 + 800*cs**8 - 1248*cs**7 + 3212*cs**6 - 2512*cs**5 + 108*cs**4 + 688*cs**3 - 345*cs**2 + 64*cs - 4
print("p10 roots in [1, oo):", sp.Poly(p10, cs).count_roots(1, sp.oo))
pickle.dump({'P1': P1, 'P2': P2}, open('P12.pkl', 'wb'))
