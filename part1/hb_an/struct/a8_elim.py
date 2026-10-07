"""Family P2: elimination of W from E1 = E2 = 0 (resultant) and factorization; the main factor is
the curve D2(X, Y) = 0, which is the polynomial P2 of equation p2poly in (cosh t, sinh v)."""
import sympy as sp, pickle, time
from common import *
d = pickle.load(open('E12.pkl', 'rb'))
def mainfactor(e, deg):
    fs = sp.factor_list(e)[1]
    return [f for f, k in fs if sp.degree(f, W) == deg][0]
e1 = mainfactor(d['E1'], 2); e2 = mainfactor(d['E2'], 4)
t0 = time.time()
Res = sp.resultant(e1, e2, W)
print("resultant computed", time.time() - t0)
fl = sp.factor_list(Res)
print("factors (degree in X, Y, mult):")
for f, k in fl[1]:
    print(sp.degree(f, X), sp.degree(f, Y), k, str(f)[:150])
pickle.dump({'e1': e1, 'e2': e2, 'Res': fl}, open('D2.pkl', 'wb'))
