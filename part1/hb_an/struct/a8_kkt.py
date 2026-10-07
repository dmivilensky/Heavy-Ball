"""Optimality systems of the families P2 and P1 (Theorem 12, equations p1system and p2system), exact.

R(c) is the numerator of the symbol L1 + W U2 + (1 - W) U3, of degree five in c, and
R = (2Xc - X^2 - 1) R1(c). Family P2: R1 is a constant times a square, which gives two equations
E1 = E2 = 0 in (X, Y, W). Family P1 (W = 0): R1 has a double root, disc_c R1 = 0."""
import sympy as sp, pickle, time
from common import *
d = pickle.load(open('thr.pkl', 'rb'))
R1 = sp.Poly(sp.expand(d['R1']), c)
A = R1.all_coeffs()            # A4..A0, linear in W
A4, A3, A2, A1, A0 = [sp.expand(x) for x in A]
# perfect-square conditions (multiplied out, polynomial):
#   al = A3/(2A4), be = (A2 - A4 al^2)/(2 A4)  ->  A1 = 2 A4 al be,  A0 = A4 be^2
E1 = sp.expand(8*A4**2*A1 - A3*(4*A4*A2 - A3**2))            # 8A4^2 * (A1 - 2A4 al be)
E2 = sp.expand(64*A4**3*A0 - (4*A4*A2 - A3**2)**2)           # 64A4^3 * (A0 - A4 be^2)
t0 = time.time()
E1f = sp.factor(E1); E2f = sp.factor(E2)
print("E1 =", E1f); print("E2 =", E2f); print(time.time() - t0)
pickle.dump({'E1': E1f, 'E2': E2f}, open('E12.pkl', 'wb'))
