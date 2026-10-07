"""Proposition 42 (threshold strip) and Theorem 2, exact.

Adding the minimizer constraint to l_1 gives int h dnu >= 0 with
  h(c) = (1 - c^2 - u^2)(E - c) - 2u(1 - c^2), u = sinh v, E = e^{|t|} (varrho in the paper).
At u = sinh(|t| + log phi) and c = 1 - 2x, E = 1 + y, the polynomial -4 phi^2 E^2 h is the cubic
A3 x^3 + A2 x^2 + A1 x + A0 of the proof; the script prints A3, A0 and the coefficients of
4 A1 A3 - A2^2 in y (all nonnegative), and checks dh/du <= 0."""
import sympy as sp
c, s, E, x, y = sp.symbols('c s E x y')
r5 = sp.sqrt(5); phi = (1 + r5)/2
h = (1 - c**2 - s**2)*(E - c) - 2*s*(1 - c**2)
# monotonicity in s:
print("dh/ds =", sp.factor(sp.diff(h, s)), "  (= -2 s (E - c) - 2(1-c^2) <= 0 for s>=0, c in [-1,1], E>=1)")
sE = (phi*E - 1/(phi*E))/2
P = sp.expand(-h.subs(s, sE)*4*phi**2*E**2)
P = sp.expand(sp.simplify(P))
Pxy = sp.expand(P.subs({c: 1 - 2*x, E: 1 + y}))
Pxy = sp.expand(sp.radsimp(Pxy))
poly = sp.Poly(Pxy, x, y)
print("monomials of -h*4phi^2E^2 in x=(1-c)/2, y=E-1:")
neg = []
for mon, co in zip(poly.monoms(), poly.coeffs()):
    co = sp.nsimplify(sp.radsimp(co))
    print("  x^%d y^%d :" % mon, co, "  ~", sp.N(co, 8))
    if sp.N(co) < 0: neg.append((mon, co))
print("negative coefficients:", neg)
# --- positivity: P = A3 x^3 + A2 x^2 + A1 x + A0 ;  P = A3 x (x + A2/(2A3))^2 + (A1 - A2^2/(4A3)) x + A0
Px = sp.Poly(Pxy, x)
A3, A2, A1, A0 = [sp.expand(sp.radsimp(co)) for co in Px.all_coeffs()]
print("A3 =", sp.factor(A3)); print("A0 =", sp.factor(A0))
G = sp.expand(sp.radsimp(4*A1*A3 - A2**2))
Gp = sp.Poly(G, y)
print("4 A1 A3 - A2^2 coefficients in y (must be >= 0):")
for k, co in zip(range(Gp.degree(), -1, -1), Gp.all_coeffs()):
    print(k, co)
# exact sign check of coefficients a + b sqrt5 (a,b rational): a + b sqrt5 >= 0  certified exactly
def nonneg(co):
    co = sp.expand(co); b_ = co.coeff(r5); a_ = sp.expand(co - b_*r5)
    assert a_.is_rational and b_.is_rational
    if a_ >= 0 and b_ >= 0: return True
    if a_ >= 0 and b_ < 0: return a_**2 >= 5*b_**2
    if a_ < 0 and b_ >= 0: return 5*b_**2 >= a_**2
    return False
assert all(nonneg(sp.radsimp(co)) for co in Gp.all_coeffs())
assert all(nonneg(sp.radsimp(co)) for co in sp.Poly(A0, y).all_coeffs())
assert all(nonneg(sp.radsimp(co)) for co in sp.Poly(A3, y).all_coeffs())
print("EXACT: A3 > 0, A0 >= 0 (A0 > 0 for y > 0), 4A1A3 - A2^2 >= 0 for y >= 0  =>  -h >= 0, and -h > 0 for E > 1.")
# --- explicit bound s_1(E): smallest s with max_{[-1,1]} h <= 0 ; double root of the cubic h(c)
D = sp.factor(sp.discriminant(sp.expand(h), c))
print("disc_c h =", D)
