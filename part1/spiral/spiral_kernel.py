"""Lemma 39 (positivity of the spiral kernel) and Proposition 10, exact.

E_n(t) = cosh(nt) + 2G^2 cos(n omega), G = cosh t, cos(omega) = -1/(2G); F_n = G^(n-2) E_n.
With X = G^2 = 1 + s the claims are identities in Q[s][[r]]:
 (1) F_1 = F_2 = 0 and sum_{n>=3} F_n r^n = (2s+3) r^3 ((1-r) + s(2-r)) / D,
     D = (1+r+Xr^2)(1-2Xr+Xr^2) = D0 - s D1 - s^2 D2,
     D0 = (1-r)(1-r^3), D1 = r(2+3r^2-2r^3), D2 = r^3(2-r);
 (2) ((1-r)+s(2-r))/D = (1-r^3)^(-1) (1 + s(2-r)/(1-r)) / (1 - (s D1 + s^2 D2)/((1-r)(1-r^3))),
     and D1/(1-r), D2/(1-r), (2-r)/(1-r), 1/(1-r^3) have nonnegative coefficients;
 (3) direct cross-check for 3 <= n <= 200: F_n has nonnegative coefficients in s,
     F_n(0) = 3 if 3 | n and 0 otherwise, and the coefficient of s is positive.
All arithmetic is exact (sympy, rationals)."""
import sympy as sp
s, r, Xs = sp.symbols('s r X')
NMAX = 200
# recurrences: A_n = G^n T_n(G), B_n = G^n T_n(-1/(2G)), as polynomials in X = G^2
A = [sp.Poly(1, Xs), sp.Poly(Xs, Xs)]
B = [sp.Poly(1, Xs), sp.Poly(sp.Rational(-1, 2), Xs)]
for n in range(1, NMAX):
    A.append(2 * sp.Poly(Xs, Xs) * A[n] - sp.Poly(Xs, Xs) * A[n - 1])
    B.append(-B[n] - sp.Poly(Xs, Xs) * B[n - 1])
F = [None]
for n in range(1, NMAX + 1):
    num = A[n] + 2 * sp.Poly(Xs, Xs) * B[n]
    q, rem = sp.div(num, sp.Poly(Xs, Xs))
    assert rem.is_zero, n
    F.append(sp.Poly(q.as_expr().subs(Xs, 1 + s), s))
assert F[1].is_zero and F[2].is_zero
print('F_1 = F_2 = 0: True')
# (3) direct check
for n in range(3, NMAX + 1):
    co = F[n].all_coeffs()
    assert all(c >= 0 for c in co), n
    assert F[n].eval(0) == (3 if n % 3 == 0 else 0), n
    assert F[n].coeff_monomial(s) > 0, n
print(f'F_n >= 0 coefficientwise, F_n(0) = 3[3|n], positive s-coefficient: 3 <= n <= {NMAX}')
# (1) generating function: compare power series coefficients up to order N1
N1 = 40
X = 1 + s
D = sp.Poly(sp.expand((1 + r + X * r**2) * (1 - 2 * X * r + X * r**2)), r)
D0 = (1 - r) * (1 - r**3); Q1 = r * (2 + 3 * r**2 - 2 * r**3); Q2 = r**3 * (2 - r)
assert sp.expand(D.as_expr() - (D0 - s * Q1 - s**2 * Q2)) == 0
num = sp.expand((2 * s + 3) * r**3 * ((1 - r) + s * (2 - r)))
assert sp.expand((2 * X + 1) * r**3 * (2 * X - 1 - X * r) - num) == 0
# series division num/D in Q[s][[r]]
dco = [sp.expand(D.as_expr().coeff(r, k)) for k in range(5)]
nco = [sp.expand(num.coeff(r, k)) for k in range(N1 + 1)]
ser = []
for k in range(N1 + 1):
    v = nco[k] - sum(dco[j] * ser[k - j] for j in range(1, min(k, 4) + 1))
    ser.append(sp.expand(v / dco[0]))
for n in range(N1 + 1):
    target = F[n].as_expr() if n >= 3 else 0
    assert sp.expand(ser[n] - target) == 0, n
print(f'generating function identity: coefficients 0..{N1} match')
# (2)
assert sp.simplify(Q1 / (1 - r) - (2 * r + 2 * r**2 + 5 * r**3 + 3 * r**4 / (1 - r))) == 0
assert sp.simplify(Q2 / (1 - r) - (2 * r**3 + r**4 / (1 - r))) == 0
assert sp.simplify((2 - r) / (1 - r) - (2 + r / (1 - r))) == 0
Z = (s * Q1 + s**2 * Q2) / ((1 - r) * (1 - r**3))
lhs = ((1 - r) + s * (2 - r)) / D.as_expr()
rhs = (1 + s * (2 - r) / (1 - r)) / (1 - r**3) / (1 - Z)
assert sp.simplify(sp.together(lhs - rhs)) == 0
print('expansion identity and nonnegative factors: True')
print('ALL OK')
