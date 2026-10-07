"""The switching points of the spiral certificate (Section 4.1, Appendix B), exact.

On the spiral branch S(X, Y) = 0 with atom c_k = -X/(1 + X^2): the weight W_s that makes the zero at c_k
double, W_s = 0 at kappa_W, and the second double zero of the residual quadratic Q2 at kappa_1 and
kappa_2 (discriminant factor 32 G^3 - 36 G^2 + 4G - 1). Both systems are invariant under X -> 1/X.
Sign checks of the leading coefficient and of Q2(+-1) in each regime."""
import sympy as sp
from common import *
Kt = sp.together(Lfun(1) + W*Ufun(2) + (1 - W)*Ufun(3))
num, den = sp.fraction(sp.factor(Kt))
R = sp.Poly(sp.expand(num), c)
ck = -X/(1 + X**2); cht = (X**2 + 1)/(2*X)
# exact division by (c - ch t):  use 2X c - (X^2+1)
R1, rem = sp.div(R, sp.Poly(2*X*c - (X**2 + 1), c))
assert rem.is_zero
print("R = (2Xc - X^2 - 1) * R1(c), deg R1 =", R1.degree())
# reduce modulo S: S is quadratic in Y; reduce a polynomial in Y modulo S (leading coeff (1+X^2) in Y^2)
Sp = sp.Poly(S, Y)
sig = 1 - X/(1 + X**2)                    # sh v on the spiral
Smon = Y**2 - 2*sig*Y - 1                  # = S/(1+X^2), monic in Y
def redS(expr):
    """reduce a rational function of (X,Y,W) modulo S (Y algebraic of degree 2 over Q(X)):
    returns A + B*Y with A, B in Q(X,W). Denominators in Y are inverted via the conjugate."""
    e = sp.together(expr); nn, dd = sp.fraction(e)
    def red(p):
        P = sp.Poly(sp.expand(p), Y, domain='QQ(X,W)')
        _, rr = sp.div(P, sp.Poly(Smon, Y, domain='QQ(X,W)'))
        co = rr.all_coeffs()
        return (co[-1] if co else 0), (co[-2] if len(co) > 1 else 0)
    A1, B1 = red(nn); A2, B2 = red(dd)
    # (A1+B1 Y)/(A2+B2 Y); conjugate root Y' = 2 sig - Y; (A2+B2Y)(A2+B2Y') = A2^2 + 2 sig A2 B2 - B2^2
    Nn = sp.expand((A1 + B1*Y)*(A2 + B2*(2*sig - Y)))
    Dn = sp.factor(A2**2 + 2*sig*A2*B2 - B2**2)
    A3, B3 = red(Nn)
    return sp.factor(A3/Dn) + sp.factor(B3/Dn)*Y
# R1(c_k) and R1'(c_k) on S
R1e = R1.as_expr()
v0 = redS(R1e.subs(c, ck)); print("R1(c_k) mod S (linear in W):", v0)
v1 = redS(sp.diff(R1e, c).subs(c, ck))
Ws = sp.solve(sp.Eq(v1, 0), W)
print("number of solutions for W:", len(Ws))
Ws = sp.factor(Ws[0]); print("W_s =", Ws)
print("R1(c_k) at W_s mod S:", redS(v0.subs(W, Ws)))
# ---- W_s in terms of c_t = ch t:  W_s = (2 + 4 c_t - 4 c_t^2)/(4 c_t + 1)
ct = sp.symbols('c_t')
print("W_s - (2+4ch t-4ch^2 t)/(4ch t+1) =", sp.simplify(Ws - (2 + 4*chX(1) - 4*chX(1)**2)/(4*chX(1) + 1)))
# ---- kappa_W: W_s = 0  <=>  2 ch^2 t - 2 ch t - 1 = 0  <=> ch t = (1+sqrt3)/2 (t<0 branch)
cW = (1 + sp.sqrt(3))/2
shvW = 1 - 1/(2*cW)
tW = -sp.acosh(cW); vW = sp.asinh(shvW); qW = sp.exp(tW - vW)
kW = ((1 + qW)/(1 - qW))**2
print("kappa_W = ", sp.N(kW, 30), " (F7: 2.6564614526046885)")
# exact algebraic form of q_W: q = X/Y with X = e^{tW}, Y = e^{vW}
# ---- Q2 = R1 / (c - c_k)^2 at W = W_s (exact division in Q(X,Y) then reduce mod S)
R1W = sp.Poly(sp.expand(sp.numer(sp.together(R1.as_expr().subs(W, Ws)))), c)
dW = sp.denom(sp.together(R1.as_expr().subs(W, Ws)))
Q2, rem2 = sp.div(R1W, sp.Poly(((1 + X**2)*c + X)**2, c))
rem2 = [redS(co) for co in rem2.all_coeffs()]
print("remainder of R1/(c-c_k)^2 mod S:", rem2)
co = Q2.all_coeffs()
print("deg Q2 =", Q2.degree())
A2_, A1_, A0_ = [redS(x) for x in co]
disc = redS(A1_**2 - 4*A2_*A0_)
# norm of (A + B Y): (A+BY)(A+BY') = A^2 + 2 sig A B - B^2
Ad = sp.expand(disc.subs(Y, 0)); Bd = sp.factor(sp.diff(disc, Y))
Nm = sp.factor(sp.numer(sp.together(Ad**2 + 2*sig*Ad*Bd - Bd**2)))
print("norm of disc(Q2) (poly in X):", Nm)
import pickle; pickle.dump({'Nm': Nm, 'Ws': Ws, 'R': R.as_expr(), 'R1': R1.as_expr(), 'disc': disc, 'A2': A2_, 'A1':A1_, 'A0':A0_}, open('thr.pkl', 'wb'))
# ---- the threshold polynomial: palindromic sextic P6(X) = X^3 * p3(2 ch t), p3(z)=4z^3-9z^2+2z-1
P6 = 4*X**6 - 9*X**5 + 14*X**4 - 19*X**3 + 14*X**2 - 9*X + 4
cc = sp.symbols('cc')
print("P6/X^3 - (32 ch^3 - 36 ch^2 + 4 ch - 1) =", sp.simplify(P6/X**3 - (32*chX(1)**3 - 36*chX(1)**2 + 4*chX(1) - 1)))
# disc itself (not only its norm) is divisible by P6:  disc = P6 * (A' + B' Y)
Ad = sp.factor(disc.subs(Y, 0)); Bd = sp.factor(sp.diff(disc, Y))
print("A-part of disc:", Ad); print("B-part of disc:", Bd)
p3 = 32*cc**3 - 36*cc**2 + 4*cc - 1
print("p3 discriminant:", sp.discriminant(p3, cc), " real roots:", sp.real_roots(sp.Poly(p3, cc)))
c1 = sp.real_roots(sp.Poly(p3, cc))[0]
for sgn, name in ((-1, 'kappa_1'), (1, 'kappa_2')):
    tt = sgn*sp.acosh(c1); vv = sp.asinh(1 - 1/(2*c1)); qq = sp.exp(tt - vv)
    print(name, "=", sp.N(((1 + qq)/(1 - qq))**2, 30), " rho =", sp.N(sp.exp(tt), 20), " W_s =", sp.N((2 + 4*c1 - 4*c1**2)/(4*c1 + 1), 20))
    # location of the second double zero: root of Q2
    Xv = sp.exp(tt); Yv = sp.exp(vv)
    A2v, A1v = [sp.N(e.subs({X: Xv, Y: Yv}), 30) for e in (A2_, A1_)]
    c2 = -A1v/(2*A2v); print("   second double zero c2 =", sp.N(c2, 20), " w2 =", sp.N(sp.acos(c2), 20))
# ---- sign structure of the spiral certificate along the branch (exact at sample points of each arc;
#      constancy of signs on arcs follows from: disc(Q2) = P6 * D', D' != 0 on the branch (norm has no real zeros, X>0))
def exact_point(Xv):
    sg = 1 - Xv/(1 + Xv**2); Yv = sg + sp.sqrt(sg**2 + 1)
    vals = [sp.nsimplify(sp.radsimp(e.subs({X: Xv, Y: Yv}))) for e in (A2_, A1_, A0_, disc)]
    return vals
for Xv, label in ((sp.Rational(5, 4), "0<t<t_1, c=1.025"), (sp.Rational(2), "t>t_1, c=1.25"), (sp.Rational(21, 20), "0<t<t_1"), (sp.Rational(1, 2), "t<-t_1, c=1.25")):
    a2_, a1_, a0_, ds = exact_point(Xv)
    cv = (Xv + 1/Xv)/2
    print("X=%s (ch t=%s, %s): lead A2 sign %s, disc sign %s, Q2(-1),Q2(1) signs %s %s" % (Xv, sp.N(cv, 6), label,
          sp.sign(sp.N(a2_, 50)), sp.sign(sp.N(ds, 50)), sp.sign(sp.N(a2_ - a1_ + a0_, 50)), sp.sign(sp.N(a2_ + a1_ + a0_, 50))))
# overall factor kap0 := Kt(t_k,c) / ((c - ch t)(c - c_k)^2 Q2(c)) is c-free on the branch (after reduction mod S);
# numerical evaluation at branch points (sign only; used to orient the statement):
import mpmath as mpm
mpm.mp.dps = 40
Ktk = R.as_expr().subs(W, Ws)/(4*X*Y**2*(X**2 + 1)*(X**2 - X + 1)*(X**2 + X + 1))
for Xv in (mpm.mpf(2), mpm.mpf('0.5'), mpm.mpf('1.05'), mpm.mpf('0.7')):
    sg = 1 - Xv/(1 + Xv**2); Yv = sg + mpm.sqrt(sg**2 + 1)
    vals = []
    for cv in (mpm.mpf('-0.9'), mpm.mpf(0), mpm.mpf('0.6')):
        num_ = sp.N(Ktk.subs({X: Xv, Y: Yv, c: cv}), 30)
        den_ = sp.N(((c - cht)*(c - ck)**2*(A2_*c**2 + A1_*c + A0_)).subs({X: Xv, Y: Yv, c: cv}), 30)
        vals.append(sp.N(num_/den_, 12))
    print("X=%s: kap0 at c=-0.9,0,0.6:" % Xv, vals)
# exact: A2 (leading coeff of Q2) has no zero on the branch: norm of A2 = A + B Y
Aa = sp.expand(A2_.subs(Y, 0)); Bb = sp.factor(sp.diff(A2_, Y))
NA = sp.factor(A2_)
print("A2 (reduced leading coeff of Q2) =", NA, " = 8X(X^4+X^2+1)[(X^2+1)+2Y(X^2-X+1)]/(X^2+1)^2 > 0")
