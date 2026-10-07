"""Proposition 41 (the spiral outside its regime): exact identities at the spiral point.

Raw active slacks f = (l_1, u_2, u_3) of one atom at frequency omega with constant C at exponent t and
twist u: a_m = g (T_m(G) - cos m omega), b_m = 2u sin(omega) sin(m omega), g = 1 - cos^2 omega - u^2,
l_1 = a_1 - b_1 + C sinh t, u_m = a_m + b_m - C sinh(mt), sinh(mt) = U_{m-1}(G) sinh t.
At the spiral point cos(omega) = -1/(2G), u = 1 - 1/(2G), with W(G) = 2(1 + 2G - 2G^2)/(4G + 1):
  (1) y = (1, W/U_1(G), (1 - W)/U_2(G)) annihilates df/dC and df/domega;
  (2) d/domega (u_2/U_1(G) - u_3/U_2(G)) = -sin(omega)(4G + 1)(2G^2 + 1)/(2G^3 (2G + 1)).
All quantities are rational functions of G and sin(omega) with sin^2 omega = 1 - 1/(4G^2)."""
import sympy as sp

G, w, C, sht = sp.symbols('G omega C sht', positive=True)
u = 1 - 1 / (2 * G)
c = sp.cos(w)
g = 1 - c ** 2 - u ** 2
T = lambda m: sp.chebyshevt(m, G)
U = lambda m: sp.chebyshevu(m, G)
a = lambda m: g * (T(m) - sp.cos(m * w))
b = lambda m: 2 * u * sp.sin(w) * sp.sin(m * w)
f = [a(1) - b(1) + C * sht, a(2) + b(2) - C * U(1) * sht, a(3) + b(3) - C * U(2) * sht]
W = 2 * (1 + 2 * G - 2 * G ** 2) / (4 * G + 1)
y = [1, W / U(1), (1 - W) / U(2)]


def at_spiral(e):
    e = sp.expand_trig(sp.expand(e))
    sw = sp.Symbol('sw', positive=True)
    e = e.subs(sp.sin(w), sw).subs(sp.cos(w), -1 / (2 * G))
    e = sp.expand(e)
    e = sp.Poly(e, sw)
    red = 0
    for (k,), co in e.terms():
        red += co * (1 - 1 / (4 * G ** 2)) ** (k // 2) * sw ** (k % 2)
    return sp.simplify(red), sw


def main():
    dC = sum(yi * sp.diff(fi, C) for yi, fi in zip(y, f))
    ok1 = sp.simplify(dC) == 0
    dw, sw = at_spiral(sum(yi * sp.diff(fi, w) for yi, fi in zip(y, f)))
    ok2 = sp.simplify(dw) == 0
    dd, sw = at_spiral(sp.diff(f[1], w) / U(1) - sp.diff(f[2], w) / U(2))
    target = -sw * (4 * G + 1) * (2 * G ** 2 + 1) / (2 * G ** 3 * (2 * G + 1))
    ok3 = sp.simplify(dd - target) == 0
    print('OK    ' if ok1 else 'FAIL  ', 'y annihilates df/dC')
    print('OK    ' if ok2 else 'FAIL  ', 'y annihilates df/dw at the spiral point')
    print('OK    ' if ok3 else 'FAIL  ', 'derivative difference = -sin w (4G+1)(2G^2+1)/(2G^3(2G+1))')
    if not (ok1 and ok2 and ok3):
        raise SystemExit('FAILED')
    print('ALL OK')


if __name__ == '__main__':
    main()
