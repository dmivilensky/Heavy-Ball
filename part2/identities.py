"""Polynomial identities of the general theory (Sections general and regions, Appendix general).

Notation of the paper: x = alpha mu, y = alpha L, A = 1 + beta - (x + y)/2, B = (y - x)/2,
a(z) = z^2 - (1 + beta - x) z + beta, b(z) = z^2 - (1 + beta - y) z + beta,
Xi(c) = (B^2 - (1-beta)^2 - A^2 + 2A(1+beta)c - 4 beta c^2) / (2B(1-beta)),
L(c) = 4 beta c^2 + 2(y - 1 - 2beta - beta^2 + beta x)c + (2 - 2x + xy + 2beta^2 - 2beta y).
Every check is an exact symbolic identity (sympy) unless marked as a random-point evaluation.
"""
import cmath
import math
import random
import sympy as sp

x, y, b, c, r, w, t, h, q = sp.symbols('x y beta c r omega t h q', real=True)
A = 1 + b - (x + y) / 2
B = (y - x) / 2
a = lambda z: z ** 2 - (1 + b - x) * z + b
bb = lambda z: z ** 2 - (1 + b - y) * z + b
Xi = (B ** 2 - (1 - b) ** 2 - A ** 2 + 2 * A * (1 + b) * c - 4 * b * c ** 2) / (2 * B * (1 - b))
Lc = 4 * b * c ** 2 + 2 * (y - 1 - 2 * b - b ** 2 + b * x) * c + (2 - 2 * x + x * y + 2 * b ** 2 - 2 * b * y)
H = (-4 * b * r ** 3 * c ** 3 - 2 * r ** 2 * (r ** 2 * y - r ** 2 - b * r ** 2 + b * x - 3 * b - b ** 2) * c ** 2
     - r * (r ** 4 + r ** 2 * x * y - 2 * r ** 2 * x - 2 * r ** 2 * y + 3 * r ** 2 + b ** 2 * r ** 2 + 3 * b ** 2
            - b * r ** 2 * x - b * r ** 2 * y + 2 * b * r ** 2 - b * x - b * y + 2 * b) * c
     + r ** 4 * (1 + y - x) + r ** 2 * (1 - x) * (1 - y) - 2 * b * r ** 2 * y + b ** 2 * r ** 2 + b ** 2)
Kmu = (4 * b * r ** 3 * c ** 2 + 2 * r * (r ** 3 * (y - 1) - b * r ** 3 + b * r ** 2 - 2 * b * r - b ** 2) * c
       + r ** 5 + r ** 4 * y - r ** 3 * y + r ** 3 - b * r ** 4 + b * r ** 3 - b * r ** 2 * y - b * r ** 2 - b * r * y + b * r
       + b ** 2 * r ** 2 + b ** 2)
KL = (-4 * b * r ** 3 * c ** 2 + 2 * r ** 2 * (b ** 2 + 3 * b + b * r - b * x - r ** 3) * c
      + r ** 5 - r ** 4 * x + 2 * r ** 4 + r ** 3 * x - r ** 3 + b * r ** 4 - b * r ** 3 + b * r ** 2 * x - 3 * b * r ** 2
      + b * r * x - b * r - b ** 2 * r ** 2 - 2 * b ** 2 * r + b ** 2)

FAILED = []


def check(name, expr):
    ok = sp.simplify(sp.expand(expr)) == 0
    print(('OK    ' if ok else 'FAIL  ') + name)
    if not ok:
        FAILED.append(name)


def main():
    check('normalization: A + B = 1 + beta - x, A - B = 1 + beta - y', (A + B - (1 + b - x)) ** 2 + (A - B - (1 + b - y)) ** 2)
    z = sp.exp(t + sp.I * w); zc = sp.exp(t - sp.I * w)
    m = (z + b / z - A) / B; mc = (zc + b / zc - A) / B
    check('eq. factor: (1 - m)(1 + conj m) = -a(z) conj b(z) / (B^2 |z|^2)',
          sp.expand((1 - m) * (1 + mc) + a(z) * bb(zc) / (B ** 2 * z * zc)))
    z0 = sp.exp(sp.I * w); z0c = sp.exp(-sp.I * w)
    m0 = (z0 + b / z0 - A) / B; m0c = (z0c + b / z0c - A) / B
    g0 = 1 - m0 * m0c; s0 = (m0 - m0c) / (2 * sp.I)
    check('Lemma angle: g sin(w) = 2 s Xi(cos w) at t = 0',
          sp.expand((g0 * sp.sin(w) - 2 * s0 * Xi.subs(c, sp.cos(w))).rewrite(sp.exp)))
    check('Lemma angle: s = (1 - beta) sin(w) / B at t = 0', sp.expand((s0 - (1 - b) * sp.sin(w) / B).rewrite(sp.exp)))
    check('Lemma angle: Xi(1) = (B^2 - (1 + beta - A)^2)/(2B(1 - beta))',
          Xi.subs(c, 1) - (B ** 2 - (1 + b - A) ** 2) / (2 * B * (1 - b)))
    check('Lemma angle: Xi(-1) = (B^2 - (1 + beta + A)^2)/(2B(1 - beta))',
          Xi.subs(c, -1) - (B ** 2 - (1 + b + A) ** 2) / (2 * B * (1 - b)))
    zr = r * sp.exp(sp.I * w); zrc = r * sp.exp(-sp.I * w)
    Hfull = ((1 - zr) * a(zr) * bb(zrc) + (1 - zrc) * a(zrc) * bb(zr)) / 2
    check('Lemma cubic: H_r(cos w) = Re[(1 - z) a(z) conj b(z)], z = r e^{iw}',
          sp.simplify(sp.expand((Hfull - H.subs(c, sp.cos(w))).rewrite(sp.exp))))
    check('Lemma cubic: H_1 = (1 - c) L(c)', H.subs(r, 1) - (1 - c) * Lc)
    check('eq. E: L(c) = (y - x)(1 - beta)(1 + c - Xi(c))', Lc - (y - x) * (1 - b) * (1 + c - Xi))
    check('Lemma cubic: L(-1) = (2 + 2beta - x)(2 + 2beta - y)', Lc.subs(c, -1) - (2 + 2 * b - x) * (2 + 2 * b - y))
    check('Lemma cubic: L(1) = 2(1 - beta)(y - x) + x y', Lc.subs(c, 1) - (2 * (1 - b) * (y - x) + x * y))
    kp = x / y
    gtd = x ** 2 - 2 * (b - c + kp * (1 - b * c)) * x + 2 * kp * (1 - c) * (1 + b ** 2 - 2 * b * c)
    check('eq. gtd: GTD Theorem 3.5 polynomial = (mu/L) L(c)', sp.simplify(gtd - kp * Lc))
    hh = sp.Symbol('hh')
    check('Corollary C(b): L(c) - (y-x)(1-beta)(1+c-h) = -(y-x)(1-beta)(Xi(c)-h)',
          (Lc - (y - x) * (1 - b) * (1 + c - hh)) + (y - x) * (1 - b) * (Xi - hh))
    check('Lemma Kmu: H_r = (1 - c) K_mu when a(r) = 0', sp.expand(r ** 2 * (H.subs(x, (1 - r) * (1 - b / r)) - (1 - c) * Kmu)))
    check('Lemma Kmu: H_r = (1 + c) K_L when b(-r) = 0', sp.expand(r ** 2 * (H.subs(y, (1 + r) * (1 + b / r)) - (1 + c) * KL)))
    check('Lemma Kmu: K_mu(-1) = (beta + r^2)(1 + r) b(-r)', sp.expand(Kmu.subs(c, -1) - (b + r ** 2) * (1 + r) * bb(-r)))
    check('Lemma Kmu: K_L(1) = (beta + r^2)(1 - r) a(r)', sp.expand(KL.subs(c, 1) - (b + r ** 2) * (1 - r) * a(r)))
    check('Lemma Kmu: leading coefficients 4 beta r^3 and -4 beta r^3',
          (sp.Poly(Kmu, c).LC() - 4 * b * r ** 3) ** 2 + (sp.Poly(KL, c).LC() + 4 * b * r ** 3) ** 2)
    check('Theorem B(b): a(z) - b(z) = (x - y) z', sp.expand(a(r) - bb(r) - (x - y) * r))
    # Lemma near1: gamma(1)/2 + gamma'(1) = L(1)/(2 y^2), gamma = a/b
    zz = sp.Symbol('zz')
    gam = a(zz) / bb(zz)
    check("Lemma near1: gamma(1)/2 + gamma'(1) = L(1)/(2y^2)",
          sp.simplify(gam.subs(zz, 1) / 2 + sp.diff(gam, zz).subs(zz, 1) - Lc.subs(c, 1) / (2 * y ** 2)))
    # Corollary D(a): Polyak's tuning
    u0 = (1 / q - q) / 2
    sub = {x: (1 - q) ** 2, y: (1 + q) ** 2, b: q ** 2}
    Xp = sp.simplify(Xi.subs(sub))
    check('Corollary D(a): Xi = (1 - c^2 - u0^2)/(2 u0)', sp.simplify(Xp - (1 - c ** 2 - u0 ** 2) / (2 * u0)))
    check('Corollary D(a): 2u0(1 + c - Xi) = (c + u0)^2 + 2u0 - 1', sp.simplify(2 * u0 * (1 + c - Xp) - ((c + u0) ** 2 + 2 * u0 - 1)))
    check('Corollary D(a): value at c = -1/2 equals (u0 + 3/2)(u0 - 1/2)',
          sp.expand(((c + u0) ** 2 + 2 * u0 - 1).subs(c, -sp.Rational(1, 2)) - (u0 + sp.Rational(3, 2)) * (u0 - sp.Rational(1, 2))))
    check('Corollary D(a): a(z) = (z - q)^2, b(z) = (z + q)^2',
          sp.expand(a(r).subs(sub) - (r - q) ** 2) ** 2 + sp.expand(bb(r).subs(sub) - (r + q) ** 2) ** 2)
    # Corollary D(b): gradient descent
    check('Corollary D(b): L(c) = 2(y - 1)c + 2 - 2x + xy at beta = 0', Lc.subs(b, 0) - (2 * (y - 1) * c + 2 - 2 * x + x * y))
    check('Corollary D(b): K_mu = r^3{2r(y-1)c + r^2 + ry - y + 1} at beta = 0',
          Kmu.subs(b, 0) - r ** 3 * (2 * r * (y - 1) * c + r ** 2 + r * y - y + 1))
    check('Corollary D(b): K_mu(-1) = r^3(1 + r)(1 + r - y) at beta = 0', Kmu.subs({b: 0, c: -1}) - r ** 3 * (1 + r) * (1 + r - y))
    check('Corollary D(b): K_mu(1) = r^3{(1-r)^2 + y(3r-1)} at beta = 0', Kmu.subs({b: 0, c: 1}) - r ** 3 * ((1 - r) ** 2 + y * (3 * r - 1)))
    check('Corollary D(b): K_mu(1) - 4r^5 = r^3 (1+r-y)(1-3r) at beta = 0',
          Kmu.subs({b: 0, c: 1}) - 4 * r ** 5 - r ** 3 * (1 + r - y) * (1 - 3 * r))
    check('Corollary D(b): K_L(-1) = r^3{3r^2 + 2r - 1 + x(1-r)} at beta = 0',
          KL.subs({b: 0, c: -1}) - r ** 3 * (3 * r ** 2 + 2 * r - 1 + x * (1 - r)))
    check('Corollary D(b): K_L(-1) - 4r^5 = r^3 (1-r)(x - 1 + r) at beta = 0',
          KL.subs({b: 0, c: -1}) - 4 * r ** 5 - r ** 3 * (1 - r) * (x - 1 + r))
    # Lemma symbols: random-point evaluation of both representations of the symbol
    random.seed(1)
    worst = 0.0
    for _ in range(40):
        xv = random.uniform(0.01, 0.5); yv = random.uniform(xv + 0.1, 2.5); bv = random.uniform(0.05, 0.9)
        if yv >= 2 * (1 + bv):
            continue
        Av = 1 + bv - (xv + yv) / 2; Bv = (yv - xv) / 2
        tv = random.uniform(0.0, 0.6); rho = complex(sp.exp(tv)); wv = random.uniform(0.1, 3.0)
        lam = [random.uniform(0, 1) for _ in range(4)]
        zv = rho * complex(sp.exp(sp.I * wv))
        mv = (zv + bv / zv - Av) / Bv
        gv = 1 - abs(mv) ** 2; sv = mv.imag
        wts = [lam[n - 1] * rho.real ** (-n) for n in range(1, 5)]
        gam_minus = sum(wts[n - 1] * math.sinh(n * tv) for n in range(1, 5))
        sym = sum(wts[n - 1] * (gv * (math.cosh(n * tv) - math.cos(n * wv)) - 2 * sv * math.sin(n * wv)) for n in range(1, 5)) + gam_minus * gv
        u = cmath.exp(1j * wv) / rho.real
        a_sh = lambda uu: 1 - (1 + bv - xv) * uu + bv * uu ** 2
        b_sh = lambda uu: 1 - (1 + bv - yv) * uu + bv * uu ** 2
        Mu = sum(lam[n - 1] * (1 - u ** n) for n in range(1, 5))
        rhs = -(rho.real ** 2 / Bv ** 2) * abs(a_sh(u)) ** 2 * (Mu * b_sh(u) / a_sh(u)).real
        worst = max(worst, abs(sym - rhs))
    ok = worst < 1e-9
    print(('OK    ' if ok else 'FAIL  ') + 'Lemma symbols (a): both forms of the symbol agree at 40 random points (max diff %.1e)' % worst)
    if not ok:
        FAILED.append('Lemma symbols')
    if FAILED:
        raise SystemExit('FAILED: ' + ', '.join(FAILED))
    print('ALL OK')


if __name__ == '__main__':
    main()
