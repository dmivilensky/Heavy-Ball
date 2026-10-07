"""Exact tensor Bernstein positivity certificates on rational boxes (used for equation p1slackpoly).

A polynomial in n variables is rewritten exactly in the tensor Bernstein basis of its multidegree on a
box prod [a_i, b_i] with rational endpoints. If every Bernstein coefficient is positive, the polynomial
is positive on the box, because the Bernstein basis polynomials are nonnegative and sum to one.
Otherwise the box is bisected. All arithmetic is exact (fractions.Fraction)."""
from fractions import Fraction as Fr
from math import comb
import itertools, sympy as sp

def to_dense(expr, gens):
    P = sp.Poly(sp.expand(expr), *gens)
    deg = [P.degree(g) for g in gens]
    D = {}
    for mon, co in P.terms():
        D[mon] = Fr(int(sp.numer(co)), int(sp.denom(co)))
    return D, deg

def _axis_affine(D, axis, a, h, n):
    """substitute x_axis = a + h*s  (exact), D: dict monomial->Fraction"""
    out = {}
    for mon, co in D.items():
        i = mon[axis]
        for j in range(i + 1):
            m2 = list(mon); m2[axis] = j; m2 = tuple(m2)
            out[m2] = out.get(m2, 0) + co*comb(i, j)*a**(i - j)*h**j
    return out

def _axis_bern(D, axis, n):
    """power basis in s (on [0,1]) -> Bernstein basis of degree n along axis"""
    out = {}
    for mon, co in D.items():
        j = mon[axis]
        for i in range(j, n + 1):
            m2 = list(mon); m2[axis] = i; m2 = tuple(m2)
            out[m2] = out.get(m2, 0) + co*Fr(comb(i, j), comb(n, j))
    return out

def bernstein_coeffs(D, deg, box):
    for ax, ((a, b), n) in enumerate(zip(box, deg)):
        D = _axis_affine(D, ax, Fr(a), Fr(b) - Fr(a), n)
        D = _axis_bern(D, ax, n)
    return D

def all_index(deg):
    return itertools.product(*[range(n + 1) for n in deg])

def certify_positive(expr, gens, box, maxdepth=40, verbose=False):
    """returns (True, number_of_boxes) if expr > 0 on box (closed), proven by exact Bernstein coefficients"""
    D, deg = to_dense(expr, gens)
    stack = [(tuple((Fr(a), Fr(b)) for a, b in box), 0)]
    nbox = 0; maxd = 0
    while stack:
        B, d = stack.pop()
        Bc = bernstein_coeffs(D, deg, B)
        mn = min(Bc.get(I, Fr(0)) for I in all_index(deg))
        if mn > 0:
            nbox += 1; maxd = max(maxd, d); continue
        # a corner value <= 0 refutes positivity
        for corner in itertools.product(*[(0, n) for n in deg]):
            if Bc.get(corner, Fr(0)) <= 0:
                return False, ('corner value <= 0', B, corner)
        if d >= maxdepth:
            return False, ('maxdepth', B)
        ax = max(range(len(B)), key=lambda i: (B[i][1] - B[i][0])*max(deg[i], 1))
        a, b = B[ax]; m = (a + b)/2
        B1 = list(B); B1[ax] = (a, m); B2 = list(B); B2[ax] = (m, b)
        stack += [(tuple(B1), d + 1), (tuple(B2), d + 1)]
    return True, (nbox, maxd)
