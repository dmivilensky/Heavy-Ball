"""Taylor-dominance sign certificates (Lemma 48).

A polynomial R with rational coefficients has the constant sign sgn R(zc) on [zc - delta, zc + delta]
if |R(zc)| > sum_{k >= 1} |R^(k)(zc)/k!| delta^k; the comparison is exact in rational arithmetic."""
from fractions import Fraction as Fr
import sympy as s
def taylor(p, zc):
    """coefficients of p(zc + e) in e (exact)"""
    e=s.symbols('e'); z=list(p.free_symbols)[0]
    q=s.Poly(s.expand(p.subs(z,zc+e)),e)
    return [Fr(int(s.numer(cf)),int(s.denom(cf))) for cf in reversed(q.all_coeffs())]
def dominance(p, lo, hi):
    """returns (sign, margin) if p has constant sign on [lo,hi] certified by |p0| > sum |pk| d^k; else (0,None)"""
    lo,hi=Fr(lo),Fr(hi); zc=(lo+hi)/2; d=(hi-lo)/2
    co=taylor(p, s.Rational(zc.numerator,zc.denominator))
    tail=sum(abs(ck)*d**k for k,ck in enumerate(co) if k>0)
    if abs(co[0])>tail: return (1 if co[0]>0 else -1), float((abs(co[0])-tail)/abs(co[0]))
    return 0,None
def certify_sign(p, lo, hi, maxsplit=12):
    """split interval until dominance holds on all pieces; returns sign and number of pieces"""
    stack=[(Fr(lo),Fr(hi))]; pieces=0; sg=None
    while stack:
        a,b=stack.pop(); r,mg=dominance(p,a,b)
        if r==0:
            if (b-a) < Fr(hi-lo)/2**maxsplit: return None
            m=(a+b)/2; stack+= [(a,m),(m,b)]; continue
        if sg is None: sg=r
        if r!=sg: return None
        pieces+=1
    return sg,pieces
def value(p,x):
    z=list(p.free_symbols)[0]; v=p.subs(z,s.Rational(Fr(x).numerator,Fr(x).denominator)); return v
