"""Family P2 in exact arithmetic over Q(G)[u]/(P2): perfect-square identity of the dual quartic,
the atom polynomial c^2 + a c + b, the first moment, C sinh t, and the slack u_5 (Lemma 46, Definition 45)."""
import sympy as s, time, json
G,u,c,W=s.symbols('G u c W')
P2=256*G**4*u**2-256*G**4-512*G**3*u**3+768*G**3*u-32*G**2*u**4-480*G**2*u**2-16*G*u**3-u**4-16*u**2
p=s.Poly(P2,u,domain='QQ(G)')
def red(e):
 n,d=s.fraction(s.cancel(e))
 return s.rem(s.Poly(n,u,domain='QQ(G)')*s.invert(s.Poly(d,u,domain='QQ(G)'),p),p).as_expr()
def mul(a,b):return red(a*b)
def div(a,b):return red(a/b)
k=1-c*c
M=8*G**3+4*G**2*c-G+(1-W)*(4*G*c*c-3*G)-W*(G+c)
N=W*(4*G*G-1)+4*(1-W)*G*(G+c)
co=s.Poly(-((k-u*u)*M-2*u*k*N),c).all_coeffs();A,B,C,D,E=co
e1=s.cancel((8*A*A*D-B*(4*A*C-B*B))/(4*G*G-W));e2=s.expand(64*A**3*E-(4*A*C-B*B)**2)
lin=s.Poly(s.subresultants(e2,e1,W)[-2],W)
ww=div(-lin.coeff_monomial(1),lin.coeff_monomial(W))
print('W computed',flush=True)
aa=div(red(B.subs(W,ww)),2*red(A.subs(W,ww)))
bb=div(red(C.subs(W,ww))-mul(red(A.subs(W,ww)),mul(aa,aa)),2*red(A.subs(W,ww)))
print('atom polynomial computed',flush=True)
mom=[s.Integer(1),u] # temporary second item overwritten
# Reduce any polynomial of c modulo c^2+aa*c+bb over the algebraic field.
def coeffs(expr):
 pol=s.Poly(expr,c); v=[s.Integer(0),s.Integer(0)]
 for coeff in pol.all_coeffs():
  v=[red(coeff-mul(bb,v[1])),red(v[0]-mul(aa,v[1]))]
 return v
def lf(m):return ((1-c*c-u*u)*(s.chebyshevt(m,G)-s.chebyshevt(m,c))-2*u*(1-c*c)*s.chebyshevu(m-1,c))/s.chebyshevu(m-1,G)
def uf(m):return ((1-c*c-u*u)*(s.chebyshevt(m,G)-s.chebyshevt(m,c))+2*u*(1-c*c)*s.chebyshevu(m-1,c))/s.chebyshevu(m-1,G)
l=coeffs(lf(1));v2=coeffs(uf(2))
mc=div(-l[0]-v2[0],l[1]+v2[1]);d=red(-l[0]-mul(l[1],mc))
print('moment and C*sinh(t) computed',flush=True)
v5=coeffs(uf(5));sl5=red(v5[0]+mul(v5[1],mc)-d)
print('U5 expression computed',flush=True)

assert red(e1.subs(W,ww))==0
assert red(e2.subs(W,ww))==0
ac=red(A.subs(W,ww))
expected=s.Poly((c*c+aa*c+bb)**2,c).all_coeffs()
for actual,want in zip(co,expected):
    assert red(actual.subs(W,ww)-mul(ac,want))==0
assert red(l[0]+mul(l[1],mc)+d)==0
assert red(v2[0]+mul(v2[1],mc)-d)==0
print('Perfect-square and primal identities: OK',flush=True)

from pathlib import Path
import pickle
Path('p2_field.pkl').write_bytes(pickle.dumps(dict(P=P2,W=ww,a=aa,b=bb,mc=mc,d=d,s5=sl5)))
print('field cached')
