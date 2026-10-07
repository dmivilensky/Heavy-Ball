"""Proposition 10 (spiral) and Lemma 38 (lag identities on the spiral), exact.

Checks, with rho = e^t, G = cosh t, c = cos(omega):
  (1) the active system l_1 = u_2 = u_3 = 0 of one atom, its determinant, and its solution
      C = g tanh t, u = 1 - 1/(2G), g = u/G, c = -1/(2G) (equation spiralsimple);
  (2) sinh(t - log q) = 1 - 1/(2G) is equivalent to F(rho, q) = 0 (equation spiralpoly);
  (3) F(1, q) = 0 exactly at q = 1/phi, that is, at kappa_c = 9 + 4 sqrt(5);
  (4) Lemma 38: l_m = g E_{m+1}/G and u_m = g E_{m-1}/G for m = 1..12, E_1 = E_2 = 0.
All checks are identities in Q(rho, c) (sympy)."""
import sympy as sp
r,c,q,t,g,C,B,G = sp.symbols('rho c q t g C B G', real=True)
ch=lambda x:(sp.exp(x)+sp.exp(-x))/2
# --- (1) spiral lemma: L1, U2, U3 active; b_m = B*U_{m-1}(c), B = 2 Im s sin w ; rho=e^t
U=lambda m:sp.chebyshevu(m-1,c); T=lambda m:sp.chebyshevt(m,c)
chm=lambda m:(r**m+r**-m)/2; shm=lambda m:(r**m-r**-m)/2
eqs=[g*(chm(1)-T(1))+C*shm(1)-B*U(1),          # L1: a1 + x1 = 0, x1 = C sh t - b1
     g*(chm(2)-T(2))-C*shm(2)+B*U(2),          # U2: a2 - x2 = 0
     g*(chm(3)-T(3))-C*shm(3)+B*U(3)]          # U3
Mx=sp.Matrix([[sp.diff(e,v) for v in (g,C,B)] for e in eqs])
det=sp.factor(sp.simplify(Mx.det()))
print('det =',det)
target=(r**2-1)*(r**2-2*c*r+1)**2*(c*(r**2+1)+r)/r**4
print('det identity:',sp.simplify(det-target)==0)
cw=-r/(1+r**2)
sol=sp.solve([e.subs(c,cw) for e in eqs[:2]],[C,B],dict=True)[0]
print('C =',sp.factor(sol[C]),' B =',sp.factor(sol[B]))
print('C==g th t:',sp.simplify(sol[C]-g*(r**2-1)/(r**2+1))==0,
      ' B formula:',sp.simplify(sol[B]-g*(r**4+r**2+1)/(r*(r**2+1)))==0,
      ' U3 residual:',sp.simplify(eqs[2].subs(c,cw).subs(sol))==0)
# --- (2) HB branch: Im s = sh v sin w => B = 2 sh v sin^2 w ; derive g = sh v / ch t and sh v = 1 - 1/(2 ch t)
cc=(r+1/r)/2; sin2=1-cw**2
shv=sp.symbols('shv')
gsol=sp.solve(sp.Eq(2*shv*sin2, g*(r**4+r**2+1)/(r*(r**2+1))),g)[0]
print('g = sh v / ch t :',sp.simplify(gsol-shv/cc)==0)
quad=sp.expand(sin2-shv**2-shv/cc)   # g = sin^2 w - sh^2 v
print('roots sh v:',[sp.simplify(x) for x in sp.solve(quad,shv)])
F=r**2*(r**2+1)-2*q*r*(r**2-r+1)-q**2*(r**2+1)
H=(r-q)**2*(r**2+1)-2*q*(q-(1-q)*r**2)
print('H==F:',sp.expand(H-F)==0)
shv_r=(r/q-q/r)/2
print('2q rho (1+rho^2)(sh v -1 + 1/(2ch t)) == F :',sp.simplify(2*q*r*(1+r**2)*(shv_r-1+1/(2*cc))-F)==0)
qr=r*(sp.sqrt((r**2-r+1)**2+(r**2+1)**2)-(r**2-r+1))/(r**2+1)
print('F(rho,q(rho))==0:',sp.simplify(sp.expand(F.subs(q,qr)))==0)
print('q(1)=',sp.nsimplify(sp.simplify(qr.subs(r,1))), ' kappa_c=',sp.simplify(((1+qr.subs(r,1))/(1-qr.subs(r,1)))**2))
print('g formula (rho^2-q^2)/(q(1+rho^2)) on F=0: ', sp.simplify(sp.numer(sp.together(shv_r/cc-(r**2-q**2)/(q*(1+r**2))))))
# --- (3) E_n lemma, symbolic in c for m=1..12  (t>0 param by rho; w with cos w=-1/(2c))
ok=True
for m in range(1,13):
    am=g*(chm(m)-T(m)); bm=B*U(m); xm=C*shm(m)-bm
    sub={C:g*(r**2-1)/(r**2+1),B:g*(r**4+r**2+1)/(r*(r**2+1)),c:cw}
    E=lambda n:chm(n)+2*cc**2*sp.chebyshevt(n,cw) if n>0 else 1+2*cc**2
    d1=sp.simplify(((am+xm)/g).subs(sub)-E(m+1)/cc)
    d2=sp.simplify(((am-xm)/g).subs(sub)-E(m-1)/cc)
    ok=ok and d1==0 and d2==0
print('E_n lemma m=1..12:',ok, '  E1,E2:',sp.simplify(E(1)),sp.simplify(E(2)))
assert ok

