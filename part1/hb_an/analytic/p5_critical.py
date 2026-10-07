"""The central family at the critical point: Theorem 50, Proposition 51, Proposition 52, and the limits in
the proof of Theorem 54, exact.

Every decision is a polynomial identity over Q, an identity in Q(y0) with q5(y0) = 0, or a
Taylor-dominance sign certificate on [0.23188, 0.23189]."""
import sys, os
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE), str(HERE.parent/'regimes')]
import sympy as s
from fractions import Fraction as Fr
from cert import certify_sign
from branch_systems import central, h, A,B,y,P,D,w2,w3, G,u,c
def check(name,cond):
    print(('OK   ' if cond else 'FAIL ')+name, flush=True)
    if not cond: sys.exit(1)
hh,vs,F,Q,su=central()
F0=[s.factor(f.subs(h,0)) for f in F]
q5=32*y**5-40*y**4+16*y**3-6*y**2-12*y+3; Q5=s.Poly(q5,y)
def red(e):
    n,d=s.fraction(s.together(e)); n=s.Poly(s.expand(n),y); d=s.Poly(s.expand(d),y)
    return s.rem(n*s.invert(d,Q5),Q5).as_expr()
lo,hi=Fr('0.23188'),Fr('0.23189')
def sgn(e):
    r=certify_sign(s.expand(red(e)),lo,hi); assert r is not None, 'uncertified sign'; return r[0]
prop=lambda a,b: s.simplify(s.cancel(a/b)).is_number and s.cancel(a/b)!=0
# --- primal equations at h=0 (as displayed in the paper, up to nonzero constants)
disp=[-12*A+4*D+P*(4*y**3-3*y-1)+2,
      -6*A+6*B-4*D+P*(4*y**4-4*y**3-7*y**2+4*y+3)+4,
      -12*B-12*D+P*(16*y**5-16*y**4-24*y**3+16*y**2+9*y-1)+18,
      -12*A+24*B-20*D+50+P*(64*y**7-64*y**6-128*y**5+112*y**4+80*y**3-56*y**2-15*y+7)]
check('primal equations at h=0 equal the displayed ones (up to constants)', all(prop(F0[i],disp[i]) for i in range(4)))
dual=F0[4:]
check('dual equations contain no A,B,D,P', all(not (e.free_symbols & {A,B,D,P}) for e in dual))
# normalization used in the paper: d1=-4*F0[4]/(2y+1), d2=F0[5], d3=-4*F0[6]
nd=[s.cancel(-4*dual[0]/(2*y+1)), dual[1], s.expand(-4*dual[2])]
lin=[s.Poly(s.expand(e),w2,w3) for e in nd]
check('dual equations linear in (w2,w3)', all(p.total_degree()<=1 for p in lin))
check('one dual equation is -3w2+22w3=12', any(prop(s.expand(s.numer(s.together(e))),-3*w2+22*w3-12) for e in dual))
Mx=s.Matrix([[p.coeff_monomial(w2),p.coeff_monomial(w3),p.coeff_monomial(1)] for p in lin])
det=s.factor(Mx.det())
check('consistency determinant = 900(2y+1)^3 q5(y)', s.expand(det-900*(2*y+1)**3*q5)==0 or s.expand(det+900*(2*y+1)**3*q5)==0)
print('    determinant =',det)
# --- solve
sol=s.solve([s.numer(s.together(e)) for e in dual[:2]],[w2,w3],dict=True)
sol=[t for t in sol if red(s.numer(s.together(dual[2])).subs(t))==0][0]
W2,W3=red(sol[w2]),red(sol[w3]); W5=red(1-W2-W3)
psol=s.solve(F0[:4],[A,B,D,P],dict=True)[0]
V={k:red(v) for k,v in psol.items()}; V[w2]=W2; V[w3]=W3
R=s.Rational
paper={A:-R(34,151)*y**4-R(77,302)*y**3+R(27,604)*y**2-R(75,302)*y+R(693,1208),
 B:-R(2275,151)*y**4+R(16769,604)*y**3-R(29287,1208)*y**2+R(19977,1208)*y-R(7471,2416),
 D:-R(243,151)*y**4+R(2097,604)*y**3-R(1959,1208)*y**2+R(1773,1208)*y+R(4381,2416),
 P:-R(396330,1057)*y**4+R(1393863,2114)*y**3-R(2294121,4228)*y**2+R(1529793,4228)*y-R(501309,8456),
 w2:(36571392*y**4+15609088*y**3+5194992*y**2+22511544*y+10388100)/R(43243231),
 w3:(4987008*y**4+2128512*y**3+708408*y**2+3069756*y+25003776)/R(43243231)}
check('closed forms of A,B,D,P,w2,w3 as in Theorem (critical point)', all(s.expand(V[k]-paper[k])==0 for k in paper))
check('w5 closed form', s.expand(W5-(7851355-25581300*y-5903400*y**2-17737600*y**3-41558400*y**4)/R(43243231))==0)
check('w2=2/5-22/25 w5, w3=3/5-3/25 w5', s.expand(red(W2-R(2,5)+R(22,25)*W5))==0 and s.expand(red(W3-R(3,5)+R(3,25)*W5))==0)
qv=lambda t: q5.subs(y,R(t))
check('q5(0.23188)>0>q5(0.23189), q5\'<0 on [0.2318,0.232]', qv('0.23188')>0>qv('0.23189') and certify_sign(s.diff(q5,y),Fr('0.2318'),Fr('0.232'))[0]==-1)
check('A0,D0,P0,w2,w3,w5>0 and B0<0', [sgn(V[A]),sgn(V[D]),sgn(V[P]),sgn(W2),sgn(W3),sgn(W5),sgn(V[B])]==[1,1,1,1,1,1,-1])
# --- Jacobian
J=s.Matrix(F).jacobian([A,B,y,P,D,w2,w3]).subs(h,0)
sub={A:V[A],B:V[B],P:V[P],D:V[D],w2:W2,w3:W3}
check('Jacobian block triangular (dual rows do not depend on A,B,P,D)', all(red(e)==0 for e in J.extract([4,5,6],[0,1,3,4]).subs(sub)))
dd=red(J.extract([4,5,6],[2,5,6]).subs(sub).det()); dp=red(J.extract([0,1,2,3],[0,1,4,3]).subs(sub).det())
check('both block determinants nonzero in Q(y0)', sgn(dd)!=0 and sgn(dp)!=0)
# --- dual quotient
Qc=s.fraction(s.cancel(Q.subs({G:1,u:R(1,2)})))[0]
Qc=s.Poly(s.expand(Qc.subs({w2:W2,w3:W3})),c)
quo,rem=s.div(Qc,s.Poly((c+R(1,2))**2*(c-y)**2,c))
check('Q(c) divisible by (c+1/2)^2(c-y0)^2', all(red(r)==0 for r in rem.all_coeffs()))
q2,q1,q0=[red(e) for e in quo.all_coeffs()]
check('quotient: q2>0, q1<0, q0-|q1|-|q2|>0', sgn(q2)==1 and sgn(q1)==-1 and sgn(q0+q1-q2)==1)
# normalization: at h=0 the polynomial Q of the system is 15 times -Phi/(G-c) of eq:dualpoly
wd={2:w2,3:w3,5:1-w2-w3}; uh=R(1,2); gh=1-c**2-uh**2
S1=1+sum(wd[m]*s.cancel((1-s.chebyshevt(m,c))/(m*(1-c))) for m in wd)
S2=sum(wd[m]*(m-s.chebyshevu(m-1,c))/m for m in wd)
Qpaper=-(gh*S1-2*uh*(1+c)*S2)
check('Q at h=0 equals 15 times the dual polynomial of eq:dualpoly', s.expand(s.cancel(Q.subs({G:1,u:uh}))-15*Qpaper)==0)
y0n=s.nsolve(q5,y,0.231889,prec=30)
qn=[s.N(e.subs(y,y0n)/15,20) for e in (q2,q1,q0)]
print('normalized quotient coefficients q2, q1, q0 =', [s.N(v,8) for v in qn])
check('normalized quotient: q2=0.09328..., q1=-0.05001..., q0=0.65052... (Proposition p5h0dual)',
      all(abs(v-t)<1e-5 for v,t in zip(qn,(0.093280,-0.050019,0.650521))))
# --- first-order slacks
T=s.chebyshevt; U=s.chebyshevu
Y=lambda m,sg: (R(3,4)-y**2)*(1-T(m,y))+sg*(1-y**2)*U(m-1,y)
def divslack(m,sg):
    r=m%3
    if r==0: return R(m*m,2)-sg*m*(V[B]+V[D])+V[P]*Y(m,sg)
    if r==1 and sg==-1: return (V[B]-3*V[A])+R(m*m,2)-m*V[B]+m*V[D]+V[P]*Y(m,-1)
    if r==2 and sg==1: return (V[B]-3*V[A])+R(m*m,2)+m*V[B]-m*V[D]+V[P]*Y(m,1)
active={(1,-1),(2,1),(3,1),(5,1)}
for m in range(1,8):
    for sg in (-1,1):
        e=divslack(m,sg)
        if e is None: continue
        r=red(e)
        if (m,sg) in active: check('first-order slack (%d,%+d) vanishes (active)'%(m,sg), r==0)
        else: check('first-order slack (%d,%+d) > 0'%(m,sg), sgn(r)==1)
check('m>=8 constants: |B0|+D0<2.37, |B0-3A0|<1.81, P0<2.58', sgn(R(237,100)+V[B]-V[D])==1 and sgn(R(181,100)-3*V[A]+V[B])==1 and sgn(R(258,100)-V[P])==1)
check('m=8: 32-2.37*8-1.81-6.45>0', 32-R(237,100)*8-R(181,100)-R(645,100)>0)
check('limits: |B0-A0|<1, |B0-A0+P0(1/4-y0^2)|<2, 2|B0|/sqrt3<1/3 (4B0^2<1/3)',
      sgn(1-V[A]+V[B])==1 and sgn(2-abs(0)+V[B]-V[A]+V[P]*(R(1,4)-y**2))==1 and sgn(R(1,3)-4*V[B]**2)==1)
check('limit |sqrt3 A0+B0/sqrt3|<2  (3A0+B0<2 sqrt3, via (3A0+B0)^2<12)', sgn(12-(3*V[A]+V[B])**2)==1)
check('strict limits of (centralbounds): A0<0.55, B0>-0.27, P0<3, D0<2.2', sgn(R(55,100)-V[A])==1 and sgn(V[B]+R(27,100))==1 and sgn(3-V[P])==1 and sgn(R(22,10)-V[D])==1)
print('ALL CHECKS PASSED')
