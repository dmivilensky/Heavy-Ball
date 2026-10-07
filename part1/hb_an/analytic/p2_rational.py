"""The arc P2 in closed form: Lemmas 44, 46, 47, 48, Definition 45, Theorem 49, Tables 3 and 4, exact.

Every decision is a polynomial identity over Q or a comparison of rational numbers. Sign certificates
use Taylor dominance at the single center zc = 1.69725 with radius 0.00175 (Lemma 48)."""
import sympy as s, sys
from fractions import Fraction as Fr
z,c,Gs,us,S=s.symbols('zeta c G u s')
T=s.chebyshevt; U=s.chebyshevu
def R(x): f=Fr(x); return s.Rational(f.numerator,f.denominator)
ok=lambda e: s.cancel(s.together(e))==0
def check(name,cond):
    print(('OK   ' if cond else 'FAIL ')+name, flush=True)
    if not cond: sys.exit(1)
qm,qp=z**2-2*z+2,z**2+2*z+2; rm,rp=z**2-2,z**2+2; em,ep=z**2-2*z-2,z**2+2*z-2
K6=z**6-4*z**5+4*z**4+8*z**2-16*z+8
G=(z**4+4)/(8*z*rm); u=qm**2/(z**4-4)
W=qm**2*qp*(7*z**2-18*z+14)/(8*rp*K6)
a=qm**3/(2*z*rm*rp); b=-2*z*(3*z**4-10*z**3+16*z**2-20*z+12)/(rm**2*rp)
M10=5*z**10-44*z**9+170*z**8-368*z**7+744*z**6-1376*z**5+1488*z**4-1472*z**3+1360*z**2-704*z+160
D14=5*z**14-34*z**13+42*z**12+648*z**11-3460*z**10+8296*z**9-14984*z**8+23616*z**7-29968*z**6+33184*z**5-27680*z**4+10368*z**3+1344*z**2-2176*z+640
mc=4*z**2*(15*z**8-92*z**7+224*z**6-440*z**5+760*z**4-880*z**3+896*z**2-736*z+240)/(rm*M10)
d=qm**3*D14/(rm**3*rp**2*M10)
N10=5*z**10-24*z**9+74*z**8-512*z**7+1448*z**6-2240*z**5+2896*z**4-2048*z**3+592*z**2-384*z+160
N8=5*z**8-24*z**7+48*z**6-104*z**5+184*z**4-208*z**3+192*z**2-192*z+80
N12=25*z**12-314*z**11+1764*z**10-6092*z**9+14988*z**8-28736*z**7+44832*z**6-57472*z**5+59952*z**4-48736*z**3+28224*z**2-10048*z+1600
f1=z**3-z**2-2
Ha=162529280000*Gs**12-238685519872*Gs**10+77668708352*Gs**8-6246545408*Gs**6-58821872*Gs**4+12933112*Gs**2+229441
P2=256*Gs**4*us**2-256*Gs**4-512*Gs**3*us**3+768*Gs**3*us-32*Gs**2*us**4-480*Gs**2*us**2-16*Gs*us**3-us**4-16*us**2

print('--- identities')
check('Lemma curve: P2(G(z),u(z))=0', ok(P2.subs({Gs:G,us:u})))
al=32*S**2*(8-16*S-S**2); be=-(S**4+16*S**3+480*S**2-768*S+256); ga=-16*S**2
check('P2(G,sG)=G^2(alpha G^4+beta G^2+gamma)', ok(P2.subs(us,S*Gs)-Gs**2*(al*Gs**4+be*Gs**2+ga)))
check('discriminant factorization', ok(be**2-4*al*ga-(S-4)**2*(S**2+32*S-16)**2*(S**2-24*S+16)))
check('u/G formula', ok(u/G-8*z*qm/(rp*qp)))
check('K6 = z^4(z-2)^2+8(z-1)^2', ok(K6-(z**4*(z-2)**2+8*(z-1)**2)))
check('1-W = e-^2 e+^2/(8 r+ K6)', ok(1-W-em**2*ep**2/(8*rp*K6)))
k=1-c*c
al_=1+Gs and None
A2=(Gs+c)/Gs; B2=1/Gs; A3=(4*Gs**2+4*Gs*c+4*c**2-3)/(4*Gs**2-1); B3=4*(Gs+c)/(4*Gs**2-1)
Wg=s.symbols('W')
Mg=8*Gs**3+4*Gs**2*c+4*Gs*c**2-4*Gs-Wg*(4*Gs*c**2-2*Gs+c); Ng=4*Gs**2+4*Gs*c-Wg*(4*Gs*c+1)
check('M = G(4G^2-1) alpha, N = G(4G^2-1) beta', ok(Mg-Gs*(4*Gs**2-1)*(1+Wg*A2+(1-Wg)*A3)) and ok(Ng-Gs*(4*Gs**2-1)*(Wg*B2+(1-Wg)*B3)))
M=Mg.subs({Gs:G,Wg:W}); N=Ng.subs({Gs:G,Wg:W})
check('Lemma p2ids(b): -{(k-u^2)M-2ukN} = 4G(1-W)(c^2+ac+b)^2', ok(-((k-u*u)*M-2*u*k*N)-4*G*(1-W)*(c*c+a*c+b)**2))
def integ(poly_c):
    p=s.Poly(s.expand(poly_c),c); r0,r1=s.Integer(0),s.Integer(0)
    for co in p.all_coeffs(): r0,r1=co-r1*b, r0-r1*a
    return r0+r1*mc
def raw(m,sg): return (1-c*c-u*u)*(T(m,G)-T(m,c))+sg*2*u*(1-c*c)*U(m-1,c)
def slack(m,sg): return integ(raw(m,sg))-sg*d*U(m-1,G)
check('Lemma p2ids(c): l1=u2=u3=0', ok(slack(1,-1)) and ok(slack(2,1)) and ok(slack(3,1)))
check('Lemma p2ends(a): 32G^3-36G^2+4G-1 factorization', ok(32*G**3-36*G**2+4*G-1-(z**3+2*z+4)*f1*(z**6-8*z**5-2*z**4+16*z**3+4*z**2-32*z-8)/(16*z**3*rm**3)))
check('Lemma p2ends(a): u-1+1/(2G) = 8 f1^2/(r- r+ q- q+)', ok(u-1+1/(2*G)-8*f1**2/(rm*rp*qm*qp)))
num=s.numer(s.together(Ha.subs(Gs,G)))
check('Lemma p2ends(b): N12 divides numerator of Ha(G(z))', s.rem(s.Poly(s.expand(num),z),s.Poly(N12,z)).is_zero)
Glo,Ghi=G.subs(z,R('1.699')),G.subs(z,R('1.6989'))
check('Lemma p2ends(b): Ha has exactly one root in [G(1.699),G(1.6989)], an interval containing (1.02358,1.02359)',
      s.Poly(Ha,Gs).count_roots(Glo,Ghi)==1 and Glo<R('1.02358') and R('1.02359')<Ghi)
check('dG/dz = r+ e- e+/(8 z^2 r-^2)', ok(s.diff(G,z)-rp*em*ep/(8*z**2*rm**2)))
check('1+a+b = q-^3 e+/(2z r-^2 r+)', ok(1+a+b-qm**3*ep/(2*z*rm**2*rp)))
check('1-a+b = -q-^3 e-/(2z r-^2 r+)', ok(1-a+b+qm**3*em/(2*z*rm**2*rp)))
check('2-a formula', ok(2-a+(z**6-10*z**5+18*z**4-32*z**3+36*z**2-8*z+8)/(2*z*rm*rp)))
check('2+a formula', ok(2+a-(z**6-2*z**5+18*z**4-32*z**3+36*z**2-40*z+8)/(2*z*rm*rp)))
check('q(mc) formula', ok(mc**2+a*mc+b-32*z**3*em*ep*(z**3+2*z-4)*f1*N10/(rm**2*rp*M10**2)))
gbar=1-u**2+a*mc+b
check('gbar = int g dnu = 1-u^2+a mc+b', ok(integ(1-c*c-u*u)-gbar))
check('gbar = 8 z q-^4 N8/(r-^2 r+^2 M10)', ok(gbar-8*z*qm**4*N8/(rm**2*rp**2*M10)))
uG=s.diff(u,z)/s.diff(G,z)
check('du/dG formula', ok(uG-32*z**2*qm*(z**4-2*z**3-4*z+4)/(rp**3*em*ep)))
# Table p2lags
rows={(1,1):(2*qm**3*(z**4-2*z**3-4*z+4), rp**2*rm**3),
      (2,-1):((z**2-4*z+2)*qp*qm**4*(z**4+2*z**3-8*z**2+4*z+4), 2*z*rp**3*rm**4),
      (3,-1):(em*ep*qm**3, 8*z**2*rp**4*rm**5),(4,-1):(em*ep*qm**4, 32*z**3*rp**5*rm**6),
      (4,1):(-qm**4, 2*z*rp**2*rm**6),(5,-1):(qm**5, 128*z**4*rp**6*rm**7),
      (5,1):(em*ep*qm**3, 8*z**2*rp**2*rm**7),(6,-1):(em*ep*qm**4, 512*z**5*rp**7*rm**8),
      (6,1):(-em*ep*qm**4, 32*z**3*rp**2*rm**8),(7,-1):(em*ep*qm**3, 2048*z**6*rp**8*rm**9),
      (7,1):(qm**5, 128*z**4*rp**2*rm**9),(8,-1):(qm**4, 8192*z**7*rp**9*rm**10),
      (8,1):(-em*ep*qm**4, 512*z**5*rp**3*rm**10),(9,-1):(em*ep*qm**5, 32768*z**8*rp**10*rm**11),
      (9,1):(em*ep*qm**3, 2048*z**6*rp**4*rm**11)}
resid={}
for (m,sg),(nf,df) in rows.items():
    q_=s.cancel(s.together(slack(m,sg)*M10*df/(N10*nf)))
    isp=s.denom(q_)==1 or s.Poly(s.denom(q_),z).degree()==0
    check('Table p2lags row %s_%d: residual is a polynomial'%('l' if sg<0 else 'v',m), isp)
    resid[(m,sg)]=s.expand(q_)
check('row v5: residual = N12', ok(resid[(5,1)]-N12))
check('rows l1? (v1,l2 have no residual): residual constants', all(s.Poly(resid[r],z).degree()==0 for r in [(1,1),(2,-1)]))

print('--- signs on J=[1.6955,1.699] by single-centre Taylor dominance')
zc=Fr(169725,100000); dl=Fr(175,100000)
def dom(p):
    e=s.symbols('e'); co=s.Poly(s.expand(p.subs(z,R(zc)+e)),e).all_coeffs()[::-1]
    co=[Fr(int(s.numer(x)),int(s.denom(x))) for x in co]
    tail=sum(abs(x)*dl**i for i,x in enumerate(co) if i>0)
    return co[0],tail
cert={}
def sgn_poly(p,label=''):
    p=s.expand(p)
    key=s.srepr(p)
    if key not in cert:
        v,t=dom(p)
        check('dominance %-34s deg %2d ratio %.4f'%(label[:34],s.Poly(p,z).degree(),float(t/abs(v))), abs(v)>t)
        cert[key]=(1 if v>0 else -1, p, label, v, t)
    return cert[key][0]
def sign_of(expr,label):
    n,dn=s.fraction(s.factor(s.together(expr)))
    sg=1
    for part in (n,dn):
        cst,fl=s.factor_list(part); sg*=(1 if cst>0 else -1)
        for f,e in fl:
            if s.Poly(f,z).degree()==0: continue
            sg*=sgn_poly(f,label+' factor')**e
    return sg
# table polynomials explicitly
for nm,p,exp_ in [('z',z,1),('q-',qm,1),('q+',qp,1),('r-',rm,1),('r+',rp,1),('e+',ep,1),('7z^2-18z+14',7*z**2-18*z+14,1),('K6',K6,1),
                  ('e-',em,-1),('z^2-4z+2',z**2-4*z+2,-1),('3z^4-..',3*z**4-10*z**3+16*z**2-20*z+12,1),
                  ('sextic(2-a)',z**6-10*z**5+18*z**4-32*z**3+36*z**2-8*z+8,-1),('sextic(2+a)',z**6-2*z**5+18*z**4-32*z**3+36*z**2-40*z+8,1),
                  ('z^4-2z^3-4z+4',z**4-2*z**3-4*z+4,-1),('z^4+2z^3-8z^2+4z+4',z**4+2*z**3-8*z**2+4*z+4,1),('z^3+2z-4',z**3+2*z-4,1),
                  ('N10',N10,1),('M10',M10,-1),('N8',N8,-1),('D14',D14,-1),("N12'",s.diff(N12,z),-1)]:
    check('table sign %s'%nm, sgn_poly(p,nm)==exp_)
for (m,sg),r in resid.items():
    if (m,sg) in [(1,1),(2,-1),(5,1)]: continue
    sgn_poly(r,'R_%d%s'%(m,'-' if sg<0 else '+'))
print('--- root isolation')
v=lambda f,x: f.subs(z,R(x))
check('f1(1.6956)<0<f1(1.6957), f1 increasing for z>2/3', v(f1,'1.6956')<0<v(f1,'1.6957'))
check('N12(1.6989)>0>N12(1.699)', v(N12,'1.6989')>0>v(N12,'1.699'))
check('J contains [zeta_1,zeta_a]', R('1.6955')<R('1.6956') and R('1.699')<=R('1.699'))
print('--- Theorem p2')
check('(i) G decreasing', sign_of(s.diff(G,z),'dG')==-1)
check('(i) G>1.0233', sign_of(G-R('1.0233'),'G-1.0233')==1)
check('(i) G<1.0337', sign_of(R('1.0337')-G,'1.0337-G')==1)
check('(i) u>0.511', sign_of(u-R('0.511'),'u-0.511')==1)
check('(i) u<0.517', sign_of(R('0.517')-u,'0.517-u')==1)
check('(ii) b<0', sign_of(b,'b')==-1)
check('(ii) 1+a+b>0, 1-a+b>0', sign_of(1+a+b,'1+a+b')==1 and sign_of(1-a+b,'1-a+b')==1)
check('(ii) |a|<2', sign_of(2-a,'2-a')==1 and sign_of(2+a,'2+a')==1)
q=lambda cc: cc*cc+a*cc+b
check('(ii) atom ranges', [sign_of(q(R(t)),'q(%s)'%t) for t in ('-0.505','-0.48','0.25','0.285')]==[1,-1,-1,1])
check('(iii) mc+a/2<0', sign_of(mc+a/2,'mc+a/2')==-1)
check('(iii) q(mc)/f1 < 0 on J', sign_of(s.cancel((mc**2+a*mc+b)/f1),'q(mc)/f1')==-1)
check('(iv) W>0.362', sign_of(W-R('0.362'),'W-0.362')==1)
check('(iv) W<0.366', sign_of(R('0.366')-W,'0.366-W')==1)
check('(v) d>0', sign_of(d,'d')==1)
check('(v) gbar>0.28', sign_of(gbar-R('0.28'),'gbar-0.28')==1)
check('(v) (gbar-0.28)^2(G^2-1)>d^2', sign_of((gbar-R('0.28'))**2*(G**2-1)-d**2,'eta')==1)
check('(vi) gbar+2u^2<1.05', sign_of(R('1.05')-gbar-2*u*u,'1.05-gbar-2u^2')==1)
check('(vi) T_10(1.0233)>4.369', T(10,R('1.0233'))>R('4.369'))
check('(vi) 0.28*4.369-1.05>0.17', R('0.28')*R('4.369')-R('1.05')>R('0.17'))
for (m,sg) in rows:
    if (m,sg)==(5,1):
        check('(vi) v5 = N12*(positive)', sign_of(s.cancel(slack(5,1)/N12),'v5/N12')==1); continue
    check('(vi) slack %s_%d>0'%('l' if sg<0 else 'v',m), sign_of(slack(m,sg),'slack')==1)
check('(vii) 0<du/dG', sign_of(uG,'uG')==1)
check('(vii) du/dG<1', sign_of(1-uG,'1-uG')==1)
# write residual file
from pathlib import Path
with open(Path(__file__).resolve().parent/'p2_residuals.txt','w') as fh:
    fh.write('# All polynomial factors whose sign on J=[1.6955,1.699] is used in the proof of Theorem B.12.\n')
    fh.write('# Certificate: |R(zc)| > sum_{k>=1}|R_k| delta^k, zc=%s, delta=%s (Taylor coefficients R_k at zc).\n\n'%(zc,dl))
    for key,(sg_,p,label,v_,t_) in cert.items():
        fh.write('%s (degree %d), sign %s:\n%s\n  R(zc) = %s\n  sum_k |R_k| delta^k = %s\n  ratio = %.6f\n\n'%(label,s.Poly(p,z).degree(),'+' if sg_>0 else '-',p,v_,t_,float(t_/abs(v_))))
print('number of certified factors:',len(cert),' max ratio:',max(float(t/abs(v)) for (_,_,_,v,t) in cert.values()))
print('ALL CHECKS PASSED')
