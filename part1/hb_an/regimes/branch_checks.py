"""Per-slab interval checks of Proposition 43 (P1), Proposition 55 (P5) and the closed-form arc P2.

P1: B < 0, g/z > 1/2, 0 < u < 1, u_z < 0, (1/2 + zB)' < 0, the factorization of Q1 with q0 - |q1| > 0,
and u_z^2 z^2 (1 - z^2) < 1 + u^2. P5: the bounds of equation centralbounds and of the proof of
Proposition 55, the dual quotient q0 - |q1| - |q2| > 0, and the raw slacks of all lags m < 24, evaluated
after division by h near h = 0. P2: an independent interval check of the closed form of Theorem 49."""
import sympy as s
from flint import arb
from interval_core import funcs
from branch_systems import G,u,c,z,h,A,B,y,P,D,w2,w3,L,U

def check_p1(Q):
    xv=z*z*B;uv=1-z/2+z*z*A
    co=s.Poly(Q,c).all_coeffs()
    b2=co[0];b1=co[1]+2*xv*b2;b0=co[2]+2*xv*co[1]+3*xv*xv*b2
    raw=[-B,(s.expand(1-xv*xv-uv*uv)/z).cancel(),b0.subs(u,uv),b1.subs(u,uv),s.Rational(1,2)-z*A]
    fun=funcs(raw,(z,[A,B]))
    def check(zz,X,der):
        bneg,gz,q0,q1,positive=fun(zz,X)
        assert bneg>0 and gz>arb('0.5') and q0-abs(q1)>0 and positive>0
        # c >= -z/2 follows from the exact endpoint and this derivative.
        assert X[1]+zz*der[1]<0
        uz=-arb('0.5')+2*zz*X[0]+zz*zz*der[0]
        assert uz<0
        uv=1-zz/2+zz*zz*X[0]
        assert uz*uz*zz*zz*(1-zz*zz)<1+uv*uv
    return check

def check_p2(vs):
    uv,xv,yv,pv,dv,W=vs
    subs={u:uv}
    def mean(e):return (1-pv)*e.subs(c,xv)+pv*e.subs(c,yv)
    e5=mean(U(5))-dv
    df=funcs([s.diff(e5,G)]+[s.diff(e5,v) for v in vs],(G,vs))
    def check(gv,X,der):
        uv,xv,yv,pv,dv,W=X
        assert uv>arb('.50') and uv<arb('.52') and xv>-1 and xv<0 and yv>0 and yv<1
        assert W>0 and W<1 and pv<1
        assert der[0]>0 and der[0]<1
        assert der[3]<0
        ds=df(gv,X);assert ds[0]+sum(ds[i+1]*der[i] for i in range(6))>0
        tt=gv.acosh();cc=[xv,yv];ww=[1-pv,pv]
        gg=[1-ci*ci-uv*uv for ci in cc]
        gbar=sum(ww[i]*gg[i] for i in range(2));eta=gbar-abs(dv)/(gv*gv-1).sqrt()
        assert eta>arb('.24')
        angles=[ci.acos() for ci in cc]
        for m in range(1,16):
            aa=sum(ww[i]*gg[i]*((m*tt).cosh()-(m*angles[i]).cos()) for i in range(2))
            bb=sum(ww[i]*2*uv*(1-cc[i]*cc[i]).sqrt()*(m*angles[i]).sin() for i in range(2))
            rr=tt.sinh()/(m*tt).sinh()
            if m!=1:assert (aa-bb)*rr+dv>0,('L',m,str(gv))
            if m not in (2,3,5):assert (aa+bb)*rr-dv>0,('U',m,str(gv))
    return check

class Jet:
    """f(h)=a+h*b; the divided part is enclosed without dividing by h."""
    def __init__(self,a,b=0):self.a=arb(a);self.b=arb(b)
    def __add__(self,o):
        if not isinstance(o,Jet):o=Jet(o)
        return Jet(self.a+o.a,self.b+o.b)
    __radd__=__add__
    def __neg__(self):return Jet(-self.a,-self.b)
    def __sub__(self,o):return self+-o if isinstance(o,Jet) else self+(-o)
    def __rsub__(self,o):return -self+o
    def __mul__(self,o):
        if not isinstance(o,Jet):o=Jet(o)
        return Jet(self.a*o.a,self.a*o.b+self.b*o.a+Jet.h*self.b*o.b)
    __rmul__=__mul__

def check_central(Q):
    xv=-s.Rational(1,2)+h*B;uv=s.Rational(1,2)+h*A
    co=s.Poly(Q,c).all_coeffs()
    # Quotient by [(c-x)(c-y)]^2; only its leading three terms are needed.
    a=-(xv+y);b=xv*y
    q2=co[0];q1=co[1]-2*a*q2;q0=co[2]-2*a*q1-(a*a+2*b)*q2
    fq=funcs([e.subs({G:1+h,u:uv}) for e in (q2,q1,q0)],(h,[A,B,y,P,D,w2,w3]))
    # Shifted Chebyshev polynomials have nonnegative coefficients.
    tt={m:[int(q) for q in s.Poly(s.cancel((s.chebyshevt(m,1+h)-1)/h),h).all_coeffs()] for m in range(1,24)}
    uu={m:[int(q) for q in s.Poly(s.chebyshevu(m-1,1+h),h).all_coeffs()] for m in range(1,24)}
    def horner(co,x):
        v=arb(0)
        for a in co:v=v*x+a
        return v
    def check(hh,X,der):
        aa,bb,yy,pp,dd,ww2,ww3=X
        uv=arb('.5')+hh*aa;xx=-arb('.5')+hh*bb;pv=hh*pp;dv=hh*dd
        assert aa>0 and aa<arb('.55') and bb>-arb('0.27') and bb<0
        assert pp>0 and pp<3 and dd>0 and dd<arb('2.2')
        assert yy>arb('.23') and yy<arb('.29') and ww2>arb('.36') and ww2<arb('.38')
        assert ww3>arb('.59') and ww3<arb('.65') and der[5]+der[6]>0
        assert 1-ww2-ww3<arb('.031')
        q2,q1,q0=fq(hh,X);assert q0-abs(q1)-abs(q2)>0
        g1=1-xx*xx-uv*uv;g2=1-yy*yy-uv*uv;gb=(1-pv)*g1+pv*g2
        assert gb>arb('.49')
        assert gb-dd*(hh/(2+hh)).nonnegative_part().sqrt()>arb('.25')
        assert abs(bb)/(1-xx*xx).sqrt()<arb(1)/3
        assert abs(bb-aa-hh*(bb*bb+aa*aa))<1
        assert abs(bb-aa-hh*(bb*bb+aa*aa)+pp*(g2-g1))<2
        assert abs(aa+hh*der[0])<1
        # Main-atom phase difference: omega-2pi/3=h*delta_h.
        omega_sin=arb(3).sqrt()/2
        hull=arb(-arb('.5')+hh*bb/2,abs(hh*bb)/2)
        delta_h=-bb/(1-hull*hull).sqrt();delta=hh*delta_h
        Jet.h=hh
        uj=Jet(arb('.5'),aa);cj=Jet(-arb('.5'),bb)
        gj=1-cj*cj-uj*uj
        sj=Jet(omega_sin,(bb-hh*bb*bb)/((1-xx*xx).sqrt()+omega_sin))
        py=Jet(0,pp);dj=Jet(0,dd)
        a2=yy.acos();s2=(1-yy*yy).sqrt()
        for m in range(1,24):
            cm=arb(1) if m%3==0 else -arb('.5')
            sm=arb(0) if m%3==0 else (omega_sin if m%3==1 else -omega_sin)
            md=m*delta;mdh=m*delta_h
            cdiff=-mdh*mdh*hh/2*(md/2).sinc()**2
            sdiff=mdh*md.sinc()
            cmj=Jet(cm,cm*cdiff-sm*sdiff);smj=Jet(sm,sm*cdiff+cm*sdiff)
            chj=Jet(1,horner(tt[m],hh));um=horner(uu[m],hh)
            aj=(1-py)*gj*(chj-cmj)+py*g2*(chj-(m*a2).cos())
            bj=(1-py)*2*uj*sj*smj+py*2*uj*s2*(m*a2).sin()
            for sign in (-1,1):
                if sign==-1 and m==1 or sign==1 and m in (2,3,5):continue
                sl=aj+sign*bj-sign*dj*um
                # At h = 0 every raw slack is 0 (critical pairs) or 3/2 (noncritical pairs).
                critzero=(m%3==0 or (sign==-1 and m%3==1) or (sign==1 and m%3==2))
                lower=sl.b if critzero else sl.b+(arb(3)/2)/arb(hh.upper())
                assert lower>0,('central slack',sign,m,str(hh),str(lower))
    return check
