"""Endpoint identities and branch identification (Theorem 12, Lemma 47, Proposition 55).

Exact: the endpoint z_W = sqrt(3) - 1 of P1 (spiral data), G_1 (root of 32G^3 - 36G^2 + 4G - 1) and G_a
(root of H_a) isolated to 10^-65; Arb: the exact endpoint data lie in the certified boxes of the adjacent
families, which identifies the branches at G_W, G_1 and G_a."""
import sympy as s,mpmath as mp,pickle,json
from pathlib import Path
from flint import arb,ctx
from interval_core import funcs,inside
from branch_systems import p1,p2,central,G,u,c,z
ctx.prec=256;mp.mp.dps=90
dat=pickle.loads(Path('p2_field.pkl').read_bytes())
Ha=162529280000*G**12-238685519872*G**10+77668708352*G**8-6246545408*G**6-58821872*G**4+12933112*G**2+229441
H1=32*G**3-36*G**2+4*G-1
def rational(q):return arb(int(s.numer(q)))/int(s.denom(q))
def isolate(poly,lo,hi):
    roots=s.Poly(poly,G).intervals(eps=s.Rational(1,10**65),inf=lo,sup=hi)
    assert len(roots)==1 and roots[0][1]==1
    a,b=map(rational,roots[0][0]);return arb((a+b)/2,(b-a)/2)
ga=isolate(Ha,s.Rational(102358,100000),s.Rational(102359,100000))
g1=isolate(H1,s.Rational(103329,100000),s.Rational(103330,100000))
fp=funcs(dat['P'],(G,u));fpu=funcs(s.diff(dat['P'],u),(G,u));fpm=s.lambdify((G,u),dat['P'],'mpmath')
def physical(g):
    gm=mp.mpf(g.mid().str(90,radius=False));um=mp.findroot(lambda vv:fpm(gm,vv),mp.mpf('.512'))
    uc=arb(str(um));uu=arb(uc,arb('1e-55'))
    inv=1/fpu(arb(g.mid()),uc)
    kk=uc-inv*fp(g,uc)+(1-inv*fpu(g,uu))*(uu-uc)
    assert inside(kk,uu)
    vals={k:funcs(v,(G,u))(g,uu) for k,v in dat.items() if k!='P'}
    disc=vals['a']**2-4*vals['b'];assert disc>0
    x=(-vals['a']-disc.sqrt())/2;y=(-vals['a']+disc.sqrt())/2
    p=(vals['mc']-x)/(y-x)
    return [uu,x,y,p,vals['d'],vals['W']],vals['s5']
def identify(name,par,point):
    boxes=json.loads(Path(name+'_boxes.json').read_text())
    for b in boxes:
        if par.lower()>=arb(str(b['a'])) and par.upper()<=arb(str(b['b'])):
            X=[arb(arb(str(x)),arb(str(r))) for x,r in zip(b['center'],b['radii'])]
            assert all(inside(x,xx) for x,xx in zip(point,X)),(name,[str(x) for x in point],X)
            print('identified',name,str(par),flush=True);return
    raise AssertionError(('no box',name,str(par)))
# P1--spiral endpoint, exactly.
par,vs,F,Q,su=p1();zw=s.sqrt(3)-1
assert all(s.simplify(f.subs({z:zw,vs[0]:0,vs[1]:-1/(2*zw)}))==0 for f in F)
zwi=arb(3).sqrt()-1
identify('p1',zwi,[arb(0),-1/(2*zwi)])
# P2--spiral endpoint, exactly modulo H1.
par,vs,F,Q,su=p2();us=1-1/(2*G);xs=-1/(2*G);ys=(2*G-1)/(4*G)
gs=(2*G-1)/(2*G*G);ds=gs*(G*G-1)/G;ws=2*(1+2*G-2*G*G)/(4*G+1)
exact=[us,xs,ys,0,ds,ws]
for f in F:
    num=s.cancel(f.subs(dict(zip(vs,exact)))).as_numer_denom()[0]
    assert s.rem(s.Poly(num,G),s.Poly(H1,G)).is_zero
point=[funcs(e,(G,))(g1) if not isinstance(e,int) else arb(e) for e in exact]
identify('p2',g1,point)
# U5 changes sign on the physical P2 branch. The exact resultant
# Lemma p2ends(b): N12 divides the numerator of Ha(G(zeta)) (checked in analytic/p2_rational.py).
lo=arb(1023581427007)/10**12;hi=arb(1023581427009)/10**12
assert physical(lo)[1]<0 and physical(hi)[1]>0
aa,sl5=physical(ga);identify('p2',ga,aa)
uv,xv,yv,pv,dv,W=aa;hh=ga-1
point5=[(uv-arb('.5'))/hh,(xv+arb('.5'))/hh,yv,pv/hh,dv/hh,W,1-W]
identify('central',hh,point5)
assert aa[3]>0 and aa[3]<1 and aa[5]>0 and aa[5]<1
print('ALL ENDPOINT IDENTITIES AND BRANCH IDENTIFICATIONS PASSED')
