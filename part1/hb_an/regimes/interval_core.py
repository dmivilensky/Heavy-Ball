"""Outward-rounded Arb enclosures and the parametric Krawczyk test (equation krawczyk, Appendix B).

For a slab I and a box X with center Z0: K = Z0 - R F(I, Z0) + (Id - R F_Z(I, X))(X - Z0) must lie in the
interior of X, and a weighted infinity norm of Id - R F_Z(I, X) must be smaller than one. F(I, Z0) is
enclosed by the mean-value form F(s0, Z0) + (I - s0) F_s(I, Z0). Derivative enclosures of the branch use
Z' = (Id - M)^{-1}(-R F_s), M = Id - R F_Z, with the Neumann remainder bound. Centers and the matrices R
are computed in floating point; every acceptance decision is an Arb interval comparison."""
import sympy as s
from sympy.printing.pycode import PythonCodePrinter
import numpy as np
from scipy.optimize import root
from flint import arb,ctx
ctx.prec=160
class Printer(PythonCodePrinter):
    def _print_Rational(self,expr):return '(arb(%s)/%s)'%(expr.p,expr.q)
def funcs(expr,args):
    return s.lambdify(args,expr,modules=[{'arb':arb}],printer=Printer(),cse=True)
def ball(lo,hi):
    a,b=arb(str(lo)),arb(str(hi))
    return arb((a+b)/2,(b-a)/2)
def center(x):return float(x.mid())
def inside(a,b):return a.lower()>b.lower() and a.upper()<b.upper()
def run(F,par,vs,grid,guess,check=None,name='branch',radius_factor=4):
    n=len(vs);jac=s.Matrix(F).jacobian(vs)
    ff=s.lambdify((par,vs),F,'numpy',cse=True)
    jf=s.lambdify((par,vs),jac,'numpy',cse=True)
    fi=funcs(F,(par,vs));ji=funcs(jac.tolist(),(par,vs))
    dp=funcs([s.diff(e,par) for e in F],(par,vs))
    accepted=[];prev=None
    for a,b in zip(grid[:-1],grid[1:]):
        mid=(a+b)/2
        sol=root(lambda x:ff(mid,x),guess,jac=lambda x:jf(mid,x),tol=1e-11)
        z=sol.x
        if np.max(np.abs(ff(mid,z)))>1e-7:raise ValueError(('Newton',a,b,sol.message))
        guess=z;A=np.linalg.inv(jf(mid,z));Af=[[arb(str(t)) for t in row] for row in A]
        pp=ball(a,b);zz=[arb(str(t)) for t in z]
        # Parameter derivatives are used only to propose the box radius.
        pval=np.array([float(s.N(s.diff(e,par).subs({par:mid,**dict(zip(vs,z))}))) for e in F]) if False else np.array([float(t.mid()) for t in dp(arb(str(mid)),zz)])
        slope=-A@pval
        rad=radius_factor*(b-a)/2*np.maximum(np.abs(slope),.05)+1e-10
        X=[arb(zz[i],arb(str(rad[i]))) for i in range(n)]
        J=ji(pp,X); f0=fi(arb(str(mid)),zz)
        fp=dp(pp,zz)
        B=[[arb(int(i==j))-sum(Af[i][k]*J[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
        K=[zz[i]-sum(Af[i][k]*f0[k] for k in range(n))-(pp-arb(str(mid)))*sum(Af[i][k]*fp[k] for k in range(n))+sum(B[i][j]*(X[j]-zz[j]) for j in range(n)) for i in range(n)]
        if not all(inside(K[i],X[i]) for i in range(n)):
            raise ValueError(('Krawczyk',a,b,[(str(K[i]),str(X[i])) for i in range(n) if not inside(K[i],X[i])]))
        # Bound the implicit derivative using a Neumann-series enclosure.
        # A weighted infinity norm accounts for differently scaled variables.
        rownorm=[sum(abs(B[i][j])*arb(str(rad[j]))/arb(str(rad[i])) for j in range(n)) for i in range(n)]
        eta_exact=arb(max(r.upper() for r in rownorm));eta=float(eta_exact)
        if not eta_exact<1:raise ValueError(('contraction',a,b,str(eta_exact)))
        pv=dp(pp,X)
        v=[-sum(Af[i][k]*pv[k] for k in range(n)) for i in range(n)]
        vnorm=arb(max((abs(v[i])/arb(str(rad[i]))).upper() for i in range(n)))
        err=eta_exact*vnorm/(1-eta_exact)
        deriv=[arb(v[i],err*arb(str(rad[i]))) for i in range(n)]
        if prev is not None:
            assert all(inside(prev[i],X[i]) for i in range(n)),('branch connection',a)
        # The root at the right endpoint is in this tighter enclosure. Its
        # inclusion in the next box identifies the same analytic branch.
        endp=arb(str(b));endf=fi(endp,zz)
        prev=[zz[i]-sum(Af[i][k]*endf[k] for k in range(n))+sum(B[i][j]*(X[j]-zz[j]) for j in range(n)) for i in range(n)]
        if check:check(pp,X,deriv)
        accepted.append(dict(a=float(a),b=float(b),center=z.tolist(),radii=rad.tolist(),contraction=eta))
        if len(accepted)%50==0:print(name,'certified',len(accepted),'through',b,flush=True)
    return accepted
