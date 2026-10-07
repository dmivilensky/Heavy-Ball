"""Polynomial systems of the three continued families (Appendix B, Table 2).

L(m), U(m): the normalized lag functionals in (G, u, c). p1(): equations p1quartic, p1system with the
coordinates p1blowup; p2(): equation p2system; central(): equation p5system with the coordinates
blowup, divided by h before clearing denominators. Each function returns the parameter, the unknowns,
the polynomial system F, the dual polynomial Q and the substitution."""
import sympy as s
G,u,c,z,h=s.symbols('G u c z h'); A,B,y,P,D,w2,w3=s.symbols('A B y P D w2 w3')
def L(m):return ((1-c*c-u*u)*(s.chebyshevt(m,G)-s.chebyshevt(m,c))-2*u*(1-c*c)*s.chebyshevu(m-1,c))/s.chebyshevu(m-1,G)
def U(m):return ((1-c*c-u*u)*(s.chebyshevt(m,G)-s.chebyshevt(m,c))+2*u*(1-c*c)*s.chebyshevu(m-1,c))/s.chebyshevu(m-1,G)
def p1():
    Q=-((1-c*c-u*u)*(2+z*c+z*z*(c*c-1))-2*u*z*(1-c*c)*(1+z*c))
    su={u:1-z/2+z*z*A,c:z*z*B}
    F=[s.cancel(e.subs(su)/z**2).expand() for e in (Q,s.diff(Q,c))]
    return z,[A,B],F,Q,su
def central():
    # All equations are polynomial and regular at h=G-1=0.
    x=-s.Rational(1,2)+h*B;uv=s.Rational(1,2)+h*A
    def mean(e):return (1-h*P)*e.subs(c,x)+h*P*e.subs(c,y)
    es=[mean(L(1))+h*D]+[mean(U(m))-h*D for m in (2,3,5)]
    Q=s.cancel(-(L(1)+w2*U(2)+w3*U(3)+(1-w2-w3)*U(5))/(G-c))
    Q=s.fraction(Q)[0]
    F=[s.cancel(e.subs({G:1+h,u:uv})/h).as_numer_denom()[0].expand() for e in es]
    F += [s.expand(e.subs({G:1+h,u:uv})) for e in (Q.subs(c,y),s.diff(Q,c).subs(c,x),s.diff(Q,c).subs(c,y))]
    return h,[A,B,y,P,D,w2,w3],F,Q,{G:1+h,u:uv,c:x}
def p2():
    x=s.symbols('x');p=s.symbols('p');d=s.symbols('d');W=s.symbols('W')
    def mean(e):return (1-p)*e.subs(c,x)+p*e.subs(c,y)
    es=[mean(L(1))+d,mean(U(2))-d,mean(U(3))-d]
    Q=s.cancel(-(L(1)+W*U(2)+(1-W)*U(3))/(G-c)).as_numer_denom()[0]
    F=[s.cancel(e).as_numer_denom()[0].expand() for e in es]
    F += [s.expand(e) for e in (Q.subs(c,y),s.diff(Q,c).subs(c,x),s.diff(Q,c).subs(c,y))]
    return G,[u,x,y,p,d,W],F,Q,{}
if __name__=='__main__':
    import pickle
    from pathlib import Path
    for name,fun in [('p1',p1),('central',central),('p2',p2)]:
        data=fun();Path(name+'_system.pkl').write_bytes(pickle.dumps(data))
        print(name,[len(s.Poly(f,*([data[0]]+data[1])).terms()) for f in data[2]],flush=True)
