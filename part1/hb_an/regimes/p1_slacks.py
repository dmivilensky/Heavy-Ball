"""Proposition 43 (family P1): positivity of the finite-lag raw slacks, equation p1slackpoly.

For 1 <= m <= 12 the slack polynomials, after c_1 = -theta z/2, have positive tensor Bernstein
coefficients on [0, 3/4] x [0, 1] (bernstein.py); the tail m >= 13 uses eta >= z^3/4 and T_m(1/z) >= z^{-m}."""
import sys,sympy as s,json
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'upper'))
from bernstein import certify_positive
z,c,a=s.symbols('z c a');DD=2+z*c+z*z*(c*c-1);EE=1+z*c+z*z*(2*c*c-1)
for m in range(1,13):
 for sign in (-1,1):
  if (m,sign) in ((1,-1),(2,1),(3,1)):
    # (B.8): l_1 and v_3 vanish identically; v_2 = -2z(2c_1+z)(c_1 z-1) >= 0 for -z/2 <= c_1 <= 0
    p=s.cancel(z**(m-1)*(z*(1+z*c)*(s.chebyshevt(m,1/z)-s.chebyshevt(m,c))+sign*DD*s.chebyshevu(m-1,c)-sign*EE*s.chebyshevu(m-1,1/z)))
    ok=s.expand(p-(-2*z*(2*c+z)*(c*z-1) if (m,sign)==(2,1) else 0))==0
    print('P1 slack',m,sign,'closed form:',ok,flush=True);assert ok
    continue
  p=s.cancel(z**(m-1)*(z*(1+z*c)*(s.chebyshevt(m,1/z)-s.chebyshevt(m,c))+sign*DD*s.chebyshevu(m-1,c)-sign*EE*s.chebyshevu(m-1,1/z)))
  if sign==1 and m>=4:p=s.cancel(p/z**2)
  p=s.expand(p.subs(c,-a*z/2))
  ok=certify_positive(p,(z,a),[(0,s.Rational(3,4)),(0,1)])
  print('P1 slack',m,sign,ok,flush=True);assert ok[0]
# The infinite tail uses eta >= z^3/4 and T_m(1/z)>=z^{-m}.
assert s.Rational(3,4)**(-10)/4>4
