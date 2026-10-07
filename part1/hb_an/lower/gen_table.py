"""Generates the closed forms A(rho) + B(rho) q of Table 5 (exact reduction modulo F)."""
import sympy as sp
from exactlib import *
lam, mu = L, M
s = (lam + q**2/lam)/(2*q); sb = swap(s)
g = 2*r*(r**2 - r + 1)/(r**2 + 1)**2; c = (r + 1/r)/2
Dl = (q/2)*(1 - 1/lam); Dlb = swap(Dl)
alp = Dlb*(s + 1)/2; alm = Dlb*(s - 1)/2
em = mu/r; ep = lam/r
b_m = g*c*em + alp; b_p = g*c*ep + alm
items = [(r'\gamma_-', symred(g*r/(4*c) - (alp + swap(alp))/2 - Dl*Dlb/4)),
         (r'\gamma_+', symred(g/(4*c*r) - (alm + swap(alm))/2 - Dl*Dlb/4)),
         (r'|\beta_-|^2', Re(b_m*swap(b_m))), (r'|\beta_+|^2', Re(b_p*swap(b_p)))]
for sg, b in [('-', b_m), ('+', b_p)]:
    for rr in range(3): items.append((r'B^{%s}_{%d}' % (sg, rr), Re(b*ep**rr)))
out = []
for nm, e in items:
    A, B = modF(e)
    out.append(r'%s&=%s\\&\quad+\Big(%s\Big)\,q,\\' % (nm, sp.latex(sp.factor(A)), sp.latex(sp.factor(B))))
open('closed_forms.tex', 'w').write('\\begin{align*}\n' + '\n'.join(out).rstrip('\\,') + '\n\\end{align*}\n')
print(open('closed_forms.tex').read())
