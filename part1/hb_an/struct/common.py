"""Common exact objects for the structural checks (sympy, exact rationals).

X = e^t = rho, Y = e^v (v = t - log q, so q = X/Y), c = cos(omega),
a_m = g (cosh mt - cos m omega), b_m = 2 sinh v sin(omega) sin(m omega), g = sin^2 omega - sinh^2 v.
Normalized C-free functionals with r_m = sinh t/sinh mt (the functionals L_m, U_m of the paper):
  lower m: r_m (a_m - b_m), coefficient of C equal to +sinh t;
  upper m: r_m (a_m + b_m), coefficient of C equal to -sinh t."""
import sympy as sp
X, Y, c, W, q = sp.symbols('X Y c W q')
shX = lambda m: (X**m - X**-m)/2
chX = lambda m: (X**m + X**-m)/2
shv = (Y - 1/Y)/2
g = (1 - c**2) - shv**2
def Tch(m): return sp.chebyshevt(m, c)
def Uch(m): return sp.chebyshevu(m, c)      # U_{m}
def a(m): return g*(chX(m) - Tch(m))
def b(m): return 2*shv*(1 - c**2)*Uch(m-1)
def r(m): return shX(1)/shX(m)
def Lfun(m): return sp.together(r(m)*(a(m) - b(m)))
def Ufun(m): return sp.together(r(m)*(a(m) + b(m)))
# spiral curve (F2 in (X,Y)):  (1+X^2)(Y^2-2Y-1) + 2XY = 0   <=>  sh v = 1 - 1/(2 ch t)
S = (1 + X**2)*(Y**2 - 2*Y - 1) + 2*X*Y
F = lambda rho, qq: rho**2*(rho**2 + 1) - 2*qq*rho*(rho**2 - rho + 1) - qq**2*(rho**2 + 1)
H = lambda rho, qq: (rho - qq)**2*(rho**2 + 1) - 2*qq*(qq - (1 - qq)*rho**2)
def numer(e): return sp.numer(sp.factor(sp.together(e)))
