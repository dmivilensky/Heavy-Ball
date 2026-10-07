"""Figure fig_spectra: worst-case spectra and trajectories for Polyak's tuning, one column per regime.

For each condition number the optimal family of Theorem classification is evaluated at
G = cosh T, T = log rho_HB (part1/rate.py). Row 1: the iterates e_k, 0 <= k <= 18, in the plane of
the main atom, e_k = sqrt(nu_1) z_1^k with z_1 = rho e^{i omega_1} (Lemma startup); the
angle is k omega_1 and the radius is drawn as rho^{k/12}, so that growth and decay are visible on
one scale (the dotted circle is radius one). Row 2: the same for the second atom. Row 3: the
certificate polynomial Q(c) = -Phi(c)/(G - c) of equation dualpoly with the optimal weights, and
the atoms of the optimal spectrum at its double zeros (area proportional to the weight; an open
circle is an atom of weight zero).
"""
import math
import numpy as np
import sympy as sp
from style import plt, save, FAMILY, FAMILY_TEX, INK, INK2, MUTED, RULE, TEXTWIDTH, tint
import rate

KC = 9 + 4 * 5 ** 0.5
COLUMNS = [(2.0, r'$\kappa=2$'), (5.0, r'$\kappa=5$'), (8.3, r'$\kappa=8.3$'), (KC, r'$\kappa=\kappa_c$'),
           (30.0, r'$\kappa=30$'), (100.0, r'$\kappa=100$')]
NSTEP = 18
SQUASH = 12.0

cs, Gs, us, w2s, w3s, w5s = sp.symbols('c G u w2 w3 w5')


def _fun(m, sign):
    g = 1 - cs ** 2 - us ** 2
    return (g * (sp.chebyshevt(m, Gs) - sp.chebyshevt(m, cs)) + sign * 2 * us * (1 - cs ** 2) * sp.chebyshevu(m - 1, cs)) / sp.chebyshevu(m - 1, Gs)


PHI = _fun(1, -1) + w2s * _fun(2, 1) + w3s * _fun(3, 1) + (1 - w2s - w3s) * _fun(5, 1)
QEXPR = sp.cancel(sp.together(-PHI / (Gs - cs)))
assert sp.degree(sp.denom(QEXPR), cs) == 0
_Q = sp.lambdify((cs, Gs, us, w2s, w3s), QEXPR, 'numpy')
QPOLY = lambda c, G, u, w2, w3, w5: _Q(c, G, u, w2, w3)


def spectrum(kappa):
    r = rate.evaluate(kappa)
    T = r['log_rate']
    G = math.cosh(T)
    fam = r['family']
    if fam == 'P1':
        z = 1 / G
        A, B = rate.polynomial_solution('p1', z)
        u = 1 - z / 2 + z * z * A
        atoms = [(z * z * B, 1.0)]
        w = (0.0, 1.0, 0.0)
    elif fam == 'S':
        u = 1 - 1 / (2 * G)
        W = 2 * (1 + 2 * G - 2 * G * G) / (4 * G + 1)
        atoms = [(-1 / (2 * G), 1.0)]
        w = (W, 1 - W, 0.0)
    elif fam == 'P2':
        u, x, y, p, d, W = rate.polynomial_solution('p2', G)
        atoms = [(x, 1 - p), (y, p)]
        w = (W, 1 - W, 0.0)
    else:
        h = G - 1
        A, B, y, P, D, w2, w3 = rate.polynomial_solution('central', h)
        u = 0.5 + h * A
        atoms = [(-0.5 + h * B, 1 - h * P), (y, h * P)]
        w = (w2, w3, 1 - w2 - w3)
    return dict(kappa=kappa, T=T, rho=math.exp(T), G=G, u=u, fam=fam, atoms=atoms, w=w)


def plane(ax, rho, omega, color, first):
    k = np.arange(NSTEP + 1)
    r = rho ** (k / SQUASH)
    x, y = r * np.cos(k * omega), r * np.sin(k * omega)
    th = np.linspace(0, 2 * np.pi, 200)
    ax.plot(np.cos(th), np.sin(th), color=RULE, lw=0.5, ls=(0, (1, 1.5)), zorder=1)
    for i in range(NSTEP):
        ax.plot(x[i:i + 2], y[i:i + 2], color=tint(color, 0.25 + 0.6 * i / NSTEP), lw=0.45, zorder=2)
    shades = [tint(color, 0.3 + 0.7 * i / NSTEP) for i in k]
    ax.scatter(x, y, s=6, c=shades, edgecolors='none', zorder=3)
    ax.scatter([x[0]], [y[0]], s=22, facecolors='none', edgecolors=INK, linewidths=0.6, zorder=4)
    lim = max(1.0, r.max()) * 1.12
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.text(0, -lim * 1.0, r'$\omega=%.3f\pi$' % (omega / np.pi), ha='center', va='top', fontsize=6.5, color=INK2)


def main():
    fig = plt.figure(figsize=(TEXTWIDTH, 3.6))
    gs = fig.add_gridspec(3, len(COLUMNS), height_ratios=[1, 1, 0.78], hspace=0.55, wspace=0.14)
    for j, (kappa, title) in enumerate(COLUMNS):
        sp_ = spectrum(kappa)
        color = FAMILY[sp_['fam']]
        atoms = sorted(sp_['atoms'], key=lambda a: -a[1])
        ax = fig.add_subplot(gs[0, j])
        plane(ax, sp_['rho'], math.acos(max(-1, min(1, atoms[0][0]))), color, j == 0)
        ax.text(0.5, 1.30, title, transform=ax.transAxes, ha='center', va='bottom', fontsize=7.5, color=INK)
        name = FAMILY_TEX[sp_['fam']] + (' (cycle)' if kappa == KC else '')
        ax.text(0.5, 1.15, name, transform=ax.transAxes, ha='center', va='bottom', fontsize=7, color=color)
        ax.text(0.5, 1.0, r'$\rho_{\rm HB}=%.3f$' % sp_['rho'], transform=ax.transAxes, ha='center', va='bottom',
                fontsize=6.5, color=INK2)
        ax2 = fig.add_subplot(gs[1, j])
        if len(atoms) > 1 and atoms[1][1] > 1e-9:
            plane(ax2, sp_['rho'], math.acos(atoms[1][0]), color, False)
        elif len(atoms) > 1:
            ax2.axis('off')
            ax2.text(0.5, 0.5, 'second atom\nof weight zero', ha='center', va='center', fontsize=6.5, color=MUTED,
                     transform=ax2.transAxes)
        else:
            ax2.axis('off')
            ax2.text(0.5, 0.5, 'one frequency', ha='center', va='center', fontsize=6.5, color=MUTED,
                     transform=ax2.transAxes)
        ax3 = fig.add_subplot(gs[2, j])
        c = np.linspace(-1, 1, 600)
        Q = QPOLY(c, sp_['G'], sp_['u'], *sp_['w'])
        ax3.plot(c, Q, color=color, lw=1.1, zorder=3)
        ax3.axhline(0, color=RULE, lw=0.6, zorder=1)
        for (ca, wa) in atoms:
            if wa > 1e-9:
                ax3.scatter([ca], [0], s=10 + 40 * wa, color=color, edgecolors=INK, linewidths=0.4, zorder=5)
            else:
                ax3.scatter([ca], [0], s=14, facecolors='white', edgecolors=INK, linewidths=0.6, zorder=5)
        top = Q.max()
        ax3.set_ylim(-0.12 * top, 1.08 * top)
        ax3.set_xlim(-1, 1)
        ax3.set_yticks([])
        ax3.set_xticks([-1, 0, 1])
        if j:
            ax3.set_xticklabels([])
        for side in ('top', 'right', 'left'):
            ax3.spines[side].set_visible(False)
        if j == 0:
            ax3.set_ylabel(r'$Q(c)$', fontsize=7.5)
        print(kappa, sp_['fam'], 'rho=%.6f' % sp_['rho'], 'atoms', [(round(a, 5), round(b, 5)) for a, b in atoms],
              'w', [round(x, 5) for x in sp_['w']])
    fig.text(0.5, 0.035, r'$c=\cos\omega$', ha='center', va='bottom', fontsize=7.5)
    save(fig, 'fig_spectra')


if __name__ == '__main__':
    main()
