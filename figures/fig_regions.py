"""Figure fig_regions: the regions of the phase diagram at kappa = 100, each at its own scale.

(a) Omega_1 and its parts (data/phase_100.npz). (b) The tongues of K-cycles, labeled by the smallest
K (Corollary C(b)). (c) The belt along the boundary of Omega_1 around the phase-locked points
P_3, P_4, P_5, classified on a fine grid by the same closed forms. (d) The gaps unrolled along the
boundary of Omega_1 (data/gap_100.npz): horizontal coordinate the rotation number of the tangency,
vertical coordinate the distance in alpha L beyond the boundary.
"""
import math
import numpy as np
from scipy.optimize import brentq
from matplotlib.patches import Patch
from style import plt, save, INK, INK2, MUTED, TEXTWIDTH, DATA, AQUA, YELLOW, BLUE, RED, tint
import hbgeom as H
from gap_data import boundary
from fig_phase import colorize, C_DIV, C_Q, C_UD, C_REST, C_T1, C_T2, C_GAP

KAPPA = 100.0


def phase_locked(K):
    f = lambda b: H.tangency_c(KAPPA, boundary(KAPPA, b), b) - math.cos(2 * math.pi / K)
    b = brentq(f, 0.25, 0.95, xtol=1e-12)
    return boundary(KAPPA, b), b


def classify_window(a0, a1, b0, b1, n=900):
    aL = np.linspace(a0, a1, n)
    be = np.linspace(b0, b1, n)
    Y, Bt = np.meshgrid(aL, be)
    lab = np.full(Y.shape, 4, np.int16)
    div = Y >= 2 * (1 + Bt)
    om = H.in_omega1(KAPPA, Y, Bt)
    lab[om] = 3
    lab[om & H.in_Qregion(KAPPA, Y, Bt)] = 1
    lab[om & H.underdamped(KAPPA, Y, Bt)] = 2
    x = Y / KAPPA
    K = np.zeros(Y.shape, np.int16)
    for k in range(3, 300):
        sel = ~om & ~div & (K == 0) & (H.E(x, Y, Bt, math.cos(2 * math.pi / k)) <= 0)
        K[sel] = k
    lab[K > 0] = 10 + K[K > 0]
    lab[div] = 0
    return aL, be, lab


def smallest_K(beta, aL):
    x = aL / KAPPA
    for k in range(3, 300):
        if H.E(x, aL, beta, math.cos(2 * math.pi / k)) <= 0:
            return k
    return 0


def main():
    fig, axes = plt.subplots(2, 2, figsize=(TEXTWIDTH, 4.6), gridspec_kw=dict(hspace=0.36, wspace=0.26))
    d = np.load(DATA / 'phase_100.npz')
    aL, be, lab = d['aL'], d['beta'], d['lab']
    # (a) Omega_1
    ax = axes[0, 0]
    lab_a = np.where(np.isin(lab, (1, 2, 3)), lab, -1)
    rgb = colorize(lab_a)
    rgb[lab_a == -1] = 1.0
    ax.imshow(rgb, origin='lower', extent=(0, aL[-1], be[0], be[-1]), aspect='auto', interpolation='nearest')
    ax.contour(aL, be, np.isin(lab, (1, 2, 3)).astype(float), levels=[0.5], colors=INK, linewidths=0.6)
    b = np.linspace(0, 0.999, 500)
    ax.plot(KAPPA * (1 - np.sqrt(b)) ** 2, b, color=INK, lw=0.6, ls=(0, (3, 1.5)))
    ax.plot((1 + np.sqrt(b)) ** 2, b, color=INK, lw=0.6, ls=(0, (1, 1.2)))
    ax.text(1.0, 0.35, r'$\mathcal{Q}$', fontsize=9, ha='center')
    ax.annotate(r'$\rho_{\rm q}<\rho\leq r_1$', xy=(0.75, 0.885), xytext=(0.45, 0.66), fontsize=6.5, ha='center',
                arrowprops=dict(arrowstyle='-', lw=0.5, color=INK2))
    ax.annotate(r'$\rho_{\rm q}\leq\rho\leq r_1$', xy=(2.42, 0.42), xytext=(2.75, 0.12), fontsize=6.5, ha='center',
                arrowprops=dict(arrowstyle='-', lw=0.5, color=INK2))
    ax.set_xlim(0, 2.8)
    ax.set_ylim(0, 1)
    ax.set_xlabel(r'$\alpha L$')
    ax.set_ylabel(r'$\beta$')
    ax.set_title(r'(a) $\Omega_1$ and its parts', fontsize=8, loc='left')
    # (b) tongues
    ax = axes[0, 1]
    lab_b = np.where(lab >= 10, lab, -1)
    rgb = colorize(lab_b)
    rgb[lab_b == -1] = 1.0
    rgb[lab == 0] = [v / 255 for v in (228, 227, 223)]
    ax.imshow(rgb, origin='lower', extent=(0, aL[-1], be[0], be[-1]), aspect='auto', interpolation='nearest')
    ax.contour(aL, be, np.isin(lab, (1, 2, 3)).astype(float), levels=[0.5], colors=INK, linewidths=0.6)
    for K, (x0, y0) in {3: (3.3, 0.86), 4: (2.75, 0.86), 5: (2.15, 0.88), 6: (1.68, 0.9), 7: (1.36, 0.92),
                        8: (1.12, 0.935)}.items():
        ax.text(x0, y0, r'$%d$' % K, fontsize=7, ha='center', va='center')
    ax.set_xlim(0.8, 4)
    ax.set_ylim(0.6, 1)
    ax.set_xlabel(r'$\alpha L$')
    ax.set_title(r'(b) tongues of $K$-cycles', fontsize=8, loc='left')
    # (c) the belt around P_3, P_4, P_5
    ax = axes[1, 0]
    P = {K: phase_locked(K) for K in (3, 4, 5)}
    a0, a1, b0, b1 = 1.85, 2.66, 0.33, 0.72
    wa, wb, wl = classify_window(a0, a1, b0, b1)
    rgb = colorize(wl)
    ax.imshow(rgb, origin='lower', extent=(a0, a1, b0, b1), aspect='auto', interpolation='nearest')
    for K, (pa, pb) in P.items():
        ax.plot([pa], [pb], 'o', ms=3.4, color=INK, zorder=5)
        ax.text(pa + 0.012, pb + 0.004, r'$P_{%d}$' % K, fontsize=7, ha='left', va='bottom')
    for K, (x0, y0) in {3: (2.635, 0.42), 4: (2.5, 0.64), 5: (2.33, 0.695)}.items():
        ax.text(x0, y0, r'$%d$' % K, fontsize=7, ha='center', va='center')
    ax.text(2.19, 0.5, r'$\Omega_1$', fontsize=8, ha='center')
    ax.annotate(r'$I_3$', xy=(2.52, 0.51), xytext=(2.58, 0.57), fontsize=7, color=RED,
                arrowprops=dict(arrowstyle='-', lw=0.5, color=RED))
    ax.annotate(r'$I_4$', xy=(2.24, 0.637), xytext=(2.17, 0.67), fontsize=7, color=RED,
                arrowprops=dict(arrowstyle='-', lw=0.5, color=RED))
    ax.set_xlim(a0, a1)
    ax.set_ylim(b0, b1)
    ax.set_xlabel(r'$\alpha L$')
    ax.set_ylabel(r'$\beta$')
    ax.set_title(r'(c) the belt along $\partial\Omega_1$', fontsize=8, loc='left')
    # (d) unrolled gaps
    ax = axes[1, 1]
    g = np.load(DATA / 'gap_100.npz')
    rot, s, gl, gb, gad = g['rot'], g['s'], g['lab'], g['beta'], g['ad']
    o = np.argsort(rot)
    rot, gl, gb, gad = rot[o], gl[:, o], gb[o], gad[o]
    K = np.zeros_like(gl)
    for j, i in np.argwhere(gl >= 100):
        K[j, i] = smallest_K(gb[i], gad[i] + s[j])
    show = np.full(gl.shape, -1, np.int16)
    show[K > 0] = 10 + K[K > 0]
    show[(gl > 0) & (gl < 100)] = 4
    show[gl >= 200] = 4
    show[gl == 0] = 0
    rgb = colorize(show)
    rgb[show == -1] = 1.0
    ax.pcolormesh(rot, s, rgb[:, :-1] if False else rgb, shading='nearest', rasterized=True)
    for k in range(4, 10):
        ax.text(1 / k, 0.046, r'$%d$' % k, fontsize=6.5, ha='center', va='center')
    ax.set_xlim(rot.min(), rot.max())
    ax.set_ylim(0, s.max())
    ax.set_xlabel(r'rotation number $\omega_*/2\pi$')
    ax.set_ylabel(r'$\alpha L-\alpha_{\partial}(\beta)$')
    ax.set_title(r'(d) the gaps, unrolled along $\partial\Omega_1$', fontsize=8, loc='left')
    for a in axes.flat:
        for side in ('top', 'right'):
            a.spines[side].set_visible(False)
    handles = [Patch(fc=C_Q, label=r'$\mathcal{Q}$'), Patch(fc=C_UD, label=r'$\Omega_1$, between damping curves'),
               Patch(fc=C_REST, label=r'$\Omega_1$, rest'), Patch(fc=C_T1, label=r'tongues of $K$-cycles'),
               Patch(fc=C_GAP, label='gaps'), Patch(fc=C_DIV, label=r'$\alpha L\geq2(1+\beta)$')]
    fig.legend(handles=handles, loc='lower center', ncol=3, bbox_to_anchor=(0.5, -0.06), fontsize=6.5)
    save(fig, 'fig_regions')
    print({K: P[K] for K in P})


if __name__ == '__main__':
    main()
