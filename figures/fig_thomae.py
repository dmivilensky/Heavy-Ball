"""Figure fig_thomae: the angle function Theta against the points (2 pi j/k, pi/k) (Theorem A).

Three tunings at kappa = 100: a tuning in Omega_1 and a tuning in the tongue of K = 5, both at the
momentum beta_5 of the phase-locked point P_5 (tangency at 2 pi/5), and the tuning of
Proposition prop:13 in the gap I_4. Theta is computed from cot Theta = Xi(cos w)/sin w (Lemma angle).
"""
import math
import numpy as np
from math import gcd
from scipy.optimize import brentq
from style import plt, save, clean, INK, INK2, MUTED, RULE, TEXTWIDTH, AQUA, BLUE, RED, tint
import hbgeom as H
from gap_data import boundary

KAPPA = 100.0


def phase_locked_beta(K):
    f = lambda b: H.tangency_c(KAPPA, boundary(KAPPA, b), b) - math.cos(2 * math.pi / K)
    return brentq(f, 0.45, 0.92, xtol=1e-12)


def main():
    b5 = phase_locked_beta(5)
    a5 = boundary(KAPPA, b5)
    tunings = [('in $\\Omega_1$', a5 - 0.06, b5, AQUA), ('tongue $K=5$', a5 + 0.05, b5, BLUE),
               ('gap $I_4$', 70471 / 31800, 6169 / 9900, RED)]
    w = np.linspace(1e-3, np.pi, 4000)
    fig = plt.figure(figsize=(TEXTWIDTH, 2.75))
    ax = fig.add_axes([0.08, 0.15, 0.6, 0.8])
    pts = [(2 * np.pi * j / k, np.pi / k, j) for k in range(2, 41) for j in range(1, k // 2 + 1)
           if gcd(j, k) == 1 and 2 * j <= k]
    P = np.array(pts)
    m1 = P[:, 2] == 1
    ax.scatter(P[~m1, 0], P[~m1, 1], s=3, color=MUTED, edgecolors='none', zorder=3)
    ax.scatter(P[m1, 0], P[m1, 1], s=7, color=INK, edgecolors='none', zorder=4)
    ax.plot([0, np.pi], [0, np.pi / 2], color=INK2, lw=0.6, ls=(0, (3, 2)), zorder=2)
    ax.text(2.7, 2.7 / 2 - 0.08, r'$\omega/2$', color=INK2, fontsize=7, ha='left', va='top')
    ax.text(np.pi - 0.03, np.pi / 2 - 0.03, r'$2$', fontsize=6, color=INK, ha='right', va='top')
    for K in range(3, 9):
        ax.text(2 * np.pi / K + 0.03, np.pi / K - 0.03, r'$%d$' % K, fontsize=6, color=INK, ha='left', va='top')
    curves = []
    for name, aL, b, col in tunings:
        th = H.theta(H.Qcoef(KAPPA, aL, b), w)
        ax.plot(w, th, color=col, lw=1.2, zorder=5)
        curves.append((name, th, col))
    ax.text(0.55, 1.62, tunings[0][0], color=AQUA, fontsize=7, ha='left')
    ax.text(0.55, 1.5, tunings[1][0], color=BLUE, fontsize=7, ha='left')
    ax.text(0.55, 1.38, tunings[2][0], color=RED, fontsize=7, ha='left')
    ax.set_xlim(0, np.pi)
    ax.set_ylim(0, 1.75)
    ax.set_xticks([0, np.pi / 4, np.pi / 2, 3 * np.pi / 4, np.pi])
    ax.set_xticklabels(['0', r'$\pi/4$', r'$\pi/2$', r'$3\pi/4$', r'$\pi$'])
    ax.set_yticks([0, np.pi / 8, np.pi / 4, 3 * np.pi / 8, np.pi / 2])
    ax.set_yticklabels(['0', r'$\pi/8$', r'$\pi/4$', r'$3\pi/8$', r'$\pi/2$'])
    ax.set_xlabel(r'frequency $\omega$')
    ax.set_ylabel(r'$\Theta(\omega)$')
    x0, x1, y0, y1 = 1.18, 1.64, 0.56, 0.86
    ax.add_patch(plt.Rectangle((x0, y0), x1 - x0, y1 - y0, fill=False, ec=INK2, lw=0.6, zorder=6))
    clean(ax)
    # inset
    ins = fig.add_axes([0.73, 0.2, 0.26, 0.72])
    ww = np.linspace(x0, x1, 1500)
    for name, aL, b, col in tunings:
        ins.plot(ww, H.theta(H.Qcoef(KAPPA, aL, b), ww), color=col, lw=1.2, zorder=5)
    ins.plot([x0, x1], [x0 / 2, x1 / 2], color=INK2, lw=0.6, ls=(0, (3, 2)), zorder=2)
    for K in (4, 5):
        ins.scatter([2 * np.pi / K], [np.pi / K], s=12, color=INK, zorder=6)
        ins.text(2 * np.pi / K + 0.01, np.pi / K - 0.012, r'$(2\pi/%d,\pi/%d)$' % (K, K), fontsize=6, ha='left', va='top')
    th = H.theta(H.Qcoef(KAPPA, 70471 / 31800, 6169 / 9900), ww)
    dip = th < ww / 2
    ins.fill_between(ww, th, ww / 2, where=dip, color=tint(RED, 0.25), lw=0, zorder=1)
    ins.text(1.43, 0.665, 'dip', color=RED, fontsize=6.5, ha='center')
    ins.set_xlim(x0, x1)
    ins.set_ylim(y0, y1)
    ins.set_xticks([2 * np.pi / 5, np.pi / 2])
    ins.set_xticklabels([r'$2\pi/5$', r'$\pi/2$'])
    ins.set_yticks([np.pi / 5, np.pi / 4])
    ins.set_yticklabels([r'$\pi/5$', r'$\pi/4$'])
    for side in ('top', 'right'):
        ins.spines[side].set_visible(False)
    save(fig, 'fig_thomae')
    print('beta_5 =', b5, 'alpha_d =', a5)


if __name__ == '__main__':
    main()
