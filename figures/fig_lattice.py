"""Figure fig_lattice: the first level in lattice coordinates (xi, eta) = (w/(2 Theta), pi/Theta).

The curve C of the tuning of Proposition prop:13 (kappa = 100, gap I_4), the lattice points (j, k),
the lag-one line xi = 1 with the points (1, K), the line U: eta = 4 xi + 1 of the upper lag 4 with the
fan points (p, 4p + 1), and the lattice lines (4m + 1) xi - m eta = 1 of the lags 4m + 1 through (1, 4)
(Theorem lattice, Theorem farey, Lemma fan).
"""
import math
import numpy as np
from math import gcd
from style import plt, save, INK, INK2, MUTED, RULE, TEXTWIDTH, RED, BLUE, ORANGE, VIOLET, tint
import hbgeom as H

KAPPA, AL, BETA, K = 100.0, 70471 / 31800, 6169 / 9900, 4


def curve(n=6000):
    w = np.linspace(0.02, np.pi, n)
    th = H.theta(H.Qcoef(KAPPA, AL, BETA), w)
    return w / (2 * th), np.pi / th


def lattice(ax, xmax, ymax):
    pts = [(j, k) for k in range(1, int(ymax) + 1) for j in range(0, int(xmax) + 1) if gcd(j, k) == 1 and j <= k]
    P = np.array(pts, float)
    ax.scatter(P[:, 0], P[:, 1], s=3, color=MUTED, edgecolors='none', zorder=2)


def main():
    xi, eta = curve()
    fig, (a0, a1) = plt.subplots(1, 2, figsize=(TEXTWIDTH, 2.75), gridspec_kw=dict(width_ratios=[1.05, 1], wspace=0.25))
    # left: the whole picture
    X, Y = 3.4, 14.5
    lattice(a0, X, Y)
    a0.axvline(1, color=INK, lw=0.7, zorder=3)
    x = np.linspace(0, X, 10)
    a0.plot(x, K * x + 1, color=VIOLET, lw=0.8, zorder=3)
    for m in (1, 2, 3):
        a0.plot(x, ((K * m + 1) * x - 1) / m, color=ORANGE, lw=0.6, ls=(0, (3, 1.5)), zorder=3)
    for k in range(2, 15):
        a0.scatter([1], [k], s=9, color=INK, zorder=5)
    for p in (1, 2, 3):
        a0.scatter([p], [K * p + 1], s=16, color=VIOLET, edgecolors='white', linewidths=0.4, zorder=6)
        if p > 1:
            a0.text(p + 0.1, K * p + 1 - 0.3, r'$(%d,%d)$' % (p, K * p + 1), fontsize=6.5, color=VIOLET, va='top')
    a0.plot(xi, eta, color=RED, lw=1.4, zorder=7)
    a0.text(0.32, 2.6, r'$\mathcal{C}$', color=RED, fontsize=9)
    a0.text(1.1, 13.9, r'$\xi=1$', fontsize=7, ha='left', va='center')
    a0.text(2.35, 11.4, r'$U$', color=VIOLET, fontsize=8, ha='right')
    a0.text(3.3, 6.0, 'dashed: lags\n$5$, $9$, $13$\nthrough $(1,4)$', color=ORANGE, fontsize=6.5, ha='right')
    a0.add_patch(plt.Rectangle((0.93, 3.85), 0.14, 1.4, fill=False, ec=INK2, lw=0.6, zorder=8))
    a0.set_xlim(0, X)
    a0.set_ylim(0, Y)
    a0.set_xlabel(r'$\xi=\omega/2\Theta$')
    a0.set_ylabel(r'$\eta=\pi/\Theta$')
    a0.text(-0.14, 1.0, '(a)', transform=a0.transAxes, va='bottom')
    # right: the crossing
    x0, x1, y0, y1 = 0.93, 1.07, 3.85, 5.25
    a1.axvline(1, color=INK, lw=0.7, zorder=3)
    x = np.linspace(x0, x1, 10)
    a1.plot(x, K * x + 1, color=VIOLET, lw=0.8, zorder=3)
    for p in (2, 3, 4):
        t = np.linspace(0, 1, 5)
        a1.plot(t * p * 10, t * (K * p + 1) * 10, color=VIOLET, lw=0.5, ls=(0, (1, 1.5)), zorder=2)
    for k in (4, 5):
        a1.scatter([1], [k], s=14, color=INK if k == 4 else VIOLET, zorder=6)
        if k == 5:
            a1.text(0.995, k, r'$(1,5)$', fontsize=6.5, ha='right', va='bottom')
        else:
            a1.text(1.004, k - 0.02, r'$(1,4)$', fontsize=6.5, ha='left', va='top')
    m = (xi > x0) & (xi < x1) & (eta > y0) & (eta < y1)
    a1.plot(xi, eta, color=RED, lw=1.4, zorder=7)
    dip = xi > 1
    a1.fill_betweenx(eta[dip], 1, xi[dip], color=tint(RED, 0.2), lw=0, zorder=1)
    a1.text(1.012, 4.62, 'dip', color=RED, fontsize=6.5)
    for p in (2, 3, 4):
        r = p / (K * p + 1)
        a1.text(x1 - 0.002, x1 * (K * p + 1) / p - 0.0, r'$%d/%d$' % (p, K * p + 1), color=VIOLET, fontsize=6,
                ha='right', va='bottom') if x1 * (K * p + 1) / p < y1 else None
    a1.text(0.94, 4.83, r'$U$', color=VIOLET, fontsize=8)
    a1.text(0.946, 4.3, r'$\mathcal{C}$', color=RED, fontsize=9)
    a1.set_xlim(x0, x1)
    a1.set_ylim(y0, y1)
    a1.set_xlabel(r'$\xi$')
    a1.set_ylabel(r'$\eta$')
    a1.text(-0.16, 1.0, '(b)', transform=a1.transAxes, va='bottom')
    for a in (a0, a1):
        for side in ('top', 'right'):
            a.spines[side].set_visible(False)
    save(fig, 'fig_lattice')
    print('max xi', xi.max(), 'eta range', eta.min(), eta.max())


if __name__ == '__main__':
    main()
