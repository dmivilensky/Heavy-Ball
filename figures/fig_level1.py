"""Figure fig_level1: the first level at kappa = 100, unrolled along the boundary of Omega_1.

Data: data/gap_100_refined.npz and data/gap_100_hi_refined.npz (gap_data.py, gap_refine.py).
Each grid point (beta, alpha L = alpha_d(beta) + s) is classified by the criteria of the first level:
tongues of K-cycles (Corollary C(b); labeled by the smallest K), the two-lag inequality (Theorem
fan-conv, with K the gap index), feasible fan programs of period pK + 1 (Theorem fan-div(a), p <= 79),
certificates built from the fan lags (Theorem fan-conv, lags up to 200), and undecided grid points.
The two propositions of the first level are marked.
"""
import math
import numpy as np
from matplotlib.patches import Patch
from style import plt, save, INK, INK2, MUTED, TEXTWIDTH, DATA, AQUA, BLUE, RED, ORANGE, VIOLET, MAGENTA, tint
import hbgeom as H
from gap_data import boundary

KAPPA = 100.0
C_TA = tint(BLUE, 0.5)
C_TB = tint(BLUE, 0.3)
C_TWO = tint(AQUA, 0.55)
C_FAN2 = RED
C_FAN3 = MAGENTA
C_FAN4 = VIOLET
C_CERT = ORANGE
C_J = VIOLET
C_NONE = INK
C_DIV = '#e4e3df'
POINTS = {'Prop. 9.4': (70471 / 31800, 6169 / 9900),
          'Prop. 9.5': (111 / 50, 31 / 50)}


def smallest_K(beta, aL):
    x = aL / KAPPA
    for k in range(3, 300):
        if H.E(x, aL, beta, math.cos(2 * math.pi / k)) <= 0:
            return k
    return 0


def load(name):
    g = np.load(DATA / f'{name}_refined.npz')
    rot, s, lab, beta, ad = g['rot'], g['s'], g['lab'], g['beta'], g['ad']
    o = np.argsort(rot)
    rot, lab, beta, ad = rot[o], lab[:, o], beta[o], ad[o]
    K = np.zeros_like(lab)
    for j, i in np.argwhere(lab >= 100):
        if lab[j, i] < 200:
            K[j, i] = smallest_K(beta[i], ad[i] + s[j])
    return rot, s, lab, K, beta, ad


def colorize(lab, K):
    rgb = np.ones(lab.shape + (3,))
    def put(mask, hexc):
        h = hexc.lstrip('#')
        rgb[mask] = [int(h[i:i + 2], 16) / 255 for i in (0, 2, 4)]
    put((K > 0) & (K % 2 == 0), C_TA)
    put((K > 0) & (K % 2 == 1), C_TB)
    put(lab == 1, C_TWO)
    put(lab == 2, C_CERT)
    put(lab == 3, C_J)
    put(lab == 202, C_FAN2)
    put(lab == 203, C_FAN3)
    put(lab >= 204, C_FAN4)
    put(lab == 4, C_NONE)
    put(lab == 0, C_DIV)
    return rgb


def mark_points(ax):
    for name, (aL, b) in POINTS.items():
        a0 = boundary(KAPPA, b)
        rot = math.acos(H.tangency_c(KAPPA, a0, b)) / (2 * math.pi)
        ax.plot([rot], [aL - a0], marker='o', ms=3.6, mfc='white', mec=INK, mew=0.8, zorder=10)
        yield name, rot, aL - a0


def emphasize_small(lab, mask_value_fn, rad=1):
    """Thicken isolated special points by one cell so that they remain visible at print size."""
    out = lab.copy()
    sel = mask_value_fn(lab)
    ii, jj = np.nonzero(sel)
    for di in range(-rad, rad + 1):
        for dj in range(-rad, rad + 1):
            i2 = np.clip(ii + di, 0, lab.shape[0] - 1)
            j2 = np.clip(jj + dj, 0, lab.shape[1] - 1)
            keep = (out[i2, j2] == 1)
            out[i2[keep], j2[keep]] = lab[ii[keep], jj[keep]]
    return out


def main():
    fig = plt.figure(figsize=(TEXTWIDTH, 4.9))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.15], hspace=0.38, wspace=0.22)
    ax = fig.add_subplot(gs[0, :])
    rot, s, lab, K, beta, ad = load('gap_100')
    lab_v = emphasize_small(lab, lambda L: (L == 4) | (L == 2))
    ax.pcolormesh(rot, s, colorize(lab_v, K), shading='nearest', rasterized=True)
    for k in range(4, 12):
        ax.text(1 / k, 0.0465, r'$%d$' % k, fontsize=6.5, ha='center', va='center')
    for name, r0, s0 in mark_points(ax):
        pass
    ax.set_xlim(rot.min(), rot.max())
    ax.set_ylim(0, s.max())
    ax.set_xlabel(r'rotation number $\omega_*/2\pi$ of the tangency')
    ax.set_ylabel(r'$\alpha L-\alpha_{\partial}(\beta)$')
    ax.set_title('(a) the first level', fontsize=8, loc='left')
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    rot, s, lab, K, beta, ad = load('gap_100_hi')
    lab_v = emphasize_small(lab, lambda L: (L == 4) | (L == 2))
    for col, (lo, hi, Kg, title) in enumerate(((0.2, 0.25, 4, r'(b) the gap $I_4$'), (0.25, 0.3268, 3, r'(c) the gap $I_3$'))):
        a = fig.add_subplot(gs[1, col])
        m = (rot >= lo - 0.004) & (rot <= hi + 0.004)
        a.pcolormesh(rot[m], s, colorize(lab_v[:, m], K[:, m]), shading='nearest', rasterized=True)
        for name, r0, s0 in mark_points(a):
            if lo <= r0 <= hi:
                dx = -0.0085 if name.endswith('9.4') else 0.0085
                a.annotate(name, xy=(r0, s0), xytext=(r0 + dx, 0.022), fontsize=6.5, ha='center',
                           bbox=dict(boxstyle='round,pad=0.15', fc='white', ec='none', alpha=0.85),
                           arrowprops=dict(arrowstyle='-', lw=0.5, color=INK))
        for pp, (dx, dy) in ((2, (-0.010, 0.016)), (3, (0.011, 0.02))):
            mm = lab[:, m] == 200 + pp
            if mm.sum() < 3:
                continue
            jj, ii = np.nonzero(mm)
            r0, s0 = rot[m][int(np.median(ii))], s[int(np.median(jj))]
            a.annotate(r'$%d/%d$' % (pp, pp * Kg + 1), xy=(r0, s0), xytext=(r0 + dx, s0 + dy), fontsize=6.5,
                       ha='center', bbox=dict(boxstyle='round,pad=0.12', fc='white', ec='none', alpha=0.85), arrowprops=dict(arrowstyle='-', lw=0.5, color=INK2))
        a.set_xlim(max(lo, rot.min()), min(hi, rot.max()))
        a.set_ylim(0, 0.16 if Kg == 3 else 0.11)
        a.set_xlabel(r'$\omega_*/2\pi$')
        if col == 0:
            a.set_ylabel(r'$\alpha L-\alpha_{\partial}(\beta)$')
        a.set_title(title, fontsize=8, loc='left')
        for side in ('top', 'right'):
            a.spines[side].set_visible(False)
    handles = [Patch(fc=C_TA, label=r'tongues of $K$-cycles'), Patch(fc=C_TWO, label='two-lag inequality'),
               Patch(fc=C_FAN2, label=r'fan program, $p=2$'), Patch(fc=C_FAN3, label=r'fan program, $p=3$'),
               Patch(fc=C_FAN4, label=r'fan program, $p\geq4$'), Patch(fc=C_CERT, label='other fan-lag certificates'),
               Patch(fc=C_NONE, label='undecided at this resolution'), Patch(fc=C_DIV, label=r'$\alpha L\geq2(1+\beta)$')]
    fig.legend(handles=handles, loc='lower center', ncol=4, bbox_to_anchor=(0.5, -0.05), fontsize=6.5)
    save(fig, 'fig_level1')


if __name__ == '__main__':
    main()
