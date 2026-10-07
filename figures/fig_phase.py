"""Figure fig_phase: the phase diagram in (alpha L, beta) for kappa = 2, 10, 100, 1000.

Each grid point of data/phase_<kappa>.npz (phase_data.py) is classified by the closed forms of
Corollary C: alpha L >= 2(1+beta); the quadratic regime Q; Omega_1 between the damping curves;
the rest of Omega_1; the tongues of K-cycles (j = 1), labeled by the smallest such K; gaps.
"""
import math
import numpy as np
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
from style import plt, save, INK, INK2, TEXTWIDTH, DATA, AQUA, YELLOW, BLUE, VIOLET, RED, tint

C_DIV = '#e4e3df'
C_Q = tint(AQUA, 0.62)
C_UD = tint(YELLOW, 0.55)
C_REST = tint(AQUA, 0.22)
C_T1 = tint(BLUE, 0.52)
C_T2 = tint(BLUE, 0.32)
C_J2 = tint(VIOLET, 0.6)
C_GAP = RED


def colorize(lab):
    rgb = np.zeros(lab.shape + (3,))
    def put(mask, hexc):
        h = hexc.lstrip('#')
        rgb[mask] = [int(h[i:i + 2], 16) / 255 for i in (0, 2, 4)]
    put(lab == 0, C_DIV)
    put(lab == 1, C_Q)
    put(lab == 2, C_UD)
    put(lab == 3, C_REST)
    put((lab >= 10) & (lab % 2 == 0), C_T1)
    put((lab >= 10) & (lab % 2 == 1), C_T2)
    put(lab == 5, C_J2)
    put(lab == 4, C_GAP)
    return rgb


def overlay(ax, kappa, aLmax=4.0):
    b = np.linspace(0, 0.999, 600)
    ax.plot(kappa * (1 - np.sqrt(b)) ** 2, b, color=INK, lw=0.6, ls=(0, (3, 1.5)))
    ax.plot((1 + np.sqrt(b)) ** 2, b, color=INK, lw=0.6, ls=(0, (1, 1.2)))
    ax.plot(2 * (1 + b), b, color=INK2, lw=0.5)
    q = (math.sqrt(kappa) - 1) / (math.sqrt(kappa) + 1)
    ax.plot([(1 + q) ** 2], [q * q], marker='*', ms=6, color=INK, mec='white', mew=0.4, zorder=10)


def label_tongues(ax, aL, be, lab, Kmax):
    for K in range(3, Kmax + 1):
        m = lab == 10 + K
        if m.sum() < 600:
            continue
        rows = np.nonzero(m.any(axis=1) & (be < 0.9))[0]
        if len(rows) == 0:
            continue
        width = m[rows].sum(axis=1)
        i0 = rows[int(np.argmax(width))]
        jj = np.nonzero(m[i0])[0]
        ax.text(aL[int(np.median(jj))], be[i0], r'$%d$' % K, fontsize=6.5, color=INK, ha='center', va='center')


def main():
    fig, axes = plt.subplots(2, 2, figsize=(TEXTWIDTH, 4.35), sharex=True, sharey=True,
                             gridspec_kw=dict(hspace=0.16, wspace=0.08))
    for ax, kappa, Kmax in zip(axes.flat, (2, 10, 100, 1000), (6, 7, 7, 7)):
        d = np.load(DATA / f'phase_{kappa:g}.npz')
        aL, be, lab = d['aL'], d['beta'], d['lab']
        ax.imshow(colorize(lab), origin='lower', extent=(0, aL[-1], be[0], be[-1]), aspect='auto',
                  interpolation='nearest')
        overlay(ax, kappa)
        label_tongues(ax, aL, be, lab, Kmax)
        ax.set_xlim(0, 4)
        ax.set_ylim(0, 0.995)
        ax.text(0.97, 0.06, r'$\kappa=%g$' % kappa, transform=ax.transAxes, fontsize=8, va='bottom', ha='right')
        for side in ('top', 'right'):
            ax.spines[side].set_visible(False)
    for ax in axes[1]:
        ax.set_xlabel(r'$\alpha L$')
    for ax in axes[:, 0]:
        ax.set_ylabel(r'$\beta$')
    handles = [Patch(fc=C_Q, label=r'$\mathcal{Q}$: $\rho=\rho_{\rm q}$'),
               Patch(fc=C_UD, label=r'$\Omega_1$, between damping curves'),
               Patch(fc=C_REST, label=r'$\Omega_1$, rest'),
               Patch(fc=C_T1, label=r'tongues of $K$-cycles'),
               Patch(fc=C_GAP, label='gaps'),
               Patch(fc=C_DIV, label=r'$\alpha L\geq2(1+\beta)$'),
               plt.Line2D([], [], color=INK, lw=0.6, ls=(0, (3, 1.5)), label=r'$\alpha\mu=(1-\sqrt{\beta})^2$'),
               plt.Line2D([], [], color=INK, lw=0.6, ls=(0, (1, 1.2)), label=r'$\alpha L=(1+\sqrt{\beta})^2$'),
               plt.Line2D([], [], color=INK, marker='*', lw=0, ms=6, label="Polyak's tuning")]
    fig.legend(handles=handles, loc='lower center', ncol=3, bbox_to_anchor=(0.5, -0.085), fontsize=6.5,
               handlelength=1.6, columnspacing=1.0)
    save(fig, 'fig_phase')


if __name__ == '__main__':
    main()
