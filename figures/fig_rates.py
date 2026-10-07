"""Figure fig_rates: the exact factor of Polyak's tuning (Theorems rate, threshold, classification).

Left: rho_HB(kappa) colored by the optimal family, the quadratic factor q, and kappa_c.
Right: the excess 100 (rho_HB / rho_sp - 1) of the exact factor over the special spiral.
rho_HB is evaluated by part1/rate.py (the scalar prescription of Theorem classification);
rho_sp is the positive root of F(rho, q) (equation spiralpoly).
"""
import math
import numpy as np
from style import plt, save, clean, FAMILY, FAMILY_TEX, INK, INK2, MUTED, RULE, TEXTWIDTH, DATA
import rate

PHI = (1 + 5 ** 0.5) / 2
KC = 9 + 4 * 5 ** 0.5


def q_of(k):
    s = math.sqrt(k)
    return (s - 1) / (s + 1)


def rho_sp(q):
    r = np.roots([1, -2 * q, 1 + 2 * q - q * q, -2 * q, -q * q])
    r = r[(abs(r.imag) < 1e-12) & (r.real > 0)].real
    return float(r.max())


def kappa_of_q(q):
    return ((1 + q) / (1 - q)) ** 2


def switching_points():
    out = {}
    for name, G in (('W', rate.GW), ('1', rate.G1), ('a', rate.GA)):
        u, _ = rate.boundary(G)
        f = math.sqrt(1 + u * u) - u
        out[name] = (kappa_of_q((G - math.sqrt(G * G - 1)) * f), kappa_of_q((G + math.sqrt(G * G - 1)) * f))
    return {'W': out['W'][0], '1': out['1'][0], 'a': out['a'][0], 'b': out['a'][1], '2': out['1'][1]}


def compute():
    cache = DATA / 'rates.npz'
    if cache.exists():
        d = np.load(cache, allow_pickle=True)
        return d['kap'], d['rho'], d['rsp'], d['fam'], d['sw'].item()
    sw = switching_points()
    kap = np.unique(np.r_[np.geomspace(1.02, 1e4, 900), np.linspace(1.02, 3.2, 200),
                          np.linspace(7.4, 9.2, 300), np.linspace(45, 80, 250),
                          [v * (1 + e) for v in sw.values() for e in (-1e-6, 1e-6)], [KC]])
    rho, fam, rsp = [], [], []
    for k in kap:
        r = rate.evaluate(float(k))
        rho.append(r['rho']); fam.append(r['family']); rsp.append(rho_sp(r['q']))
    np.savez(cache, kap=kap, rho=np.array(rho), rsp=np.array(rsp), fam=np.array(fam), sw=sw)
    return kap, np.array(rho), np.array(rsp), np.array(fam), sw


def runs(fam):
    out, start = [], 0
    for i in range(1, len(fam) + 1):
        if i == len(fam) or fam[i] != fam[start]:
            out.append((fam[start], start, i))
            start = i
    return out


def main():
    kap, rho, rsp, fam, sw = compute()
    fig = plt.figure(figsize=(TEXTWIDTH, 2.55))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.2, 1], hspace=0.55, wspace=0.34)
    a0 = fig.add_subplot(gs[:, 0])
    # left panel: the exact factor
    qv = np.array([q_of(k) for k in kap])
    a0.plot(kap, qv, color=MUTED, lw=0.8, ls=(0, (4, 2)), zorder=2)
    a0.text(2.5e3, q_of(2.5e3) - 0.06, r'$q$', color=INK2, ha='center', va='top')
    for f, i, j in runs(fam):
        sl = slice(max(i - 1, 0), j)
        a0.plot(kap[sl], rho[sl], color=FAMILY[f], lw=1.6, zorder=4)
    for key in ('W', '1', 'a', 'b', '2'):
        k0 = sw[key]
        a0.plot([k0], [np.interp(k0, kap, rho)], marker='o', ms=2.2, color=INK, zorder=5)
    a0.axhline(1, color=RULE, lw=0.6, zorder=1)
    a0.axvline(KC, color=INK2, lw=0.6, ls=(0, (1, 1.5)), zorder=1)
    a0.text(KC * 1.13, 0.05, r'$\kappa_c=9+4\sqrt{5}$', color=INK2, rotation=90, va='bottom', ha='left', fontsize=7)
    a0.axhline(1.71298, color=RULE, lw=0.5, ls=(0, (2, 2)), xmin=0.55, zorder=1)
    a0.text(1.1e4, 1.71298, r'$\rho_\ast$', color=INK2, va='center', ha='left', fontsize=7)
    a0.set_xscale('log')
    a0.set_xlim(1, 1e4)
    a0.set_ylim(0, 1.8)
    a0.set_xlabel(r'condition number $\kappa$')
    a0.set_ylabel(r'factor per step')
    a0.text(1.2e3, 1.53, r'$\rho_{\rm HB}$', color=INK, ha='center')
    lab = {'P1': (1.33, 0.40), 'P2': (6.0, 0.95), 'P5': (25.0, 1.20)}
    for f, (x, y) in lab.items():
        a0.text(x, y, FAMILY_TEX[f], color=FAMILY[f], ha='center', va='center', fontsize=8)
    a0.text(4.0, 0.73, FAMILY_TEX['S'], color=FAMILY['S'], ha='center', va='center', fontsize=8)
    a0.text(260, 1.36, FAMILY_TEX['S'], color=FAMILY['S'], ha='center', va='center', fontsize=8)
    a0.text(-0.2, 1.0, '(a)', transform=a0.transAxes, fontsize=8, va='bottom')
    clean(a0, grid=True)
    # right panels: the excess over the special spiral on its two windows
    ex = 100 * (rho / rsp - 1)
    wins = [(1.0, sw['W'] * 1.08, '(b)'), (sw['1'] / 1.2, sw['2'] * 1.2, '(c)')]
    for r, (lo, hi, tag) in enumerate(wins):
        ax = fig.add_subplot(gs[r, 1])
        m = (kap >= lo) & (kap <= hi)
        for f, i, j in runs(fam):
            idx = np.arange(max(i - 1, 0), j)
            idx = idx[m[idx]]
            if len(idx):
                ax.plot(kap[idx], ex[idx], color=FAMILY[f], lw=1.4, zorder=4)
        top = ex[m].max()
        ax.axhline(0, color=RULE, lw=0.6)
        ax.set_xscale('log')
        ax.set_xlim(lo, hi)
        ax.set_ylim(-0.25 * top, 1.18 * top)
        clean(ax, grid=True)
        ax.text(-0.2, 1.0, tag, transform=ax.transAxes, fontsize=8, va='bottom')
        if r == 0:
            ax.plot([sw['W']], [0], marker='o', ms=2.4, color=INK, zorder=5)
            ax.text(sw['W'], -0.08 * top, r'$\kappa_W$', ha='right', va='top', fontsize=7, color=INK2)
            ax.text(1.55, 0.55 * top, FAMILY_TEX['P1'], color=FAMILY['P1'], fontsize=8, ha='center')
            ax.set_xticks([1, 1.5, 2, 2.5])
            ax.set_xticklabels(['1', '1.5', '2', '2.5'])
        else:
            for key, name, ha in (('1', r'$\kappa_1$', 'right'), ('a', r'$\kappa_a$', 'left'),
                                  ('b', r'$\kappa_b$', 'right'), ('2', r'$\kappa_2$', 'left')):
                ax.plot([sw[key]], [0], marker='o', ms=2.4, color=INK, zorder=5)
                ax.text(sw[key] * (0.985 if ha == 'right' else 1.015), -0.06 * top, name, ha=ha, va='top', fontsize=7, color=INK2)
            ax.plot([KC], [0], marker='o', ms=2.4, mfc='white', mec=INK, mew=0.7, zorder=5)
            ax.text(KC, 0.12 * top, r'$\kappa_c$', ha='center', va='bottom', fontsize=7, color=INK2)
            ax.text(8.1, 0.62 * top, FAMILY_TEX['P2'], color=FAMILY['P2'], fontsize=8, ha='center')
            ax.text(62, 0.62 * top, FAMILY_TEX['P2'], color=FAMILY['P2'], fontsize=8, ha='center')
            ax.text(30, 1.0 * top, FAMILY_TEX['P5'], color=FAMILY['P5'], fontsize=8, ha='center')
            ax.set_xticks([8, 10, 20, 50, 70])
            ax.set_xticklabels(['8', '10', '20', '50', '70'])
            ax.set_xlabel(r'condition number $\kappa$')
        ax.xaxis.set_minor_formatter(plt.NullFormatter())
        ax.set_ylabel('excess (%)')
    save(fig, 'fig_rates')
    print('switching points', sw)


if __name__ == '__main__':
    main()
