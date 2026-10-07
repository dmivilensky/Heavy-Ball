"""Figure fig_boundary: the spectral boundary v = V(t) = arsinh u(cosh t) (Theorem classification).

The feasible set at the exponent t is |v| <= V(t) (Theorem spectral). The boundary is drawn in the
plane of the exponent t and the twist v and colored by family. For a condition number kappa the
actual point is (t, t - log q), a line of slope one; its last intersection with the boundary is
T = log rho_HB (marked). The right panel magnifies |t| <= 0.32, where the central family, the
two-frequency family and the spiral meet at the three-cycle (0, log phi).
"""
import math
import numpy as np
from style import plt, save, clean, FAMILY, FAMILY_TEX, INK, INK2, MUTED, RULE, TEXTWIDTH, DATA
import rate

KC = 9 + 4 * 5 ** 0.5
PHI = (1 + 5 ** 0.5) / 2


def q_of(k):
    s = math.sqrt(k)
    return (s - 1) / (s + 1)


def boundary_curve():
    cache = DATA / 'boundary_t.npz'
    if cache.exists():
        d = np.load(cache)
        return d['t'], d['v'], d['fam']
    tA, t1 = math.acosh(rate.GA), math.acosh(rate.G1)
    t = np.unique(np.r_[np.linspace(0, tA, 160), np.linspace(tA, t1, 80), np.linspace(t1, 0.4, 120),
                        np.linspace(0.4, 2.6, 300)])
    out = [rate.boundary(float(math.cosh(x))) for x in t]
    v = np.array([math.asinh(o[0]) for o in out])
    fam = np.array([o[1] for o in out])
    t = np.r_[-t[::-1], t[1:]]
    v = np.r_[v[::-1], v[1:]]
    fam = np.r_[fam[::-1], fam[1:]]
    np.savez(cache, t=t, v=v, fam=fam)
    return t, v, fam


def runs(fam):
    out, start = [], 0
    for i in range(1, len(fam) + 1):
        if i == len(fam) or fam[i] != fam[start]:
            out.append((fam[start], start, i))
            start = i
    return out


def draw_boundary(ax, t, v, fam, lw):
    for f, i, j in runs(fam):
        sl = slice(max(i - 1, 0), min(j + 1, len(t)))
        ax.plot(t[sl], v[sl], color=FAMILY[f], lw=lw, zorder=5, solid_capstyle='butt')


def draw_line(ax, k, tlo, thi, label=None, lx=None):
    q = q_of(k)
    tt = np.array([tlo, thi])
    ax.plot(tt, tt - math.log(q), color=MUTED, lw=0.6, zorder=2)
    T = rate.evaluate(k)['log_rate']
    ax.plot([T], [T - math.log(q)], 'o', ms=3.2, mfc='white', mec=INK, mew=0.8, zorder=7)
    return T


def main():
    t, v, fam = boundary_curve()
    fig, (a0, a1) = plt.subplots(1, 2, figsize=(TEXTWIDTH, 2.6), gridspec_kw=dict(width_ratios=[1.45, 1], wspace=0.28))
    # left: the whole boundary on -2.4 <= t <= 0.75
    kap_left = [1.3, 2.0, 4.0, KC, 100.0]
    m = (t >= -2.45) & (t <= 0.8)
    a0.fill_between(t[m], 0, v[m], color='#f3f2ee', zorder=0, lw=0)
    for k in kap_left:
        T = draw_line(a0, k, -2.45, 0.8)
    draw_boundary(a0, t[m], v[m], fam[m], 1.8)
    for k in (1.3, 2.0, 4.0):
        y = 0.93
        x = y + math.log(q_of(k))
        a0.text(x - 0.03, y, r'$\kappa=%g$' % k, color=INK2, fontsize=6.5, ha='right', va='center')
    for k, name in ((KC, r'$\kappa_c$'), (100.0, r'$\kappa=100$')):
        y = 0.07
        x = y + math.log(q_of(k))
        a0.text(x + 0.04, y, name, color=INK2, fontsize=6.5, ha='left', va='center')
    a0.text(-1.9, 0.81, FAMILY_TEX['P1'], color=FAMILY['P1'], ha='center')
    a0.text(-0.62, 0.66, FAMILY_TEX['S'], color=FAMILY['S'], ha='center')
    a0.text(0.62, 0.62, FAMILY_TEX['S'], color=FAMILY['S'], ha='center')
    a0.text(-1.2, 0.25, r'$|v|\leq V(t)$', color=INK2, ha='center', fontsize=7.5)
    xw = -math.acosh(rate.GW)
    yw = float(np.interp(xw, t, v))
    a0.plot([xw], [yw], 'o', ms=2.6, color=INK, zorder=8)
    a0.annotate(r'$-\mathrm{arcosh}\,G_W$', xy=(xw, yw), xytext=(-1.0, 0.44), fontsize=6.5, color=INK2,
                ha='center', va='center', arrowprops=dict(arrowstyle='-', lw=0.5, color=INK2, shrinkA=2, shrinkB=2))
    a0.add_patch(plt.Rectangle((-0.32, 0.47), 0.64, 0.045, fill=False, ec=INK2, lw=0.6, zorder=8))
    a0.set_xlim(-2.45, 0.8)
    a0.set_ylim(0, 0.98)
    a0.set_xlabel(r'exponent $t$')
    a0.set_ylabel(r'twist $v$')
    a0.text(-0.15, 1.0, '(a)', transform=a0.transAxes, va='bottom')
    clean(a0, grid=True)
    # right: the window of the three-cycle
    m = np.abs(t) <= 0.33
    a1.fill_between(t[m], 0, v[m], color='#f3f2ee', zorder=0, lw=0)
    for k in (8.3, KC, 30.0, 62.0):
        draw_line(a1, k, -0.33, 0.33)
    draw_boundary(a1, t[m], v[m], fam[m], 1.8)
    for s in (-1, 1):
        for g0, name in ((rate.GA, r'$\pm\arcosh G_a$'), (rate.G1, r'$\pm\arcosh G_1$')):
            x0 = s * math.acosh(g0)
            a1.plot([x0], [np.interp(x0, t, v)], 'o', ms=2.2, color=INK, zorder=8)
    a1.plot([0], [math.log(PHI)], 'o', ms=3.6, color=FAMILY['P5'], zorder=9)
    a1.annotate('three-cycle', xy=(0, math.log(PHI)), xytext=(0.075, 0.4745), fontsize=6.5, color=INK2,
                ha='left', va='center', arrowprops=dict(arrowstyle='-', lw=0.5, color=INK2))
    a1.text(-0.11, 0.4875, FAMILY_TEX['P5'], color=FAMILY['P5'], ha='center')
    a1.text(0.11, 0.4875, FAMILY_TEX['P5'], color=FAMILY['P5'], ha='center')
    a1.text(-0.238, 0.4985, FAMILY_TEX['P2'], color=FAMILY['P2'], ha='center')
    a1.text(0.238, 0.4985, FAMILY_TEX['P2'], color=FAMILY['P2'], ha='center')
    a1.text(-0.30, 0.5065, FAMILY_TEX['S'], color=FAMILY['S'], ha='center')
    a1.text(0.30, 0.5065, FAMILY_TEX['S'], color=FAMILY['S'], ha='center')
    for k, x in ((8.3, -0.300), (KC, -0.045), (30.0, 0.07), (62.0, 0.19)):
        q = q_of(k)
        y = 0.5135
        name = r'$\kappa_c$' if k == KC else r'$%g$' % k
        a1.text(y + math.log(q) - 0.008, y, name, color=INK2, fontsize=6.5, ha='right', va='center')
    a1.set_xlim(-0.33, 0.33)
    a1.set_ylim(0.47, 0.515)
    a1.set_xlabel(r'exponent $t$')
    a1.text(-0.2, 1.0, '(b)', transform=a1.transAxes, va='bottom')
    clean(a1, grid=True)
    save(fig, 'fig_boundary')


if __name__ == '__main__':
    main()
