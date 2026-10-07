"""Shared plotting style for all figures of the paper.

Figures are drawn at the text width of the SIAM template (370.4 pt = 5.125 in) with the
Latin Modern Roman font, so that labels appear at the size of the body text. The categorical
colors are the validated reference palette (blue, orange, aqua, violet for the four spectral
families; their pairwise CVD separation passes the all-pairs test).
"""
import sys
import os
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'data'
OUT = Path(os.environ.get('HB_FIGURES', ROOT / 'figures' / 'pdf'))
sys.path.insert(0, str(ROOT / 'part1'))
sys.path.insert(0, str(ROOT / 'part2'))

TEXTWIDTH = 370.38374 / 72.27

for f in Path('/usr/share/texmf/fonts/opentype/public/lm').glob('lmroman10-*.otf'):
    font_manager.fontManager.addfont(str(f))

plt.rcParams.update({
    'font.family': 'LM Roman 10',
    'mathtext.fontset': 'cm',
    'font.size': 8,
    'axes.titlesize': 8,
    'axes.labelsize': 8,
    'xtick.labelsize': 7,
    'ytick.labelsize': 7,
    'legend.fontsize': 7,
    'axes.linewidth': 0.5,
    'xtick.major.width': 0.5,
    'ytick.major.width': 0.5,
    'xtick.minor.width': 0.4,
    'ytick.minor.width': 0.4,
    'xtick.major.size': 2.5,
    'ytick.major.size': 2.5,
    'xtick.minor.size': 1.5,
    'ytick.minor.size': 1.5,
    'xtick.direction': 'out',
    'ytick.direction': 'out',
    'axes.edgecolor': '#3d3c39',
    'axes.labelcolor': '#0b0b0b',
    'xtick.color': '#3d3c39',
    'ytick.color': '#3d3c39',
    'lines.linewidth': 1.2,
    'lines.solid_capstyle': 'round',
    'legend.frameon': False,
    'savefig.dpi': 600,
    'pdf.fonttype': 42,
})

BLUE = '#2a78d6'
ORANGE = '#eb6834'
AQUA = '#1baf7a'
YELLOW = '#eda100'
MAGENTA = '#e87ba4'
GREEN = '#008300'
VIOLET = '#4a3aa7'
RED = '#e34948'
INK = '#0b0b0b'
INK2 = '#52514e'
MUTED = '#8a8984'
RULE = '#d4d3cf'
GRID = '#ecebe7'

FAMILY = {'P1': ORANGE, 'S': BLUE, 'P2': AQUA, 'P5': VIOLET}
FAMILY_TEX = {'P1': r'$P_1$', 'S': r'$S$', 'P2': r'$P_2$', 'P5': r'$P_5$'}


def tint(hex_color, a):
    """Mix a color with white; a=1 gives the color, a=0 white."""
    h = hex_color.lstrip('#')
    r, g, b = (int(h[i:i + 2], 16) for i in (0, 2, 4))
    mix = lambda v: int(round(255 - a * (255 - v)))
    return '#%02x%02x%02x' % (mix(r), mix(g), mix(b))


def clean(ax, grid=False):
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    if grid:
        ax.grid(True, color=GRID, lw=0.5, zorder=0)
        ax.set_axisbelow(True)


def save(fig, name):
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / f'{name}.pdf', bbox_inches='tight', pad_inches=0.02)
    fig.savefig(DATA / f'{name}.png', bbox_inches='tight', pad_inches=0.02, dpi=200)
    plt.close(fig)
