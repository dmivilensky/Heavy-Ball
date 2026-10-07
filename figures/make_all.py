"""Build all figures of the paper into figures/pdf/ (PDF; another folder can be set with the environment
variable HB_FIGURES) and data/ (PNG previews).

Usage: python make_all.py [--data]

The classification data of Part II are read from data/. They are recomputed when missing or when
--data is given: phase_data.py (about two minutes), gap_data.py for the first level at kappa = 100
(about an hour on one core), and gap_refine.py. Part I curves are computed by part1/rate.py and cached
in data/rates.npz and data/boundary_t.npz.
"""
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / 'data'


def run(*args):
    print('$', ' '.join(args), flush=True)
    subprocess.run([sys.executable, *args], cwd=HERE, check=True)


def main():
    fresh = '--data' in sys.argv
    DATA.mkdir(exist_ok=True)
    if fresh or not all((DATA / f'phase_{k}.npz').exists() for k in (2, 10, 100, 1000)):
        run('phase_data.py')
    if fresh or not (DATA / 'gap_100.npz').exists():
        run('gap_data.py', '100', '0.40', '0.93', '0.05', '900', '260')
    if fresh or not (DATA / 'gap_100_hi.npz').exists():
        run('gap_data.py', '100', '0.38', '0.66', '0.16', '700', '320', '_hi')
    if fresh or not (DATA / 'gap_100_refined.npz').exists() or not (DATA / 'gap_100_hi_refined.npz').exists():
        run('gap_refine.py')
    if fresh:
        for f in ('rates.npz', 'boundary_t.npz'):
            (DATA / f).unlink(missing_ok=True)
    for script in ('fig_rates.py', 'fig_boundary.py', 'fig_spectra.py', 'fig_thomae.py', 'fig_phase.py',
                   'fig_regions.py', 'fig_lattice.py', 'fig_level1.py'):
        run(script)
    print('ALL FIGURES BUILT')


if __name__ == '__main__':
    main()
