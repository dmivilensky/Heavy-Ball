"""Run every verification of the paper from source and record the results.

Usage:  python verify_all.py [--jobs 4] [--figures]

The Part I programs are copied into a fresh temporary directory and executed there, so that every
intermediate algebraic object (pickles, continuation boxes) is regenerated. Printed exact Boolean
checks are promoted to assertions by a small harness. The Part II programs and the spiral programs
are executed in place. With --figures the figure scripts are run afterwards (figures/make_all.py).
A failed identity, sign check, interval inclusion or branch connection gives a nonzero exit status.
Results are written to verification/manifest.json and verification/<script>.log.
"""
from pathlib import Path
import argparse
import concurrent.futures
import hashlib
import importlib.metadata
import json
import shutil
import subprocess
import sys
import tempfile
import time

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'verification'

PART1 = {
    'spiral': ['base/base_check.py'],
    'dual': ['upper/upper_sym.py'],
    'lower': ['lower/a3_identities.py', 'lower/a3_regions.py'],
    'struct': ['struct/' + s + '.py' for s in (
        'a7_involution', 'a7_thresholds', 'a7_minpolys', 'a7_cycles', 'a7_threshold_kc',
        'a8_kkt', 'a8_elim', 'a8_onewave', 'a8_curves', 'a8_regular', 'a9_asym', 'a9_kappa1')],
    'regimes': ['regimes/verify_regimes.py'],
    'analytic': ['analytic/p5_critical.py', 'analytic/p2_rational.py'],
}
DIRECT = ['part1/spiral/spiral_kernel.py', 'part1/spiral/spiral_subopt.py',
          'part2/identities.py', 'part2/convexity.py', 'part2/fan_identities.py', 'part2/gap_examples.py']

HARNESS = r'''
import builtins, runpy, sys
from sympy.logic.boolalg import BooleanFalse
original_print = builtins.print
def checked_print(*args, **kwargs):
    def check(x):
        if x is False or isinstance(x, BooleanFalse):
            raise AssertionError('A printed exact Boolean check is false: ' + repr(args))
        if isinstance(x, (list, tuple)):
            for y in x: check(y)
    for x in args: check(x)
    original_print(*args, **kwargs)
builtins.print = checked_print
runpy.run_path(sys.argv[1], run_name='__main__')
'''


def run(cmd, cwd, log):
    start = time.monotonic()
    with log.open('w') as f:
        proc = subprocess.run(cmd, cwd=cwd, stdout=f, stderr=subprocess.STDOUT)
    text = log.read_text()
    passed = proc.returncode == 0 and 'FAIL' not in text.replace('FAILED_ALLOWED', '')
    return passed, proc.returncode, round(time.monotonic() - start, 2)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--jobs', type=int, default=4)
    parser.add_argument('--figures', action='store_true')
    args = parser.parse_args()
    if not __debug__:
        raise RuntimeError('Do not run the verification with python -O')
    OUT.mkdir(exist_ok=True)
    sources = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
               for p in sorted(ROOT.glob('part*/**/*.py'))}
    versions = {}
    for name in ('sympy', 'mpmath', 'python-flint', 'numpy', 'scipy', 'matplotlib'):
        try:
            versions[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            versions[name] = None
    results = []
    with tempfile.TemporaryDirectory(prefix='hb_two_answers_') as tmp:
        tmp = Path(tmp)
        shutil.copytree(ROOT / 'part1' / 'hb_an', tmp / 'hb_an',
                        ignore=shutil.ignore_patterns('*.pkl', '*.out', '__pycache__'))

        def run_group(item):
            name, scripts = item
            out = []
            for rel in scripts:
                p = tmp / 'hb_an' / rel
                passed, code, sec = run([sys.executable, '-u', '-c', HARNESS, str(p)], p.parent, OUT / (p.stem + '.log'))
                out.append(dict(script='part1/hb_an/' + rel, passed=passed, exit_code=code, seconds=sec))
                print(('PASS  ' if passed else 'FAIL  ') + 'part1/hb_an/' + rel, f'({sec} s)', flush=True)
                if p.stem == 'verify_regimes' and passed:
                    for boxes in p.parent.glob('*_boxes.json'):
                        shutil.copy2(boxes, OUT / boxes.name)
                if not passed:
                    break
            return out

        def run_direct(rel):
            p = ROOT / rel
            passed, code, sec = run([sys.executable, '-u', str(p)], p.parent, OUT / (p.stem + '.log'))
            print(('PASS  ' if passed else 'FAIL  ') + rel, f'({sec} s)', flush=True)
            return [dict(script=rel, passed=passed, exit_code=code, seconds=sec)]

        with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
            futures = [pool.submit(run_group, it) for it in PART1.items()]
            futures += [pool.submit(run_direct, rel) for rel in DIRECT]
            for f in futures:
                results += f.result()
    expected = sum(map(len, PART1.values())) + len(DIRECT)
    passed = len(results) == expected and all(r['passed'] for r in results)
    if args.figures:
        ok, code, sec = run([sys.executable, '-u', str(ROOT / 'figures' / 'make_all.py')], ROOT / 'figures',
                            OUT / 'figures.log')
        results.append(dict(script='figures/make_all.py', passed=ok, exit_code=code, seconds=sec))
        print(('PASS  ' if ok else 'FAIL  ') + 'figures/make_all.py', f'({sec} s)', flush=True)
        passed = passed and ok
    manifest = dict(passed=passed, python=sys.version, versions=versions, sources_sha256=sources, results=results)
    (OUT / 'manifest.json').write_text(json.dumps(manifest, indent=2))
    print('ALL CHECKS PASSED' if passed else 'VERIFICATION FAILED')
    return 0 if passed else 1


if __name__ == '__main__':
    sys.exit(main())
