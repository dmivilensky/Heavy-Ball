"""Floating-point evaluation of the rate prescription of Theorem classification.

Usage: python rate.py 2 8 10 17.94427190999916 60 100

For each condition number the scalar equation t - arsinh u(cosh t) = log q is solved,
with u(G) taken from the family that Theorem classification assigns to G. The families
P1, P2, P5 are evaluated by Newton's method on their polynomial systems (branch_systems.py),
started from the centers of the certified continuation boxes in boxes/. Used by the figures.
"""
import argparse
from functools import lru_cache
import json
import math
from pathlib import Path
import sys

import numpy as np
from scipy.optimize import brentq, root
import sympy as sp

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'hb_an' / 'regimes'))
import branch_systems

GA = 1.0235814270077666593
G1 = 1.0332964166603128676
GW = (1 + math.sqrt(3)) / 2


@lru_cache(None)
def system(family):
    par, variables, equations, *_ = getattr(branch_systems, family)()
    f = sp.lambdify((par, variables), equations, 'numpy', cse=True)
    j = sp.lambdify((par, variables), sp.Matrix(equations).jacobian(variables),
                    'numpy', cse=True)
    boxes = json.loads((ROOT / 'boxes' / f'{family}_boxes.json').read_text())
    nodes = np.array([(b['a'] + b['b']) / 2 for b in boxes])
    centers = np.array([b['center'] for b in boxes])
    return f, j, nodes, centers


def polynomial_solution(family, parameter):
    f, j, nodes, centers = system(family)
    i = int(np.argmin(np.abs(nodes - parameter)))
    solution = root(lambda x: f(parameter, x), centers[i],
                    jac=lambda x: j(parameter, x), tol=1e-11)
    residual = float(np.max(np.abs(f(parameter, solution.x))))
    if not np.all(np.isfinite(solution.x)) or residual > 1e-8:
        raise RuntimeError(f'{family} solve failed: residual={residual:g}')
    return solution.x


def boundary(g):
    """Return (u(G), family) using the proved branch selection."""
    if g >= GW:
        z = 1 / g
        a, _ = polynomial_solution('p1', z)
        return 1 - z / 2 + z * z * a, 'P1'
    if g >= G1:
        return 1 - 1 / (2 * g), 'S'
    if g >= GA:
        values = polynomial_solution('p2', g)
        return float(values[0]), 'P2'
    h = g - 1
    values = polynomial_solution('central', h)
    return 0.5 + h * values[0], 'P5'


def evaluate(kappa):
    if not math.isfinite(kappa) or kappa <= 1:
        raise ValueError('kappa must be finite and strictly greater than one')
    sq = math.sqrt(kappa)
    q = (kappa - 1) / (sq + 1) ** 2
    if not 0 < q < 1:
        raise ValueError('kappa is too large for this double-precision evaluator')
    logq = math.log(q)

    def equation(t):
        u, _ = boundary(math.cosh(t))
        return t - math.asinh(u) - logq

    t = brentq(equation, logq, logq + math.asinh(1),
               xtol=1e-14, rtol=1e-14)
    _, family = boundary(math.cosh(t))
    return dict(kappa=kappa, q=q, rho=math.exp(t), log_rate=t, family=family)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('kappa', type=float, nargs='+')
    args = parser.parse_args()
    for kappa in args.kappa:
        try:
            print(json.dumps(evaluate(kappa)))
        except (ValueError, RuntimeError) as exc:
            parser.exit(1, str(exc) + '\n')


if __name__ == '__main__':
    main()
