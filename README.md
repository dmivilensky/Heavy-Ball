# Two Answers on Heavy-Ball Dynamics: verification and figures

This package contains every computation used in the paper: the computer-assisted parts of the proofs of
Part I (the exact rate of Polyak's heavy ball for all condition numbers), the identity and certificate
checks of Part II (general tunings, the language of levels, the first level), and the programs that draw
all eight figures. Statement, equation, table and figure numbers refer to the paper; for the names used
inside the programs see "Numbering in the code" below.

## Contents

```
Heavy-Ball/
  verify_all.py              runs every verification and writes verification/manifest.json
  requirements.txt           pinned versions
  part1/
    hb_an/base/              Proposition 4.1, Lemma B.1
    hb_an/upper/             Lemma 4.2, the spiral certificate (4.14), Bernstein library
    hb_an/lower/             the uniform lower bound (2.3), Lemma C.1, Table 5
    hb_an/struct/            threshold, switching points, curves of P1 and P2, regularity, expansions
    hb_an/regimes/           polynomial systems and certified continuation of P1, P2, P5 (Table 2)
    hb_an/analytic/          closed forms: the arc P2 and the central family at kappa_c
    spiral/                  Lemma B.2, Proposition B.4
    boxes/                   certified continuation boxes (output of hb_an/regimes)
    rate.py                  floating-point evaluation of the rate prescription (4.13), used by the figures
  part2/
    hbgeom.py                geometry of general tunings (angle function, regions, linear programs)
    identities.py            polynomial identities of Sections 6-7 and Appendix D
    convexity.py             Theorem 8.3
    fan_identities.py        Theorems 8.1, 8.2, 9.1, Lemmas 9.3, E.3, equation (E.1)
    gap_examples.py          Propositions 9.4 and 9.5
  figures/                   one program per figure, the data programs, make_all.py, style.py; PDFs in figures/pdf/
  data/                      figure data (.npz) and PNG previews of the figures
  verification/              logs of the last complete run and manifest.json
```

## Running

```
pip install -r requirements.txt
python verify_all.py --jobs 6            # all verifications, about 12 minutes with 6 workers
python verify_all.py --jobs 6 --figures  # the same, followed by figures/make_all.py
python figures/make_all.py               # figures from the cached data in data/
python figures/make_all.py --data        # recompute all figure data first (several hours on one core)
```

`verify_all.py` copies `part1/hb_an` into a fresh temporary directory and runs every Part I program
there, so every intermediate object (pickled eliminants, continuation boxes) is regenerated from the
polynomial systems. The groups `spiral`, `dual`, `lower`, `struct`, `regimes`, `analytic` run in
parallel; inside a group the programs run in the listed order, because later programs read objects
written by earlier ones. The Part II programs and the two spiral programs run in place. Each program
writes `verification/<name>.log`; the run writes `verification/manifest.json` with the Python and package
versions, the SHA-256 hash of every source file, and the exit status and running time of every program.
The final line is `ALL CHECKS PASSED` exactly when all 25 programs pass. The figure programs write PDF
files to `figures/pdf/` (or to the folder given by the environment variable `HB_FIGURES`, for example
`HB_FIGURES=../paper/figures`) and PNG previews to `data/`.

Python 3.13 with the versions of `requirements.txt` was used for the run recorded in `verification/`.

## Notation in the code

The code uses ASCII names for the objects of the paper.

| Code | Paper |
|---|---|
| `X = e^t`, `Y = e^v`, `G = cosh t`, `u = sinh v`, `c = cos w` | the same; `rho = e^t` in Part I |
| `l_m`, `u_m` (raw slacks); `L(m)`, `U(m)` (normalized) | the slacks of Definition 3.3; the functionals of equation (4.3) |
| `W`, `W_s` in the spiral programs | the weight W(G) of the spiral certificate (4.14) |
| `Q` in `branch_systems.central` | 15 times the dual polynomial Q of (4.4) at h = 0 (checked in `p5_critical.py`) |
| `A, B` in the continuation of P1 and P5 | the blow-up coordinates of (B.4) and (B.6) |
| `Q` in `hbgeom.py` | the quadratic Xi of Lemma 6.2 |
| `W` in `convexity.py` | the polynomial written in sans-serif type in Theorem 8.3 |
| `lam` in `lower/` | the wave multiplier z = rho e^{iw} of (C.1) |

## Methods and their correctness

Every statement of the paper that depends on a computation is decided by one of the following methods.
Floating-point arithmetic is used only to propose centers, preconditioners and starting values, and in
the figure programs.

M1. Exact identities. An identity between rational functions with rational coefficients is checked by
expanding the numerator of the difference to the zero polynomial in `sympy` with rational arithmetic.
Identities in an algebraic extension, such as Q(y0) with q5(y0) = 0, Q(rho)[q]/(F) or Q(G)[u]/(P2), are
checked by reduction modulo the defining polynomial, with denominators inverted by the extended
Euclidean algorithm.

Lemma M1. If the defining polynomial p is irreducible over the base field, then Q[y]/(p) is a field,
reduction modulo p is a ring homomorphism onto it, and an element whose reduced representative is the
zero polynomial vanishes at every root of p. Hence an identity that reduces to zero holds at the root
used in the paper.

M2. Sturm sequences. The number of distinct real roots of a squarefree polynomial in (a, b] equals
V(a) - V(b), where V(x) is the number of sign changes in its Sturm sequence at x. All sequences are
computed in exact rational arithmetic.

Lemma M2. If a polynomial P with rational coefficients has no real root in [a, b] and P(x0) > 0 for one
x0 in [a, b], then P > 0 on [a, b]. For an element A(rho) + B(rho) q of Q(rho)[q]/(F), positivity on an
interval follows from the absence of real roots of the norm polynomial N(rho, q+) N(rho, q-) and of the
denominator, together with the sign at one point and the identification of the physical root q+
(`lower/exactlib.py`).

M3. Taylor dominance. A polynomial R with rational coefficients satisfies, for |zeta - zc| <= delta,
|R(zeta) - R(zc)| <= sum_{k >= 1} |R^(k)(zc)/k!| delta^k by Taylor's formula.

Lemma M3 (Lemma B.11). If |R(zc)| > sum_{k >= 1} |R^(k)(zc)/k!| delta^k, then R has the constant sign
sgn R(zc) on [zc - delta, zc + delta]. The comparison is a finite comparison of rational numbers
(`analytic/cert.py`).

M4. Bernstein positivity. On a box prod [a_i, b_i] a polynomial of multidegree (d_1, ..., d_n) is a
combination of the tensor Bernstein polynomials, which are nonnegative and sum to one.

Lemma M4. If every Bernstein coefficient is positive, the polynomial is positive on the box. If the
coefficients on every box of a finite subdivision are positive, the polynomial is positive on the union
(`upper/bernstein.py`, exact rational arithmetic).

M5. Interval arithmetic. Arb (`python-flint`) and `mpmath.iv` evaluate every operation with outward
rounding.

Lemma M5. For every point of the input intervals, the exact value of an expression lies in the interval
computed for it. Hence a strict inequality between computed intervals implies the same inequality for
all points of the inputs.

M6. Parametric Krawczyk test. For a system F(s, Z) = 0, a parameter interval I, a box X with center Z0
and a matrix R, the program verifies, in interval arithmetic,
K = Z0 - R F(I, Z0) + (Id - R F_Z(I, X))(X - Z0) contained in the interior of X, and a weighted infinity
norm of Id - R F_Z(I, X) smaller than one, with F(I, Z0) enclosed by F(s0, Z0) + (I - s0) F_s(I, Z0)
(equation (B.7), `regimes/interval_core.py`).

Lemma M6. Under these conditions, for every s in I the map Z -> Z - R F(s, Z) maps X into itself and is
a contraction, so F(s, .) has exactly one zero in X; R and F_Z(s, Z) are nonsingular there, and the zero
depends analytically on s by the implicit function theorem. If the enclosure of the zero at the common
endpoint of two adjacent slabs, computed from the left slab, lies in the box of the right slab, the two
zeros coincide by uniqueness, and the zeros on all slabs form one analytic branch. Derivatives along the
branch are enclosed by Z' = (Id - M)^{-1}(-R F_s), M = Id - R F_Z, with the Neumann remainder bounded by
|M| |R F_s|/(1 - |M|).

M7. Formal power series. Identities between rational functions of (s, r) are checked as identities in
Q[s][[r]], and nonnegativity of coefficients is checked on the closed forms of the factors
(`spiral/spiral_kernel.py`), with a direct cross-check of F_n for 3 <= n <= 200.

M8. The harness. `verify_all.py` runs every Part I program under a wrapper that turns every printed
`False` (Python or `sympy`) into an exception. Every failed assertion and every printed exact Boolean that is false makes the program exit with a nonzero status, and a log containing `FAIL` marks the program as failed. Values that are only printed are informational and are not checked by the harness.

Lemma M8. If `verify_all.py` prints `ALL CHECKS PASSED`, then all 25 programs ran to completion, every
assertion and every printed exact Boolean check was true, and the continuation boxes in
`verification/*_boxes.json` were produced in this run.

M9. Figure data. The Part II figures classify grid points by the closed forms of Corollary 7.3 and Theorem 7.1 and by linear programs over frequency grids (`part2/hbgeom.py`, `scipy.optimize.linprog`
with HiGHS): `periodic_feasible` solves the program of Theorem 9.1(a) for a given period, and
`certificate_margin` searches a fan-lag certificate of Theorem 9.2 on a grid of 6000 frequencies
that contains all 2 pi j/k with k < 80 (12000 frequencies in `gap_refine.py`). Part I curves are evaluated by Newton's method on the polynomial systems, started at the
centers of the certified boxes (`part1/rate.py`).

## Correspondence between the paper and the code

Part I. Paths are relative to `part1/hb_an/` unless they start with `part1/`.

| Statement in the paper | Program | Methods |
|---|---|---|
| Proposition 4.1 (spiral), equations (2.1), (B.1); kappa_c is the root of F(1, q) | `base/base_check.py` | M1 |
| Lemma B.1 (lag identities on the spiral), lags 1 to 12 and E_1 = E_2 = 0 | `base/base_check.py` | M1 |
| Lemma B.2 (positivity of the spiral kernel), equation (4.2) | `part1/spiral/spiral_kernel.py` | M1, M7 |
| Lemma 4.2 (dual factorization), equations (4.4), (4.5); the spiral certificate (4.14), W(G), G_1, G_W | `upper/upper_sym.py` | M1, M2 |
| Proposition B.4 (the spiral outside its regime): the vector y and the derivative identity | `part1/spiral/spiral_subopt.py`, `upper/upper_sym.py` (S3) | M1 |
| Lemma A.6(i); the functionals (4.3) are even in t; their values at c = cosh t and at c = +-1 (proof of Lemma A.7) | `struct/a7_involution.py` | M1 |
| Proposition B.5 (threshold strip), equation (B.3), the cubic and 4 A1 A3 - A2^2 | `struct/a7_threshold_kc.py` | M1 |
| Proposition B.5 at t = 0 and the three-cycle; Theorem 2.2 | `struct/a7_cycles.py`, `struct/a7_threshold_kc.py` | M1 |
| Table 1: kappa_W, kappa_1, kappa_2, kappa_c and their minimal polynomials (kappa_a, kappa_b are evaluated from G_a in floating point by `figures/fig_rates.py`) | `struct/a7_thresholds.py`, `struct/a7_minpolys.py`, `upper/upper_sym.py` (S6) | M1, M2 |
| Theorem 4.3: the optimality systems (4.6), (4.7), (4.10) and the curves P1, P2 (4.8) | `struct/a8_kkt.py`, `struct/a8_elim.py`, `struct/a8_onewave.py` | M1 |
| Lemma B.10(a): P2 restricted to the spiral is -(32G^3 - 36G^2 + 4G - 1)^2; P1 meets the spiral at G_W | `struct/a8_curves.py` | M1 |
| Proposition B.6 (family P1), equations (B.4), (B.5), Table 2 | `regimes/branch_systems.py`, `regimes/continue_certify.py p1`, `regimes/branch_checks.py` | M5, M6 |
| Proposition B.6: finite-lag slacks (B.8), 1 <= m <= 12, including the closed forms of l_1, v_2, v_3 | `regimes/p1_slacks.py`, `upper/bernstein.py` | M1, M4 |
| Lemma B.7, Definition B.8, Lemma B.9, Lemma B.10, Theorem B.12, Tables 3 and 4 (the arc P2) | `analytic/p2_rational.py` | M1, M2, M3 |
| Lemma B.11 (sign certificates), 48 polynomials | `analytic/cert.py`, `analytic/p2_rational.py` | M3 |
| The arc P2 in the field of its curve; interval cross-check over 4000 slabs | `regimes/p2_field.py`, `regimes/continue_certify.py p2` | M1, M5, M6 |
| Theorem B.13 (closed form at kappa_c), Proposition B.14 (dual nonnegativity), Proposition B.15 (first-order slacks), the limits in the proof of Theorem B.17 | `analytic/p5_critical.py` | M1, M3 |
| Proposition B.18 (the central family), equations (B.6), (B.10), Table 2 | `regimes/continue_certify.py central`, `regimes/branch_checks.py` | M5, M6 |
| Endpoint identifications at G_W, G_1, G_a and at the critical point h = 0 (Theorem 4.3, Lemma B.10, Theorem B.13, Proposition B.18) | `regimes/endpoints.py` | M1, M2, M5 |
| Theorem 2.3 (regularity): P_u and d/dt[t - V(t)] do not vanish on the arcs | `struct/a8_regular.py` | M2 |
| Theorem 2.3 (expansions): kappa to infinity, kappa_c (equation (2.4)), q to 0 | `struct/a9_asym.py`, `struct/a9_kappa1.py` | M1 |
| Equation (2.3) (uniform lower bound): Lemma C.1, equations (C.1)-(C.3) | `lower/a3_identities.py` | M1 |
| Equation (2.3): the sign certificates of Table 5 | `lower/a3_regions.py`, `lower/exactlib.py`; `lower/gen_table.py` prints Table 5 | M1, M2 |

Part II. Paths are relative to `part2/`.

| Statement in the paper | Program | Methods |
|---|---|---|
| Lemma 3.1 (normalization), equation (6.3) | `identities.py` | M1 |
| Lemma 6.2 (angle function): g sin w = 2s Xi(cos w), Xi(+-1) | `identities.py` | M1 |
| Lemma D.2 (cubic), equations (6.6), (7.1) (the polynomial of [5, Theorem 3.5]) | `identities.py` | M1 |
| Theorem 7.2(b), Lemma D.8, Lemma D.9 (K_mu, K_L), Corollary 7.3(b), Corollary 7.4(a), (b) | `identities.py` | M1 |
| Lemma D.6(a): the two forms of the symbol | `identities.py` | random-point evaluation |
| Theorem 8.3 (convexity of the angle function) | `convexity.py` | M1 |
| Theorem 8.2(a), Lemma E.3 | `fan_identities.py` | exact integer arithmetic, k <= 80, K, p <= 40 |
| Theorem 8.1(b), Theorem 8.2(b) | `fan_identities.py` | floating-point checks at 4000 random points and 37 phases per fraction |
| Equation (E.1) (Lambda_n as a polynomial in c) | `fan_identities.py` | M1 |
| Theorem 9.1(b), (c) and Lemma 9.3 | `fan_identities.py` | consistency checks in floating point |
| Proposition 9.5 (convergence in a gap) | `gap_examples.py` | M1, M2 |
| Proposition 9.4 (a two-harmonic thirteen-cycle) | `gap_examples.py` | M5 (200-bit `mpmath.iv`) |

Figures. Paths are relative to `figures/`.

| Figure | Program | Data |
|---|---|---|
| Figure 1 (the exact factor) | `fig_rates.py` | `part1/rate.py`, `data/rates.npz` |
| Figure 2 (the spectral boundary) | `fig_boundary.py` | `part1/rate.py`, `data/boundary_t.npz` |
| Figure 3 (worst-case spectra and trajectories) | `fig_spectra.py` | `part1/rate.py` |
| Figure 4 (the angle function) | `fig_thomae.py` | `gap_data.boundary` |
| Figure 5 (phase diagrams) | `fig_phase.py` | `phase_data.py`, `data/phase_*.npz` |
| Figure 6 (the regions at kappa = 100) | `fig_regions.py` | `data/phase_100.npz`, `data/gap_100.npz` |
| Figure 7 (lattice coordinates) | `fig_lattice.py` | `part2/hbgeom.py` |
| Figure 8 (the first level) | `fig_level1.py` | `gap_data.py`, `gap_refine.py`, `data/gap_100*_refined.npz` |

## Numbering in the code

Docstrings, comments and printed messages name statements by the labels of the LaTeX source or by the
numbers of an earlier draft. The table gives the numbers of the paper.

| Name in the code | Paper |
|---|---|
| Theorem 2, Theorem 3 | Theorem 2.2, Theorem 2.3 |
| Lemma normalization, Theorem spectral, Theorem classification | Lemma 3.1, Theorem 3.6, Theorem 4.3 |
| Proposition 10, Lemma 11, Theorem 12 | Proposition 4.1, Lemma 4.2, Theorem 4.3 |
| Lemma 36, Lemma 37, Lemma startup | Lemma A.6, Lemma A.7, Lemma A.5 |
| Lemma 38, Lemma 39, Proposition 41, Proposition 42, Proposition 43 | Lemma B.1, Lemma B.2, Proposition B.4, Proposition B.5, Proposition B.6 |
| Lemma 44 (Lemma curve), Definition 45, Lemma 46 (Lemma p2ids), Lemma 47 (Lemma p2ends), Lemma 48 | Lemma B.7, Definition B.8, Lemma B.9, Lemma B.10, Lemma B.11 |
| Theorem 49 (Theorem p2), Theorem 50, Proposition 51 (Proposition p5h0dual), Proposition 52 | Theorem B.12, Theorem B.13, Proposition B.14, Proposition B.15 |
| Theorem 54, Proposition 55, Lemma 56 | Theorem B.17, Proposition B.18, Lemma C.1 |
| Lemma angle, Theorem A, Theorem B, Corollary C, Corollary D | Lemma 6.2, Theorem 7.1, Theorem 7.2, Corollary 7.3, Corollary 7.4 |
| Lemma cubic, Lemma symbols, Lemma near1 | Lemma D.2, Lemma D.6, Lemma D.8 |
| Theorem lattice, Theorem farey, Theorem convex | Theorem 8.1, Theorem 8.2, Theorem 8.3 |
| Theorem fan-div, Theorem fan-conv, Lemma pin | Theorem 9.1, Theorem 9.2, Lemma 9.3 |
| prop:13, prop:conv-gap | Proposition 9.4, Proposition 9.5 |
| Proposition periodic, Lemma fan | Proposition E.2, Lemma E.3 |
| Table p2lags; GTD Theorem 3.5 | Table 4; [5, Theorem 3.5] |
| equations spiralpoly, spirallower, critical, spiralsimple | (2.1), (2.3), (2.4), (B.1) |
| equations dualpoly, spiralcertificate | (4.4), (4.14) |
| equations p1system, p2system, p2poly, p5system | (4.6), (4.7), (4.8), (4.10) |
| equations p1blowup, p1quartic, blowup, krawczyk, p1slackpoly, centralbounds | (B.4), (B.5), (B.6), (B.7), (B.8), (B.10) |
| equations explicitdata, startupslack, startupidentities | (C.1), (C.2), (C.3) |
| eq. factor, equation (E) and eq. E, eq. gtd, equation Lambdapoly | (6.3), the second form of (6.6), (7.1), (E.1) |

## Outputs of the recorded run

`verification/` contains the log of each program, the certified continuation boxes
`p1_boxes.json`, `p2_boxes.json`, `central_boxes.json`, and `manifest.json`. `part1/boxes/` holds the
same boxes; `part1/rate.py` reads their centers as starting values.
