# Independent prototype checks

Run `.conda/bin/python prototypes/check.py` after building `bin/integrator_probe`.
No third-party solver is imported. The Python prototype independently implements
forcing, a negative analytic Jacobian, dense LU, triangular solves, all stage
recurrences, extrapolation, error norms, rejection handling and controllers.
Only the coefficient tables are shared as mathematical data printed by the C
probe. The retained tables also have Tier 1 tests; analytic linear and nonlinear
convergence tests check every experimental method's order independently.

The mechanism is 2A -> B, B -> 2A, B -> C, C -> 2A with rates 20, 2, 0.4 and
0.1. Two cells start from (1,0.1,0.01) and (0.5,0.2,0.03), and advance through
two intervals of 0.2 seconds. The invariant A + 2B + 2C is checked in C.
The prototype checks eight methods across I, PI, damped H211b, three norms,
warm/cold starts, and Jacobian/step reuse. The 48 strict comparisons require
identical operation counts and maximum absolute solution difference <3e-13.
The observed maximum with clang is 4.61e-15, consistent with accumulated
floating-point evaluation-order differences, not truncation error.

Forcing caching is additionally toggled in C in 64 comparisons: the resulting
trajectories must be bit-identical and the function-call saving must equal the
number of rejected attempts. No timing comes from Python.

## Aggressive-controller sensitivity

Sixteen additional comparisons use H211b with b=1 and k=2 or 1.7. These are
sensitivity diagnostics, not strict roundoff comparisons. Different floating-
point evaluation orders in independent dense Python and sparse C implementations
are amplified by this marginal controller into different adaptive grids. Some
runs have different attempt counts; even equal counts need not imply equal
step sizes. The maximum measured state difference is 5.40e-7 at rtol=1e-5.
`roundoff_sensitivity.csv` preserves the C/Python step counts and concentration
differences for every such case. Recreate the diagnostic file with
`prototypes/check.py --record-sensitivity` (ordinary tests do not rewrite it).
The diagnostic test bounds the absolute difference by 1e-6, while all eight
methods still must pass the separate strict roundoff comparisons above.

Dreger et al. (2025), equation (6), discuss b>1 for stability and the more
aggressive b=1 experiments; see REFERENCES.md. Our default b=4 avoids this
sensitivity in the small-case strict tests. Scenario experiments must establish
accuracy and robustness independently before an aggressive setting is declared.
This distinction is part of the evidence, not an assertion that adaptive grids
are invariant under roundoff for all controller parameters.

`tests/test_new_methods.c` also checks fixed-step convergence against exp(-3t)
and the exact Riccati solution of reversible dimerization. The latter uses
A' = -4A^2 - 0.5A + 0.5M and conserved M=A+2B. Both sets of tests require the
expected order under step halving and preserve the appropriate invariant.
