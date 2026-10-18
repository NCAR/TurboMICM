# Numerical prototypes

`adaptive.py` is a dense, independent implementation of the adaptive Rosenbrock
and SDIRK algorithms. It accepts arbitrary forcing and Jacobian functions.
It uses Python/NumPy for arrays and explicit dense LU/triangular solves, with
actual factor reuse. No native kinetics, Jacobian, LU or integration code is
called to produce its answers. Python timings are never used.

Run the C/prototype comparison through:

```sh
make test-methods CC=clang
```

`tests/test_integrator.py` supplies an independently written nonlinear
three-species mechanism with two distinct cells and changing rate parameters.
Each of seven methods is tested with seven policy configurations and three
external intervals: 147 comparisons. Concentrations agree to accumulated
double-precision roundoff and every operation counter agrees exactly.

The five retained Tier 1 mathematical tableaux are shared data in this test
(already tested by the Tier 1 suite). The new RODAS3.1 tableau is independently
transcribed in Python; SDIRK's equations are implemented separately. Tests also
check seven stiff irreversible-chain solutions, conservation, nonnegativity,
formal convergence orders on a nonstiff chain, zero/invalid intervals and
invalid options. Native scenario contract tests cover real TS1 input/output,
forcing interpolation, boundary values, snapshot timing and thread independence.

The first draft used NumPy's dense solve and counted algorithmic LU work. The
committed prototype instead performs and reuses its own explicit factors, so
there is no hidden per-stage refactorization or counter qualification.
