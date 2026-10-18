# Small-system prototypes

`integrators.py` implements Rosenbrock stages, error controllers, arbitrary-step BDF coefficients, polynomial prediction, modified Newton iteration, and dense Gaussian elimination in ordinary Python. It receives chemistry functions; the nonlinear three-species fixture and its analytical Jacobian are written independently in `tests/test_adaptive.py`. No C forcing, Jacobian, factorization, or integration routine is called by the prototype.

The inherited Rosenbrock coefficient tables are supplied as data through the public API. Their numerical values are separately covered by the retained Tier 1 regression tests. The prototypes test the new integration loops and controls, rather than re-verifying those tables.

Run:

```sh
make test-tier3 CC=clang
make test-tier3 CC=gcc-15
```

The tests compare all five Rosenbrock formulas with fixed/adaptive steps, reuse, the maximum norm, and H211b control. They compare BDF maximum orders 1–5 with the applicable control variants. Every small-case solution agrees within accumulated double-precision roundoff and every operation count matches. Reports are retained as `verification-clang.csv` and `verification-gcc15.csv`. BDF's polynomial extrapolation and residual subtraction amplify small arithmetic-order differences between dense Python and sparse C; the observed absolute differences are reported, not hidden by a nominal solver tolerance.

Additional C tests exercise three lanes with different amplitudes of a very stiff conservative transfer, finite/nonnegative output, equilibrium accuracy, weighted conservation on the nonlinear fixture, invalid parameters, cached photolysis vs full rate recalculation, partial cell groups, output snapshots, and scheduling independence. The stiff transfer check allows integration error: an arbitrarily large stable step is not expected to reproduce its exponentially small transient to roundoff.

Python runtime and its operation counts are not speed measurements. The prototypes need no third-party numerical integrator. Native C runs through `bin/ts1_scenario` supply every speed measurement in the study.
