# Structure scan protocol

The numerical implementation is mechanism-generic. Every disclosed cell retains
all 209 integrated species and all 547 TS1 reactions. Every run is driven through
`bin/ts1_scenario`. Python coordinates runs and invokes the supplied scorer;
all times are measured inside C and exclude input, setup and CSV output.

## Controlled comparisons

- **Baseline:** original three-stage engine at default tolerances, then changing
  only rtol/atol to chemically meaningful local weights. This distinguishes a
  tolerance problem from an integrator-family problem.
- **Family and order:** ROS2, ROS3, ROS4, RODAS3, RODAS4, and order-2/order-3
  linearly implicit Euler extrapolation, at rtol 1e-3, 1e-4, 1e-5 with atol 1e-18.
  Within a tolerance column, all other experimental options are identical.
  Different families have different embedded estimators and costs; the resulting
  curves, not equal local tolerances, determine work at equal global accuracy.
- **Tolerance:** additional coarse and intermediate relative tolerances locate
  the 2- and 3-digit operating regions. Absolute tolerance is varied separately
  where needed; it is never equated with the scoring threshold.
- **Control and local norm:** one-at-a-time changes to I/PI control, local norm,
  warm starts, initial h, growth limit, safety and positivity rejection. This
  examines whether the chosen stage formula is being limited by its controller.
- **Reuse:** age 1,2,4,8 Jacobians with no deadband, then deadband sweeps at fixed
  age. Count Jacobian builds, LU factorizations and linear solves separately.
  This avoids calling fewer Jacobian evaluations a speedup when extra rejected
  steps erase the savings. Old-Jacobian Rosenbrock results are empirical; their
  exact-Jacobian order does not automatically carry over.

Precision and chemistry splitting remain fixed: all arithmetic is double and
all species are integrated together. Reducing precision would introduce a second
error source into the present comparisons; operator partitioning would require
new splitting-error studies. These are untested directions, not claimed
improvements. No mechanism-specific QSSA partition, reduction or surrogate is
used.

## Records and reproducibility

`scans/*.csv` contains one row per completed run, including failures. A row gives
its configuration, scenario, compiler, threads, supplied accuracy score, C time,
and operation counts. `scans/README.md` distinguishes diagnostic, corrected and
final timing records. Counts are per independent cell group, so compare the same
number and size of groups when attributing changes to a numerical choice.
Initial C/Python agreement tests are in `prototypes/check.py`; analytic order
and invariants are checked independently in `tests/test_new_methods.c`.

A scan time is an exploratory development measurement, not an evaluator median.
Final declared points are run repeatedly without other task processes competing
for CPU; their times and score minima appear in RESULTS.md. The supplied
high-accuracy references, scorer and read-only mechanism are never modified.

## Early findings

The default absolute tolerance is unsuitable for several day/night scenarios.
Replacing it makes the original method robust but does not remove its repeated
startup cost. Carrying a step proposal across chemistry calls reduces this
cost, while invalidating matrices when the photolysis changes preserves the
scenario's piecewise-constant forcing contract.

Second-order Rosenbrock is substantially more expensive at tight accuracy in
these TS1 runs. The family comparison tests whether the extra stages of ROS3,
ROS4 or stiffly accurate RODAS can buy larger accepted steps. Negativity-aware
rejection can materially affect this answer; accuracy and operation counts
must be compared together rather than judging by stage count alone.
