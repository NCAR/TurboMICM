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

## Family results

At rtol=1e-3 and atol=1e-18, the minimum accuracy across all six disclosed
scenarios and the sum of their C development times are:

| Method | Minimum overall digits | Minimum species digits | Sum of C times (ms) |
|---|---:|---:|---:|
| ROS2 | 2.737 | 1.749 | 23207 |
| ROS3 | 3.543 | 2.346 | 12771 |
| ROS4 | 3.855 | 2.721 | 16582 |
| RODAS3 | 2.892 | 2.052 | 4928 |
| RODAS4 | 3.000 | 2.184 | 7669 |
| Euler extrapolation 2 | 2.988 | 1.982 | 11109 |
| Euler extrapolation 3 | 3.225 | 2.137 | 9976 |
| RODAS3.1 | 2.799 | 1.689 | 3630 |

The initial seven-family scan used uncached retry forcing; RODAS3.1 was added
later with the exact forcing cache. Thus the table is a development overview,
not a matched-clock ranking. The within-family tolerance curves, operation
counts, and final repeated operating-point timings are the stronger evidence.

Second order was expensive for the required accuracy. Euler extrapolation was
functional and passed its independent order/prototype checks, but its multiple
factorizations per macro step made it uncompetitive here. RODAS3 was inexpensive
near two to three digits; its optimized RODAS3.1 coefficients improved the
measured curve further. Higher-order RODAS4 was useful for the precision point.
Equal local tolerances do not imply equal global accuracy, which is why the
final choices follow the measured curves and both species/overall constraints.

## Control, tolerance and reuse findings

- Warm step proposals were the largest clear control improvement. On the
  RODAS3 rtol=3e-4 comparison, cold starts required 313504 attempts across
  the six scenarios versus 175758 with warm starts. The first chemistry call
  still starts cautiously, and every rate-change boundary invalidates matrices.
- Maximum cell RMS avoids diluting a difficult cell when adding more cells.
  Group RMS was cheaper but reduced worst-species accuracy; a component
  maximum overspent work, reaching 4.80 digits overall at the same weights.
- Jacobian age 4 with deadband 0.2 reduced factorizations from 175758 to
  156360, but minimum species accuracy fell from 2.564 to 2.236 digits. The
  modest time saving did not justify consuming that margin in the recommendation.
- Exact retry-forcing reuse preserves the trajectory and saves one evaluation
  per rejection. Unlike stale Jacobians, it introduces no approximation.
- PI and damped H211b often produced more accuracy than requested. Follow-up
  tolerance scans did not place them ahead of the selected I-controller curve.
  Aggressive H211b b=1 increased rejections and exposed roundoff sensitivity;
  it is available for research but absent from the declared points. The gains
  reported for other atmospheric systems in the literature did not transfer
  automatically to this setup.
- RODAS3.1 at rtol=3e-4, atol=1e-16 reached 3.341 overall / 2.322 worst-species
  digits. Reducing atol to 1e-18 raised these minima to 3.705 / 2.734. This
  provides a measured speed/accuracy tradeoff without changing any reactions.

## Scaling and independent checks

At 10000 cells, RODAS3.1 rtol=3e-4/atol=1e-18 took 11222 ms with groups
of 32 and guided scheduling, compared with 11338 ms for groups of 16 and
11931 ms for groups of 64 under the same scheduling policy. Static/dynamic
differences were small. We select groups of 32 with guided minimum chunk 1
as an explicit runtime operating-point choice, not a hardware-specific kernel.
These nine development runs are in `scans/scale.csv`. Their scorer accuracy
covers the disclosed first 1000 cells only; the finite/nonnegative gate covers
all 10000 cells. This is not a claim of full 10000-cell reference accuracy.

Three independent 128-cell, 48-hour cases test high-latitude summer, polar
winter and equatorial conditions, including NOx multipliers 0.02 and 50.
Reference convergence is 10.076–10.349 overall digits, with at least 8.992
digits for every scored species. All five points pass their declared levels
on these cases. Generation settings and reference work/convergence are
preserved in `scans/stress_conditions.json`, `stress-reference-work.csv` and
`stress-reference-convergence.csv`; the point checks are in `stress.csv`.

## Final selection and stopping criterion

The five points in OPERATING_POINTS.md cover two accuracy budgets, extra margin
for each, and a precision option. The recommendation uses exact current
Jacobians, the strict absolute floor and the simple I controller. The final
timings use five sequential repetitions in a deterministic shuffled order,
with no other CPU-heavy task processes running. RESULTS.md reports their
medians rather than extrapolating from Python or from operation counts.

After the coefficient improvement, the controller/reuse/growth/safety scans
gave no compelling further improvement at the chosen accuracy margins. We stop
with these verified configurations rather than claiming a global optimum over
all possible integrators. BDF/NDF, multirate splitting and reduced precision
remain untested future directions. No performance or accuracy claim is made
for them.
