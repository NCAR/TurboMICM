# Structure scans

Each scan changes one structural choice at a time and measures accuracy
(`scenarios/score.py`, digits overall and for the worst species), the
operation counts (internal steps, accepted, rejected, forcing evaluations,
Jacobian updates, factorizations, solves, factorization reuses), and, in C,
the time, on the disclosed scenarios. Every configuration tried is a row in
`scans/*.csv`; the Python prototypes supply accuracy and counts on small
cases (`prototypes/scenario.py`), and `bin/ts1_scenario` supplies all three
on the full scenarios. Times are from 1,000-cell scenarios at 10 threads
unless stated, medians of 3 runs, idle machine; they rank configurations
against each other inside one scan and are not the experimenter's numbers.

## The cost model

From `bin/kernel_bench` at 32-cell groups, per cell: forcing 225 ns,
Jacobian assembly 438 ns, factorization 676 ns, triangular solve 279 ns. A
MICM three-stage step (1 Jacobian, 1 factorization, 2 forcings, 3 solves)
is 2,400 ns, of which the matrix is 46 %, the solves 35 %, the forcings
19 %. The levers follow: fewer steps per interval (warm start, controller,
tolerance), fewer matrix builds per step (reuse), fewer stages or forcings
per step (the tableau), and a cheaper matrix (precision).

## Plan

| Scan | Choice varied | Fixed | Measures |
|---|---|---|---|
| S1 warm start | MICM's restart at 1e-6 dt vs the carried step; `h_start` fractions | ros3, rtol 1e-3, atol 1e-18 | steps per interval, time, digits |
| S2 tableau and tolerance | ros2, ros3, rodas3, ros34pw2, ros4, rodas4, extrap (orders 2 to 5) at rtol 1e-2 to 1e-4 | atol 1e-18, warm start, MICM's controller and norm | the work-precision curves: time and counts against digits, overall and per species |
| S3 absolute tolerance | atol 1e-15, 1e-16, 1e-18, 1e-20 | best of S2 at each level | digits per species (the species near the threshold), time |
| S4 controller | micm; predictive; h211b with (b, k) in {(4, p), (1, 1.7), (2, p)}; safety 0.9, 1.1, 1.3, 1.5; fac_max 6, 10 | best tableau per level | steps, rejections, digits |
| S5 error norm | rms; cellmax; max | " | digits per species against time |
| S6 factorization reuse | reuse_age 0, 2, 5, 20, 100; reuse_change 0.05, 0.1, 0.2, 0.5, inf; freeze_band 1, 1.3, 1.5, 2; quantize | ros2 and ros34pw2 (W-methods), and rodas3 as the control that loses order | reuses, Jacobian updates, digits, time |
| S7 clamping | clamp 0/1 | " | gate, digits |
| S8 extrapolation | sequence harmonic/shifted/even; order 2 to 6 | rtol series | the work-precision curves |
| S9 grouping and scheduling | group size 16, 32, 64; static/dynamic | best points | time at 10 and 1 threads |
| S10 generic kernels | generated vs generic | best points | time, identity of results |
| S11 precision (if time allows) | double; single-precision factorization and solves with double residuals | best points | digits, time |

Scenarios: all six disclosed ones for the work-precision scans (S2, S3, S6,
S8); the others on diurnal-10min and polluted-10min, with the winner
confirmed on all six. The held-out behavior is estimated by the spread
across the six.

## Results

*(One section per scan, each with its CSV, the table or plot that answers
it, and the choice it fixed for the next scans.)*
