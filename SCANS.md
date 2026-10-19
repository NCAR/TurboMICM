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
| S4b H211b grid | safety 0.7 to 1.0 against tolerance; b = 3 and 6 | RODAS3; extrapolation at safety 0.9 | the frontier against the elementary controller |
| S12 joint grid | rtol 1e-2 to 2e-4 x atol 1e-18, 1e-20, 1e-22 x controller, with the quantized step | RODAS3; extrapolation order 4 | the frontier that defines the operating points |

S11 would need the generator to emit single-precision variants of the
Jacobian, LU, and solve kernels (the forcing stays double): the matrix
halves in size (37 to 18 KB per block of two cells) and the vector units do
four lanes per instruction instead of two, so the matrix and solves, 76 %
of a RODAS3 step by the cost model, could run up to 1.5 times faster, about
1.3 times for the step. The risk is accuracy: the stage solutions carry the
single-precision solve's error (about 1e-7 relative) straight into the
state, which is harmless for two or three digits if it stays random, and
the EMAC single-precision Rosenbrock runs kept species masses within
0.004 % [SKTC26]; but species near the threshold and a W-method-like loss
of order from the inexact matrix would have to be measured on every
scenario. It is the last lever, taken only if the algorithmic ones leave
it worth a day of generator work.

Scenarios: all six disclosed ones for the work-precision scans (S2, S3, S6,
S8); the others on diurnal-10min and polluted-10min, with the winner
confirmed on all six. The held-out behavior is estimated by the spread
across the six.

## Results

### S1: the first step of each solve (`scans/S1-warm-start.csv`)

ROS3 at rtol 1e-3 and atol 1e-18 on the six disclosed scenarios, 1,000
cells, 10 threads, 3 rounds. Steps and rejections are per scenario step and
group (32 cells); the time is the sum over the scenarios of the median
elapsed time.

| First step | Lowest digits (worst species) | Steps per interval | Rejections per interval | Time, s |
|---|---|---:|---:|---:|
| MICM's: 1e-6 of the interval | 3.34 (2.12) | 13.85 | 0.75 | 6.70 |
| 1e-3 of the interval | 3.33 (2.12) | 11.66 | 1.57 | 5.42 |
| 0.1 of the interval | 3.31 (2.07) | 10.71 | 2.23 | 4.89 |
| the whole interval | 3.27 (2.06) | 9.98 | 2.71 | 4.47 |
| carried from the previous solve | 3.33 (2.12) | 9.76 | 2.44 | 4.34 |
| carried, steps tiling the interval | 3.46 (2.25) | 9.73 | 2.36 | 4.32 |

- **The carried step is adopted** (with quantization, which costs nothing
  and is what the reuse rule needs): 35 % less time than MICM's restart at
  the same accuracy. On Tier 1's 30-second steps the gain is larger (1.43
  against 9.29 internal steps per step); on the day-night scenarios the
  solver needs about ten internal steps per 10-minute interval at this
  tolerance whatever the start, so the start-up is a smaller share there.
- **Rejections are a quarter of the attempts** once the step is carried
  (2.4 of 9.7 per interval), against 5 % with MICM's tiny restart: MICM's
  controller lets the step grow six-fold after every accepted step and the
  next attempt often fails. That is the controller scan's target (S4).
- The internal step at rtol 1e-3 is about a minute of simulated time on
  every day-night scenario (17 per 30-minute interval, 10.7 per 10-minute,
  7 per 5-minute), so at this accuracy the cost is set by the chemistry's
  own time scales, and fewer, cheaper, or larger steps are the levers:
  order (S2), the controller (S4), and factorization reuse (S6).

### S2: tableau and relative tolerance (`scans/S2-tableau-tolerance.csv`)

Seventy configurations on the six disclosed scenarios, 1,000 cells, 10
threads, one round: six tableaus (ROS2, ROS3, RODAS3, ROS34PW2, ROS4,
RODAS4) and extrapolated linearly implicit Euler at orders 2 to 5, each at
relative tolerances from 1e-2 to 1e-4, atol 1e-18, the carried step without
quantization (the driver's default; S1's "warm-quantized" row and S3 to S12
quantize, which adds about 0.1 digit for the worst species at equal time),
MICM's controller and norm. "Lowest" is the minimum over the six scenarios, the
ranking's view; the time is the sum over the scenarios of the elapsed time
(`tools/summarize_scan.py scans/S2-tableau-tolerance.csv --sort time`).

The fastest configurations that reach each level on all six:

| Level | Configuration | Lowest digits (worst species) | Time, s | Steps per interval | Factorizations per interval |
|---|---|---|---:|---:|---:|
| 2 | RODAS3, rtol 5e-3 | 2.04 (1.13) | 1.20 | 2.1 | 2.1 |
| 2 | ROS3, rtol 1e-2 | 2.42 (1.46) | 1.49 | 3.4 | 3.4 |
| 2 | RODAS3, rtol 2e-3 | 2.60 (1.47) | 1.54 | 2.7 | 2.7 |
| 2 | extrapolation order 3, rtol 1e-2 | 2.34 (1.32) | 1.60 | 1.4 | 4.2 |
| 3 | RODAS3, rtol 5e-4 | 3.04 (2.07) | 2.65 | 4.7 | 4.7 |
| 3 | extrapolation order 4, rtol 5e-4 | 3.19 (2.34) | 2.73 | 1.7 | 6.8 |
| 3 | ROS4, rtol 5e-3 | 3.03 (2.04) | 2.95 | 5.4 | 5.4 |
| 3 | extrapolation order 3, rtol 5e-4 | 3.29 (2.18) | 3.10 | 3.1 | 9.2 |
| 3 | extrapolation order 4, rtol 2e-4 | 3.87 (2.73) | 3.32 | 2.2 | 8.7 |
| 4 | extrapolation order 5, rtol 1e-4 | 4.02 (3.15) | 3.75 | 1.8 | 8.8 |

- **Third-order Rosenbrock with one factorization per step is the cheapest
  way to two and three digits**, as the 1997 benchmarks found [S97b]:
  RODAS3 leads both levels, ROS3 is close at level 2, and ROS4 at level 3.
  ROS2 is never competitive (order 2 needs 14 to 40 steps per interval at
  these tolerances), nor RODAS4 (six forcing evaluations per step).
- **ROS34PW2 needs about as many steps as ROS3 for the same digits but a
  fourth forcing evaluation per step**, so without matrix reuse it is 30 %
  slower than ROS3; its value, if any, is as the W-method that may keep a
  factorization (S6).
- **Extrapolation takes the largest steps** (1.1 to 1.7 per interval at
  orders 4 and 5 until the tolerance is 2e-4 or tighter) but pays four or
  five factorizations per step, and its accuracy hardly moves with the
  tolerance between 1e-2 and 1e-3: on polluted-10min it gives 2.85 digits at
  any of those, with one step per interval. Its error estimate, the
  difference of the last two tableau entries, does not see the error that a
  single 600-second step with the Jacobian of the step's start makes on
  stiff chemistry that moves within the step, the perturbed asymptotics of
  stiff problems [HW96, IV.9]; the estimate becomes reliable once the
  tolerance forces two or more steps per interval. It overtakes RODAS3 from
  about 3.5 digits on and reaches four digits (3.15 per species) in 3.75 s,
  so it is the method for the accurate end.
- **The 30-minute scenario sets the lowest digits for every method** (RODAS3
  needs 4.3 steps per interval there at 5e-3 and 9.8 at 5e-4; the 10-minute
  scenarios take 1.2 to 2.7 and 2.7 to 5.5), and isoprene is the worst
  species almost everywhere, followed by the halogen radicals BrO and OClO
  and the organic peroxy radicals on the short-step scenario.
- **Margins:** RODAS3 at 5e-3 clears level 2 by 0.04 digit overall and 0.13
  per species, and at 5e-4 clears level 3 by 0.04 and 0.07; the undisclosed
  scenarios (other seasons, latitudes, 48-hour runs) will move the lowest
  digits by the spread seen among the six, about 0.1 to 0.2. The declared
  points will need the next tolerance down or a margin bought otherwise,
  which S4 and S6 try to make cheap: fewer rejections at a given tolerance,
  and factorizations reused.

### S6a: frozen matrices, a first look (`scans/S6a-frozen-first-look.csv`)

Three scenarios (diurnal-10min, polluted-10min, diurnal-30min), 1,000
cells, 10 threads, one round. "Reuse" keeps every factorization while the
step is unchanged, the age is below the limit, and the rate constants have
changed by less than the threshold since the factorization.

| Configuration | diurnal-10min: digits (worst), ms | polluted-10min | diurnal-30min | Factorizations per interval (diurnal-10min) |
|---|---|---|---|---:|
| extrapolation order 3, rtol 1e-2 | 2.67 (1.49), 243 | 2.60 (1.74), 275 | 2.33 (1.40), 155 | 3.5 |
| the same, reuse up to 100 steps, 5 % rate change | 2.60 (1.50), 200 | 2.24 (1.45), 219 | 2.30 (1.43), 142 | 1.6 |
| the same, 20 % | 2.53 (1.50), 197 | 2.22 (1.45), 228 | 2.33 (1.47), 147 | 1.3 |
| extrapolation order 4, rtol 1e-2 | 3.21 (2.20), 349 | 2.85 (1.62), 386 | 2.31 (1.05), 173 | 4.5 |
| the same, reuse, 5 % | 3.11 (2.17), 287 | 2.36 (1.45), 319 | 2.31 (1.05), 169 | 2.0 |
| the same, 20 % | 3.02 (1.89), 269 | 2.13 (1.18), 313 | 2.30 (1.05), 168 | 1.5 |
| RODAS3, rtol 1e-2 | 2.06 (1.26), 160 | 2.36 (1.28), 193 | 1.96 (1.04), 124 | 1.4 |
| the same, reuse, 20 % | 1.74 (0.77), 150 | 1.74 (0.85), 192 | 1.91 (0.83), 130 | 0.6 |
| ROS34PW2, rtol 1e-2, reuse, 20 % | 1.97 (1.17), 327 | 1.91 (1.21), 402 | 2.19 (1.16), 271 | 1.9 |

- **Reuse buys the extrapolation about 18 % and costs accuracy**: the
  stale Jacobian loses up to half a digit on polluted-10min, and the error
  estimate asks for more steps (1.16 to 1.24 per interval). The
  factorization count per interval cannot fall below about 1.3, because by
  day the photolysis rates change by more than 5 % from one interval to the
  next and the matrix must be rebuilt, and at 20 % the accuracy loss grows
  while the time does not improve further.
- **RODAS3 with a reused matrix loses its order**, as a method that is not
  a W-method must: below two digits at rtol 1e-2 with no time gained, since
  the error estimate then asks for more steps.
- **ROS34PW2 with reuse stays twice as slow as RODAS3 without**: it needs
  three to eight steps per interval at this tolerance.
- **Decision:** factorization reuse is not a lever for this problem at
  these levels, because the rate constants change every interval by day and
  the W-methods available need more steps than RODAS3 saves in
  factorizations. The option stays, documented and tested, for mechanisms
  or models whose rates change more slowly. The levers that remain are
  fewer steps and rejections at a given accuracy (S4), the absolute
  tolerance (S3), and the cost of a step.

### S4: the step-size controller (`scans/S4-controller.csv`)

106 configurations on diurnal-10min and polluted-10min, 1,000 cells, 10
threads, one round, the carried step without quantization as in S2: RODAS3
at rtol 3e-3 and 5e-4, ROS2 at 3e-3, ROS34PW2 at 1e-2 and 1e-3, and
extrapolation of order 4 at 2e-3, each with MICM's elementary controller
(safety 0.9, 1.2, 1.4, and the growth limit raised from 6 to 10), the
predictive controller [G94] at the same safeties, H211b [So03] with
(b, k) = (4, p), (1, 1.7) [D25], and (2, p), where p is the estimator's
order, and the filtered error estimate ((alpha I - J)^-1 alpha e, one extra
solve) at safety 0.9 with the elementary and the (1, 1.7) controllers. Time
is the sum over the two scenarios; steps and rejections are per interval
and 32-cell group, the mean of the two scenarios.

RODAS3 at rtol 3e-3:

| Controller | Safety | Lowest digits (worst species) | Time, s | Steps per interval | Rejections per interval (share) |
|---|---:|---|---:|---:|---|
| elementary | 0.9 | 2.49 (1.40) | 0.510 | 2.39 | 0.33 (14 %) |
| elementary | 1.2 | 2.42 (1.30) | 0.497 | 2.28 | 0.44 (19 %) |
| elementary | 1.4 | 2.42 (1.30) | 0.505 | 2.33 | 0.54 (23 %) |
| elementary, fac_max 10 | 0.9 | 2.49 (1.41) | 0.527 | 2.41 | 0.34 (14 %) |
| predictive | 0.9 | 2.50 (1.42) | 0.513 | 2.39 | 0.20 (8 %) |
| predictive | 1.2 | 2.42 (1.30) | 0.491 | 2.28 | 0.32 (14 %) |
| predictive | 1.4 | 2.42 (1.30) | 0.519 | 2.36 | 0.44 (18 %) |
| H211b (4, p) | 0.9 | 2.90 (1.93) | 0.674 | 3.12 | 0.31 (10 %) |
| H211b (4, p) | 1.2 | 2.61 (1.71) | 0.580 | 2.62 | 0.39 (15 %) |
| H211b (4, p) | 1.4 | 2.50 (1.43) | 0.547 | 2.48 | 0.46 (19 %) |
| H211b (1, 1.7) | 0.9 | 2.42 (1.30) | 0.501 | 2.29 | 0.51 (22 %) |
| H211b (1, 1.7) | 1.2 | 2.41 (1.30) | 0.493 | 2.26 | 0.53 (23 %) |
| H211b (1, 1.7) | 1.4 | 2.41 (1.30) | 0.478 | 2.26 | 0.54 (24 %) |
| H211b (2, p) | 0.9 | 2.53 (1.47) | 0.560 | 2.56 | 0.37 (14 %) |
| H211b (2, p) | 1.2 | 2.42 (1.30) | 0.509 | 2.38 | 0.45 (19 %) |
| H211b (2, p) | 1.4 | 2.42 (1.30) | 0.505 | 2.35 | 0.48 (21 %) |
| elementary, filtered estimate | 0.9 | 2.01 (1.21) | 0.376 | 1.58 | 0.10 (6 %) |
| H211b (1, 1.7), filtered estimate | 0.9 | 2.02 (1.21) | 0.351 | 1.48 | 0.14 (10 %) |

RODAS3 at rtol 5e-4:

| Controller | Safety | Lowest digits (worst species) | Time, s | Steps per interval | Rejections per interval (share) |
|---|---:|---|---:|---:|---|
| elementary | 0.9 | 3.04 (2.07) | 1.054 | 5.00 | 0.86 (17 %) |
| elementary | 1.2 | 2.93 (2.14) | 1.204 | 5.67 | 1.68 (30 %) |
| elementary | 1.4 | 2.92 (2.07) | 1.308 | 6.16 | 2.19 (36 %) |
| predictive | 0.9 | 3.07 (2.13) | 1.087 | 5.19 | 0.68 (13 %) |
| predictive | 1.2 | 2.92 (2.20) | 1.234 | 5.82 | 1.49 (26 %) |
| H211b (4, p) | 0.9 | 3.37 (2.52) | 1.403 | 6.68 | 0.86 (13 %) |
| H211b (4, p) | 1.2 | 3.17 (2.20) | 1.229 | 5.76 | 1.17 (20 %) |
| H211b (4, p) | 1.4 | 2.95 (2.18) | 1.270 | 6.06 | 1.63 (27 %) |
| H211b (1, 1.7) | 0.9 | 2.92 (2.14) | 1.072 | 5.20 | 1.56 (30 %) |
| H211b (2, p) | 0.9 | 3.08 (2.09) | 1.190 | 5.59 | 0.97 (17 %) |
| elementary, filtered estimate | 0.9 | 2.72 (1.91) | 0.756 | 3.18 | 0.40 (13 %) |

- **Controllers are compared on the work-precision frontier, not at a fixed
  tolerance**, since each maps the tolerance to a different step sequence.
  Against S2's RODAS3 points on the same two scenarios (the elementary
  controller at rtol 1e-2 to 1e-4), one controller lies above the curve:
  **H211b with b = 4 at safety 0.9 or 1.2**. At 3e-3 and safety 0.9 it
  reaches 2.90 (1.93) in 0.674 s, where the elementary controller needs
  rtol 1e-3 and 0.811 s for 2.76 (1.82); at safety 1.2 it reaches
  2.61 (1.71) in 0.580 s against 2.60 (1.55) in 0.622 s at rtol 2e-3; at
  5e-4 and safety 0.9 it reaches 3.37 (2.52) in 1.40 s against 3.68 (2.54)
  in 1.63 s at rtol 2e-4. That is 7 to 17 % less time at equal worst-species
  digits. Its step sequence is smoother and more conservative: 30 % more
  steps than the elementary controller at the same tolerance, half the
  rejection share (10 to 13 %), and 0.4 to 0.5 more digits for the worst
  species.
- **The predictive controller halves the rejections (8 % and 13 %) but
  gains nothing**: it takes the same number of steps for the same digits,
  smaller ones after each acceptance instead of a rejected one.
- **H211b with (1, 1.7)** [D25] behaves like the elementary controller here
  (same digits and time, more rejections), because b = 1 removes the
  smoothing term and k = 1.7 below p makes the response more aggressive.
  (2, p) sits between.
- **Safety above 1** saves 2 to 3 % at 3e-3 and costs 10 to 20 % at 5e-4
  (rejections 26 to 36 %); 0.9 stays. The growth limit (6 against 10)
  makes no difference.
- **The filtered estimate** shrinks the error estimate, and the controller
  takes about a third fewer steps at the same tolerance; the result lands on
  the elementary controller's curve (3e-3 filtered: 2.01 (1.21) in 0.376 s
  against rtol 1e-2 unfiltered: 2.05 (1.26) in 0.381 s) or just above it
  (5e-4 filtered: 2.72 (1.91) in 0.756 s against 2.76 (1.82) in 0.811 s at
  rtol 1e-3). One extra solve per step buys nothing a looser tolerance
  does not; not adopted.
- **ROS34PW2 and ROS2** confirm S2 with every controller: ROS34PW2 costs 1.4
  to 1.9 times RODAS3's time at equal digits; ROS2 needs 12 to 21 steps per
  interval, and with safety 1.2 or 1.4 the driver exited with an error on
  polluted-10min (status `exit 1` in the CSV).
- **Extrapolation of order 4 at 2e-3** gives 2.91 (1.62) in 0.77 to 0.90 s
  with every controller but H211b (4, p), which takes 2.0 steps per interval
  instead of 1.3 and reaches 3.31 (2.20) in 1.20 s, a point that order 3 at
  5e-4 with the elementary controller dominates (3.53 (2.36) in 1.21 s): no
  gain for extrapolation.
- **Decision:** H211b with b = 4, k = p, safety 0.9 is the candidate
  controller for RODAS3, checked across tolerances and safeties in S4b and
  under the final conventions (quantized step, joint with the absolute
  tolerance) in S12; the elementary controller stays for extrapolation.

### S5: the error norm (`scans/S5-norm.csv`)

21 configurations on diurnal-10min and polluted-10min, quantized step,
elementary controller: RODAS3 at rtol 1e-2 to 5e-4 with the group RMS norm
(MICM's; "rms"), the per-cell RMS maximized over the group's cells
("cellmax"), and the maximum over all components ("max"); extrapolation of
order 4 at 2e-3 to 5e-4 with rms and max.

| Norm | rtol | Lowest digits (worst species) | Time, s | Steps per interval |
|---|---|---|---:|---:|
| rms | 1e-2 | 2.06 (1.26) | 0.35 | 1.57 |
| rms | 5e-3 | 2.39 (1.28) | 0.42 | 1.88 |
| rms | 2e-3 | 2.67 (1.56) | 0.58 | 2.68 |
| rms | 1e-3 | 2.79 (2.00) | 0.78 | 3.56 |
| rms | 5e-4 | 3.04 (2.21) | 1.02 | 4.95 |
| cellmax | 1e-2 | 2.40 (1.29) | 0.40 | 1.82 |
| cellmax | 5e-3 | 2.62 (1.50) | 0.52 | 2.36 |
| cellmax | 2e-3 | 2.77 (1.98) | 0.72 | 3.37 |
| cellmax | 1e-3 | 2.90 (2.20) | 0.98 | 4.67 |
| cellmax | 5e-4 | 3.42 (2.46) | 1.35 | 6.39 |
| max | 1e-2 | 2.90 (2.20) | 1.00 | 4.66 |
| max | 5e-3 | 3.24 (2.53) | 1.34 | 6.14 |
| max | 2e-3 | 3.96 (2.94) | 1.96 | 9.03 |
| max | 1e-3 | 4.44 (3.21) | 2.62 | 12.31 |
| max | 5e-4 | 4.77 (3.60) | 3.44 | 16.03 |
| extrapolation 4, rms | 2e-3 | 2.91 (1.62) | 0.79 | 1.31 |
| extrapolation 4, rms | 5e-4 | 3.21 (2.35) | 0.98 | 1.56 |
| extrapolation 4, max | 2e-3 | 4.67 (3.52) | 1.62 | 2.67 |
| extrapolation 4, max | 5e-4 | 5.00 (3.84) | 2.96 | 4.94 |

- **The three norms lie on one curve.** The per-cell norm at rtol r matches
  the group norm at about r/2.5 (cellmax 2e-3: 2.77 (1.98) in 0.72 s; rms
  1e-3: 2.79 (2.00) in 0.78 s), 4 to 8 % less time for the same worst
  species and 0.1 digit less overall; the maximum norm at r matches the
  group norm at r/20 (max 1e-2: 2.90 (2.20) in 1.00 s; rms 5e-4: 3.04 (2.21)
  in 1.02 s). Controlling the worst component directly buys nothing,
  because the species that sets the score, isoprene near its night-time
  floor, is not the component with the largest scaled error estimate.
- **For extrapolation the maximum norm** at 2e-3 gives 4.67 (3.52) in
  1.62 s against 4.53 (3.32) in 1.63 s for the group norm at 1e-4 (S2): a
  small gain at the accurate end, not pursued.
- **Decision:** MICM's group RMS norm stays. The per-cell norm's 4 to 8 %
  is within the drift between scans (the same configuration, RODAS3 at
  rtol 1e-2 with the group norm, took 0.381 s in S2 and 0.346 s here, 10 %
  apart, in different sessions).

### S3: the absolute tolerance (`scans/S3-atol.csv`)

18 configurations on diurnal-10min and polluted-10min, quantized step,
elementary controller, RODAS3.

| rtol | atol | Lowest digits (worst species) | Time, s | Steps per interval |
|---|---|---|---:|---:|
| 1e-2 | 1e-15 | 1.98 (1.17) | 0.33 | 1.44 |
| 1e-2 | 1e-16 | 1.99 (1.17) | 0.33 | 1.50 |
| 1e-2 | 1e-18 | 2.06 (1.26) | 0.35 | 1.57 |
| 5e-3 | 1e-14 | 1.99 (1.17) | 0.34 | 1.47 |
| 5e-3 | 1e-15 | 1.99 (1.18) | 0.36 | 1.59 |
| 5e-3 | 1e-16 | 2.21 (1.28) | 0.39 | 1.70 |
| 5e-3 | 1e-18 | 2.39 (1.28) | 0.42 | 1.88 |
| 5e-3 | 1e-20 | 2.42 (1.29) | 0.44 | 2.03 |
| 5e-3 | 1e-22 | 2.46 (1.33) | 0.47 | 2.10 |
| 2e-3 | 1e-15 | 2.39 (1.28) | 0.43 | 1.92 |
| 2e-3 | 1e-16 | 2.42 (1.31) | 0.46 | 2.14 |
| 2e-3 | 1e-18 | 2.67 (1.56) | 0.57 | 2.68 |
| 5e-4 | 1e-14 | 2.50 (1.41) | 0.54 | 2.47 |
| 5e-4 | 1e-15 | 2.69 (1.63) | 0.66 | 3.00 |
| 5e-4 | 1e-16 | 2.77 (1.81) | 0.78 | 3.54 |
| 5e-4 | 1e-18 | 3.04 (2.21) | 1.08 | 4.95 |
| 5e-4 | 1e-20 | 3.39 (2.42) | 1.25 | 5.81 |
| 5e-4 | 1e-22 | 3.58 (2.60) | 1.38 | 6.42 |

- **The absolute tolerance is a lever at level 3.** The score counts
  species down to 1e-14 mol m-3, and the species that set it (isoprene at
  night, phenol) live near that floor, where an absolute tolerance above
  rtol x 1e-14 leaves them weakly controlled. At rtol 5e-3 the gain
  saturates at 1e-18 (1e-20 adds 0.03 digit for 5 % more time). At 5e-4 the
  tightening keeps paying: 1e-20 and 1e-22 add 0.2 digit each for the worst
  species, for 15 and 10 % more time, which beats buying the same digits
  with the relative tolerance (2.6 digits for the worst species by rtol
  alone needs about 2e-4 at 1e-18, 1.63 s in S2, against 1.38 s here). At
  1e-15 and 1e-16, MICM's and the Tier 1 entries' choices, level 3 is out
  of reach at this rtol (2.69 and 2.77 overall).
- **Decision:** 1e-18 for the level-2 points; the level-3 points come from
  the joint grid S12 (rtol 1e-2 to 2e-4 against atol 1e-18, 1e-20, 1e-22,
  both controllers).

### S7: clamping (`scans/S7-clamp.csv`)

RODAS3 at 5e-3 and 5e-4, quantized step, clamping negative concentrations
to zero after each accepted step on and off: identical digits (2.39 (1.28);
3.04 (2.21)) and times within 2 %. No scenario produced a value below the
gate's -1e-20 mol m-3 without clamping (every S1 to S12 run passes the
gate). **Decision:** off, as in MICM; the option stays for hosts that need
non-negative output.

### S9: group size and scheduling (`scans/S9-grouping.csv`, 10 threads)

RODAS3 at 5e-3 (level 2) and 5e-4 (level 3), quantized step, 1,000 cells,
time summed over diurnal-10min and polluted-10min:

| Groups | 8 cells | 16 | 32 | 64 | 128 | 32, static | 32, guided | 32, dynamic chunks of 4 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Level 2, s | 0.40 | 0.40 | 0.41 | 0.48 | 0.53 | 0.45 | 0.42 | 0.42 |
| Level 3, s | 1.05 | 1.03 | 1.05 | 1.19 | 1.34 | 1.10 | 1.06 | 1.07 |

- Digits do not depend on the grouping (the worst species moves by 0.02
  between 8 and 128 cells, as the group's error norm pools more cells).
- Groups of 8 to 32 are equal within noise at 10 threads; 64 and 128 cost
  15 to 30 % more, from load imbalance among 1,000 cells (8 to 16 groups
  for 10 threads) and from every cell in a group taking the steps its
  stiffest member needs. Static scheduling costs 5 to 10 %; guided and
  chunked dynamic scheduling match the default (dynamic, chunk 1).
- **Decision:** 32-cell groups and dynamic scheduling stay, as in Tier 1,
  where the kernels were tuned at 32. The 1-thread run of the same matrix
  (`scans/S9-grouping-1thread.csv`) follows below.

