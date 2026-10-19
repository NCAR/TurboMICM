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
relative tolerances from 1e-2 to 1e-4, atol 1e-18, the carried step, MICM's
controller and norm. "Lowest" is the minimum over the six scenarios, the
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
