# Verification of Tier 3 astra-max, 2026-10-18

The fourth and last Tier 3 entry, astra-max (`3cab602` on `tier3/astra/max`, 6
commits after its Tier 1 commit `0dc7f67`), was verified the same way as sol-xhigh
(see [`2026-10-18-tier3-sol-xhigh-verification`](../2026-10-18-tier3-sol-xhigh-verification/README.md)):
`admin/verify_tier3.py`, from 11:27 to 11:30 on the disclosed scenarios, then on the
held-out ones. With the last run over, the held-out scenarios are published:
[`tier3/heldout`](../../tier3/heldout/README.md).

| File | Contents |
|---|---|
| `points.txt` | the declared operating points, copied from the entry's OPERATING_POINTS.md |
| `report.md` | the verification report: build, MICM's test cases, the Tier 1 TS1 gate, the points on the disclosed scenarios, and the build audit |
| `points.json` | one row per point and disclosed scenario: exit status, elapsed_ms, digits, worst species, gate, levels 2 and 3 |
| `heldout-report.md`, `heldout-points.json` | the same for the five held-out scenarios |

## Builds, tests, and the Tier 1 gate

- **Builds:** `make env`, then `make cpu bench OPTFLAGS="-O3 -DNDEBUG"`, with
  `CC=clang` and with `CC=gcc-15`. gcc-15 gives no warnings; clang gives one, an
  optimizer warning that it could not vectorize a loop of the error norm marked
  `#pragma omp simd` (`src/integrator.c:83`).
- **MICM's test cases:** 680 PASS, 0 FAIL, 0 N/A with each compiler, every one
  of the 680 ids exactly once.
- **Tier 1 TS1 gate:** largest relative difference from MICM's final state
  2.96e-13 with clang and 2.75e-13 with gcc-15 (limit 1e-10).
- **Audit:** only `-std=c11`, warnings, `-O3 -DNDEBUG`, include paths, OpenMP
  flags, and a define of the default group size, besides the options of
  `conda env update`; every program links only the OpenMP runtime and
  `libSystem`; no banned construct in the 25 tracked and generated C files.
- **Timer:** `tools/ts1_scenario.c` times each step with `omp_get_wtime()`: the
  photolysis update, every rate constant (recomputed each step), and the
  integration. Setup, reading the input, a per-step scan for non-finite values
  (its own diagnostic; the integrator checks its states itself), snapshot
  copies, and writing the output are outside it.

## Operating points

Each point ran once on each scenario at 1,000 cells and 10 threads, built with
clang, and was scored with `tier3/scenarios/score.py`. All four use RODAS3, the
four-stage, stiffly accurate third-order formula, at relative tolerance 3e-3,
with 32-cell groups, a new Jacobian at every attempted step, and a step-size
prediction carried from one scenario step to the next. They differ in absolute
tolerance and in a cap on the internal step.

| Point | atol, step cap | Disclosed: lowest overall (worst species) | Level 2 | Level 3 | Held out: lowest overall (worst species) | Level 2 | Level 3 |
|---|---|---|---|---|---|---|---|
| economy | 1e-16, none | 2.27 (1.27) | all 6 | 1 of 6 | 2.28 (1.27) | all 5 | none |
| fast3 | 1e-18, 300 s | 3.09 (2.06) | all 6 | all 6 | 2.46 (1.57) | all 5 | **3 of 5** |
| balanced (recommended) | 1e-18, 200 s | 3.43 (2.51) | all 6 | all 6 | 3.27 (2.26) | all 5 | all 5 |
| conservative | 1e-18, 100 s | 3.87 (3.065) | all 6 | all 6 | 3.775 (2.95) | all 5 | all 5 |

- **The gate:** every point passes on all 11 scenarios.
- **Agreement with the agent's report:** the disclosed scores match its
  RESULTS.md to the third decimal.
- **fast3, declared for level 3, misses it on two held-out scenarios:**
  clean-48h-5min (2.46 overall, phenol 1.57) and winter60-10min (2.94 overall,
  isoprene 1.87). On clean-48h-5min, whose steps are 5 minutes, its 300-second
  cap never binds, and it scores the same as economy. It is the first declared
  point in Tier 3 to miss its level on a held-out scenario. The agent's
  OPERATING_POINTS.md says fast3 "approaches the three-digit acceptance
  boundary on some disclosed and independent cases"; its own three 48-hour
  cases, at 96 cells, put fast3's worst species at 2.057.
- **balanced, the recommended point, reaches level 3 on all 11.** Its lowest
  held-out score, on winter60-10min, is 0.16 digit below its lowest disclosed
  one overall and 0.25 per species. economy and conservative reach their levels
  on all 11.
- **Times are not the ranking.** Single runs at 1,000 cells took 0.23 s to
  1.5 s on the disclosed scenarios. The ranking times every entry's points at
  10,000 cells, median of 5, in one session.
