# Verification of Tier 3 sol-xhigh, 2026-10-18

The first Tier 3 entry, sol-xhigh (`63d246f` on `tier3/sol/xhigh`, 14 commits
after its Tier 1 commit `171f54a`), was rebuilt from its branch with both
compilers, retested, and its five declared operating points were run on every
scenario, by `admin/verify_tier3.py`, from 04:17 to 04:22, before sol-max was
launched.

| File | Contents |
|---|---|
| `points.txt` | the declared operating points, copied from the entry's OPERATING_POINTS.md |
| `report.md` | the verification report: build, MICM's test cases, the Tier 1 TS1 gate, the points on the disclosed scenarios, and the build audit |
| `points.json` | one row per point and disclosed scenario: exit status, elapsed_ms, digits, worst species, gate, levels 2 and 3 |
| `heldout-report.md`, `heldout-points.json` | the same for the five held-out scenarios, added when they were published after the last run ([`tier3/heldout`](../../tier3/heldout/README.md)) |

## Builds, tests, and the Tier 1 gate

- **Builds:** `make env`, then `make cpu bench OPTFLAGS="-O3 -DNDEBUG"`, with
  `CC=clang` and with `CC=gcc-15`: no warnings.
- **MICM's test cases:** 675 PASS, 0 FAIL, 5 N/A with each compiler, every one
  of the 680 ids exactly once. The five N/A cases are C++ compile-time
  type-category tests, as in the entry's Tier 1 solution.
- **Tier 1 TS1 gate:** 10,000 cells of seed 2026, 100 steps of 30 s, against
  MICM's final state: largest relative difference 0 with clang and 1.43e-13 with
  gcc-15 (limit 1e-10).
- **Audit:** only `-std=c11`, warnings, `-O3 -DNDEBUG`, include paths, a
  group-size define, and OpenMP flags; every program links only the OpenMP
  runtime and `libSystem`; no banned construct in the 57 tracked and generated
  C files.
- **Timer:** `src/scenario.c` times, with `omp_get_wtime()`, each stretch of
  steps between snapshots: the photolysis update, the rate constants, the
  integration, and a per-step check for non-finite or negative values. Reading
  the input, setup, copying snapshots, and writing output are outside it, as
  the prompt requires.

## Operating points

Each point ran once on each scenario at 1,000 cells and 10 threads, built with
clang, and was scored with `tier3/scenarios/score.py`. The entry declares four
RODAS3.1 points (coefficients from KPP's documentation) and one
six-stage RODAS4 point, all with 32-cell groups and guided scheduling.

| Point | Method, rtol, atol | Disclosed: lowest overall (worst species) | Level 2 | Level 3 | Held out: lowest overall (worst species) | Level 2 | Level 3 |
|---|---|---|---|---|---|---|---|
| L2 | RODAS3.1, 2e-3, 1e-18 | 2.54 (1.45) | all 6 | 2 of 6 | 2.46 (1.47) | all 5 | none |
| L2_margin | RODAS3.1, 1e-3, 1e-18 | 2.80 (1.69) | all 6 | 5 of 6 | 2.82 (1.70) | all 5 | 1 of 5 |
| L3 | RODAS3.1, 3e-4, 1e-16 | 3.34 (2.32) | all 6 | all 6 | 3.39 (2.31) | all 5 | all 5 |
| recommended | RODAS3.1, 3e-4, 1e-18 | 3.71 (2.73) | all 6 | all 6 | 3.93 (2.75) | all 5 | all 5 |
| L4 | RODAS4, 1e-4, 1e-18 | 4.66 (3.61) | all 6 | all 6 | 4.72 (3.61) | all 5 | all 5 |

- **The gate:** every point passes on all 11 scenarios: no NaN, no infinity,
  and no value below -1e-20 mol m-3.
- **Agreement with the agent's report:** the disclosed scores match its
  RESULTS.md to the third decimal.
- **The worst species is isoprene almost everywhere**, as for MICM's own
  three-stage Rosenbrock at loose tolerances.
- **The held-out scenarios,** made from a script fixed before any Tier 3
  result, are not named here until the last Tier 3 run ends. Each point's
  lowest held-out score, overall and per species, is within 0.1 digit of its
  lowest disclosed score or above it. So the entry's accuracy carried over to
  conditions it had not seen.
- **Times are not the ranking.** These are single runs at 1,000 cells, from
  0.24 s to 4.5 s on the disclosed scenarios. The ranking times every entry's
  points at 10,000 cells, median of 5, in one session after the last run.
