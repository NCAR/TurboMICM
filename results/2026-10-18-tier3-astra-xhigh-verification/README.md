# Verification of Tier 3 astra-xhigh, 2026-10-18

The third Tier 3 entry, astra-xhigh (`db44838` on `tier3/astra/xhigh`, 6 commits
after its Tier 1 commit `fbd365c`), was verified the same way as sol-xhigh (see
[`2026-10-18-tier3-sol-xhigh-verification`](../2026-10-18-tier3-sol-xhigh-verification/README.md)):
`admin/verify_tier3.py`, from 08:32 to 08:38, before astra-max was launched.

| File | Contents |
|---|---|
| `points.txt` | the declared operating points, copied from the entry's OPERATING_POINTS.md |
| `report.md` | the verification report: build, MICM's test cases, the Tier 1 TS1 gate, the points on the disclosed scenarios, and the build audit |
| `points.json` | one row per point and disclosed scenario: exit status, elapsed_ms, digits, worst species, gate, levels 2 and 3 |

## Builds, tests, and the Tier 1 gate

- **Builds:** `make env`, then `make cpu bench OPTFLAGS="-O3 -DNDEBUG"`, with
  `CC=clang` and with `CC=gcc-15`: no warnings.
- **MICM's test cases:** 680 PASS, 0 FAIL, 0 N/A with each compiler, every one
  of the 680 ids exactly once, as in the entry's Tier 1 solution.
- **Tier 1 TS1 gate:** largest relative difference from MICM's final state
  2.03e-13 with clang, as in Tier 1, and 2.14e-13 with gcc-15 (limit 1e-10).
- **Audit:** only `-std=c11`, warnings, `-O3 -DNDEBUG`, include paths, and
  OpenMP flags, besides the options of its kernel generator and of
  `conda env update`; every program links only the OpenMP runtime and
  `libSystem`; no banned construct in the 77 tracked and generated C files.
- **Timer:** `src/scenario.c` times, with `omp_get_wtime()`, each block of up to
  four steps, cut at snapshots: the photolysis update at the start of every
  step, the rate constants, the integration, and a per-step check for
  non-finite or negative values. Setup, reading the input, and writing
  snapshots and output are outside it, as the prompt requires. It computes
  every rate constant in the first step and, after that, only those that depend
  on photolysis, since the scenarios hold every other condition constant.

## Operating points

Each point ran once on each scenario at 1,000 cells and 10 threads, built with
clang, and was scored with `tier3/scenarios/score.py`. The two points for level
2 use RODAS3, the four-stage, stiffly accurate third-order formula whose
coefficients the entry kept from MICM, in a new integration loop. The three for
level 3 use a method new to the entry: linearly implicit Euler extrapolation to
order 4, which takes 1, 2, 3, and 4 substeps across each step with one Jacobian
and combines the four results by Neville extrapolation (four factorizations and
ten linear solves per attempted step). All five use an absolute tolerance of
1e-18 mol m-3, a new Jacobian at every step, 32-cell groups, four-step blocks,
and dynamic scheduling.

| Point | Method, rtol | Disclosed: lowest overall (worst species) | Level 2 | Level 3 | Held out: lowest overall (worst species) | Level 2 | Level 3 |
|---|---|---|---|---|---|---|---|
| fast | RODAS3, 1e-2 | 2.04 (1.035) | all 6 | 1 of 6 | 2.19 (1.13) | all 5 | none |
| balanced | RODAS3, 3e-3 | 2.65 (1.65) | all 6 | 2 of 6 | 2.46 (1.57) | all 5 | none |
| accurate | extrapolation, order 4, 4e-3 | 3.20 (2.39) | all 6 | all 6 | 3.56 (2.63) | all 5 | all 5 |
| recommended | extrapolation, order 4, 2e-3 | 3.45 (2.55) | all 6 | all 6 | 3.67 (2.79) | all 5 | all 5 |
| strict | extrapolation, order 4, 1e-3 | 4.06 (2.92) | all 6 | all 6 | 3.86 (2.98) | all 5 | all 5 |

- **The gate:** every point passes on all 11 scenarios.
- **Agreement with the agent's report:** the disclosed scores match its
  OPERATING_POINTS.md and RESULTS.md to the third decimal.
- **fast has almost no margin per species.** Its worst species, isoprene on
  diurnal-30min, reaches 1.035 digits, against a threshold of 1; on the
  held-out scenarios its worst reaches 1.13.
- **The held-out scenarios** stay unnamed until the last Tier 3 run ends. fast,
  accurate, and recommended score higher on them than on the disclosed set,
  overall and per species. balanced's lowest held-out score is 0.18 digit lower
  overall and 0.08 lower per species; strict's is 0.20 lower overall and higher
  per species. Every point still reaches its declared level on every held-out
  scenario.
- **Times are not the ranking.** Single runs at 1,000 cells took 0.15 s to
  1.2 s on the disclosed scenarios. The ranking times every entry's points at
  10,000 cells, median of 5, in one session after the last run.
