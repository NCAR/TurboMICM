# Verification of Tier 3 sol-max, 2026-10-18

The second Tier 3 entry, sol-max (`c122d12` on `tier3/sol/max`, 4 commits after
its Tier 1 commit `001fb64`), was verified the same way as sol-xhigh (see
[`2026-10-18-tier3-sol-xhigh-verification`](../2026-10-18-tier3-sol-xhigh-verification/README.md)):
`admin/verify_tier3.py`, from 06:15 to 06:19, before astra-xhigh was launched.

| File | Contents |
|---|---|
| `points.txt` | the declared operating points, copied from the entry's OPERATING_POINTS.md |
| `report.md` | the verification report: build, MICM's test cases, the Tier 1 TS1 gate, the points on the disclosed scenarios, and the build audit |
| `points.json` | one row per point and disclosed scenario: exit status, elapsed_ms, digits, worst species, gate, levels 2 and 3 |

## Builds, tests, and the Tier 1 gate

- **Builds:** `make env`, then `make cpu bench OPTFLAGS="-O3 -DNDEBUG"`, with
  `CC=clang` and with `CC=gcc-15`: no warnings.
- **MICM's test cases:** 675 PASS, 0 FAIL, 5 N/A with each compiler, every one
  of the 680 ids exactly once. The five N/A cases are C++ compile-time
  type-category tests.
- **Tier 1 TS1 gate:** largest relative difference from MICM's final state
  1.26e-13 with clang and 1.16e-13 with gcc-15 (limit 1e-10).
- **Audit:** only `-std=c11`, warnings, `-O3 -DNDEBUG`, include paths, a
  group-size define, OpenMP and shared-library flags, and its kernel
  generator's options; every program links only its own library, the OpenMP
  runtime, and `libSystem`; no banned construct in the 31 tracked and generated
  C files.

## Operating points

Each point ran once on each scenario at 1,000 cells and 10 threads, built with
clang, and was scored with `tier3/scenarios/score.py`. All five points use
RODAS3, the four-stage, stiffly accurate third-order formula whose coefficients
the entry kept from MICM, in a new integration loop with its own error norm
(each cell's RMS, maximum over the group) and a positivity check. They differ
only in relative tolerance, with an absolute tolerance of 1e-18 mol m-3,
16-cell groups, and dynamic scheduling.

| Point | Method, rtol | Disclosed: lowest overall (worst species) | Level 2 | Level 3 | Held out: lowest overall (worst species) | Level 2 | Level 3 |
|---|---|---|---|---|---|---|---|
| fast | RODAS3, 1e-2 | 2.04 (1.06) | all 6 | 1 of 6 | 2.18 (1.09) | all 5 | none |
| balanced | RODAS3, 5e-3 | 2.29 (1.39) | all 6 | 1 of 6 | 2.34 (1.39) | all 5 | none |
| accurate | RODAS3, 5e-4 | 3.27 (2.34) | all 6 | all 6 | 3.25 (2.39) | all 5 | all 5 |
| recommended | RODAS3, 3e-4 | 3.57 (2.55) | all 6 | all 6 | 3.69 (2.61) | all 5 | all 5 |
| precise | RODAS3, 1e-4 | 4.25 (3.07) | all 6 | all 6 | 4.31 (3.08) | all 5 | all 5 |

- **The gate:** every point passes on all 11 scenarios.
- **Agreement with the agent's report:** the disclosed scores match its
  OPERATING_POINTS.md and RESULTS.md to the third decimal.
- **fast has almost no margin per species.** Its worst species reaches 1.06
  digits on the disclosed scenarios and 1.09 on the held-out ones, against a
  threshold of 1.
- **The held-out scenarios** stay unnamed until the last Tier 3 run ends. Each
  point's lowest held-out score, overall and per species, is within 0.1 digit
  of its lowest disclosed score or above it.
- **Times are not the ranking.** Single runs at 1,000 cells took 0.16 s to
  2.4 s on the disclosed scenarios. The ranking times every entry's points at
  10,000 cells, median of 5, in one session after the last run.
