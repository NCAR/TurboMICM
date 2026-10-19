# Verification of Fable's Tier 3 reference solution, 2026-10-19

Fable's own answer to the Tier 3 prompt, built outside the contest (`ec34c99` on
`tier3/fable/max`, 26 commits after the Tier 1 reference entry's head `6d2fcde`;
wiki: [Fable's Tier 3 reference solution](https://github.com/NCAR/TurboMICM/wiki/Fable-Tier-3-Reference-Solution)),
was verified the same way as the four contest entries
(`admin/verify_tier3.py`; see
[`2026-10-18-tier3-sol-xhigh-verification`](../2026-10-18-tier3-sol-xhigh-verification/README.md)):
from 07:34 to 07:43 on the disclosed scenarios with both compilers, then on the
five held-out scenarios ([`tier3/heldout`](../../tier3/heldout/README.md)). It is
not a contestant: not ranked, not blind-reviewed.

| File | Contents |
|---|---|
| `points.txt` | the five declared operating points, from the entry's OPERATING_POINTS.md |
| `report.md` | the verification report: build, MICM's test cases, the Tier 1 TS1 gate, the points on the disclosed scenarios, and the build audit |
| `points.json` | one row per point and disclosed scenario: exit status, elapsed_ms, digits, worst species, gate, levels 2 and 3 |
| `heldout-report.md`, `heldout-points.json` | the same for the five held-out scenarios |

## Builds, tests, and the Tier 1 gate

- **Builds:** `make env`, then `make cpu bench OPTFLAGS="-O3 -DNDEBUG"`, with
  `CC=clang` (127 s) and `CC=gcc-15` (54 s). No warnings with either.
- **MICM's test cases:** 673 PASS, 0 FAIL, and 7 N/A (C++-only cases) on the
  680 ids with each compiler, every id exactly once, plus two PASS lines with
  ids outside the list (the entry's own additional cases), as in its Tier 1
  verification.
- **Tier 1 TS1 gate:** largest relative difference from MICM's final state
  1.04e-13 with clang and 9.97e-14 with gcc-15 (limit 1e-10); the Tier 1
  benchmark itself took 4,597 ms (clang) and 7,089 ms (gcc-15) at 10 threads.
- **Audit:** only `-std=c11`, warnings, `-O3 -DNDEBUG`, include paths, OpenMP
  flags, and two defines (the vector width and the default group size), besides
  the options of the kernel generator and of `conda env create`; every program
  links only the OpenMP runtime and `libSystem`; the 18 lines matching a banned
  construct are test names that mirror MICM's (`...ThrowsException`,
  `vector_size...`) and a comment, in `logs/banned.txt` of the work directory.
- **Timer:** `bin/ts1_scenario` reproduces MICM's scenario driver: `elapsed_ms`
  covers the per-step photolysis update, the rate constants that depend on it,
  and the integration; reading the input, the snapshots, and writing the output
  are outside it.

## Operating points

Each point ran once on each scenario at 1,000 cells and 10 threads, built with
clang, and was scored with `tier3/scenarios/score.py`. All five carry the step
from one interval to the next, quantized to the interval, with 32-cell groups
and the group RMS error norm; they differ in method, controller, and
tolerances.

| Point | Method, relative / absolute tolerance | Declared | Disclosed: lowest overall (worst species) | Level 2 | Level 3 | Held out: lowest overall (worst species) | Level 2 | Level 3 |
|---|---|---|---|---|---|---|---|---|
| L2 | RODAS3, MICM's elementary controller, 5e-3 / 1e-20 | 2 | 2.42 (isoprene 1.29) | all 6 | 1 of 6 | 2.27 (isoprene 1.16) | all 5 | none |
| balanced | RODAS3, H211b controller, 2e-3 / 1e-20 | 2 | 3.17 (isoprene 2.17) | all 6 | all 6 | 2.71 (phenol 1.79) | all 5 | **4 of 5** |
| L3 | extrapolation, order 3, shifted sequence, 1e-3 / 1e-20 | 3 | 3.61 (isoprene 2.49) | all 6 | all 6 | 3.56 (isoprene 2.47) | all 5 | all 5 |
| recommended | extrapolation, order 4, shifted sequence, 5e-4 / 1e-22 | 3 | 3.62 (OClO 2.69) | all 6 | all 6 | 3.92 (NO 3.20) | all 5 | all 5 |
| accurate | extrapolation, order 4, shifted sequence, 2e-4 / 1e-20 | 3 | 3.87 (BrO 3.00) | all 6 | all 6 | 4.27 (OClO 3.42) | all 5 | all 5 |

- **The gate:** every point passes on all 11 scenarios.
- **Every point reaches its declared level on all 11 scenarios.** The three
  level-3 points keep 0.47 to 1.2 digits of margin per species on the held-out
  set; L2 keeps 0.16 (isoprene on winter60-10min), below the 0.2 its ACCURACY.md
  asks on the disclosed scenarios but above the threshold.
- **balanced**, declared for level 2 because its three-digit margin on the
  disclosed scenarios was 0.17, misses level 3 on clean-48h-5min (2.71 overall,
  phenol 1.79) and reaches it on the other four; a ranking would credit it at
  level 2, as declared. Phenol on the clean 48-hour scenario is the species
  that cost astra-max's fast3 its level too.
- **Agreement with the entry's own report:** the disclosed scores match its
  RESULTS.md (scan S13) to the third decimal; the single-run times here are
  within 2 % of that scan's medians.
- **Times** (single runs, sums over the scenarios, ms): disclosed, L2 1,175,
  balanced 2,221, L3 2,681, recommended 3,097, accurate 3,497; held out,
  L2 1,439, balanced 2,646, L3 3,206, recommended 3,582, accurate 4,212. They
  are not comparable with the contest's ranking session (10,000 cells, five
  rounds, its own night); only a session that times this solution beside the
  contestants' points would place it among them.
