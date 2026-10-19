# Tier 3 ranking, 2026-10-19

The primary ranking of the four Tier 3 entries: each entry's fastest declared
operating point that reaches level 2 (two digits overall, one per species) and
level 3 (three and two) on every one of the 11 scenarios, six disclosed and
five held out, at 10,000 cells. The rules were fixed before the runs
(`tier3/HANDOFF.md` until this commit, "Decisions already made"; now in
`tier3/README.md`).

## Results

Time is the sum over the 11 scenarios of the point's median elapsed time over
5 rounds at 10 threads, in seconds; the sums over the six disclosed and the
five held-out scenarios are shown beside it. Speedup is against MICM's
three-stage Rosenbrock at the loosest tolerance that reaches the level
(`micm/L2`: relative 1e-2, absolute 1e-15; `micm/L3`: 3e-4 and 1e-15), with
32-cell groups and the same photolysis update per step.

**Level 2** (2 digits overall, 1 per species, on all 11 scenarios):

| Rank | Entry | Point | Method | Time (s) | Disclosed (s) | Held out (s) | Speedup over MICM |
|---:|---|---|---|---:|---:|---:|---:|
| 1 | astra-xhigh | fast | RODAS3, relative 1e-2, absolute 1e-18 | 34.72 | 16.15 | 18.58 | 13.4 |
| 2 | sol-max | fast | RODAS3, relative 1e-2, 16-cell groups | 39.80 | 18.79 | 21.01 | 11.7 |
| 3 | sol-xhigh | L2 | RODAS3.1, relative 2e-3, absolute 1e-18 | 54.52 | 25.79 | 28.74 | 8.5 |
| 4 | astra-max | economy | RODAS3, relative 3e-3, absolute 1e-16, no step cap | 62.13 | 29.02 | 33.12 | 7.5 |
| | MICM | L2 | three-stage Rosenbrock, relative 1e-2, absolute 1e-15 | 465.40 | 221.06 | 244.34 | 1 |

**Level 3** (3 digits overall, 2 per species, on all 11 scenarios):

| Rank | Entry | Point | Method | Time (s) | Disclosed (s) | Held out (s) | Speedup over MICM |
|---:|---|---|---|---:|---:|---:|---:|
| 1 | astra-xhigh | accurate | extrapolated linearly implicit Euler, order 4, relative 4e-3, absolute 1e-18 | 74.47 | 34.42 | 40.05 | 7.3 |
| 2 | sol-xhigh | L3 | RODAS3.1, relative 3e-4, absolute 1e-16 | 93.67 | 43.58 | 50.09 | 5.8 |
| 3 | sol-max | accurate | RODAS3, relative 5e-4, 16-cell groups | 106.57 | 50.65 | 55.92 | 5.1 |
| 4 | astra-max | balanced | RODAS3, relative 3e-3, absolute 1e-18, steps capped at 200 s | 109.29 | 48.52 | 60.77 | 5.0 |
| | MICM | L3 | three-stage Rosenbrock, relative 3e-4, absolute 1e-15 | 542.66 | 254.96 | 287.70 | 1 |

- **The held-out scenarios changed no rank.** The disclosed and held-out sums
  order the entries the same way at level 2 and at the top two places of
  level 3. The only close call is level 3's third place: sol-max's accurate
  point beats astra-max's balanced point by 2.6 % overall, behind it on the
  disclosed scenarios (50.65 against 48.52 s) and ahead on the held-out ones
  (55.92 against 60.77 s). That gap is inside the session's noise (below).
- **Margins.** The two level-2 winners sit just above the thresholds: 2.041
  digits overall on diurnal-5min, and isoprene at 1.035 (astra-xhigh/fast)
  and 1.064 (sol-max/fast) on diurnal-30min. The level-3 winner has 0.2
  digit overall (3.201 on clean-10min) and 0.39 per species (OClO 2.390 on
  diurnal-10min) to spare; sol-xhigh's L3, sol-max's accurate, and astra-max's
  balanced have 0.25 to 0.34 overall and 0.26 to 0.34 per species.
- **astra-max's fast3**, declared for level 3, misses it on the two held-out
  scenarios it missed at 1,000 cells (clean-48h-5min: 2.46 digits overall,
  phenol 1.57; winter60-10min: 2.94 overall), so astra-max's level-3 point is
  balanced, its recommended point, as the handoff anticipated.
- **MICM's level-3 baseline holds**: `micm/L3` reaches 2.053 digits for
  isoprene on diurnal-30min at 10,000 cells, as at 1,000, so the insurance
  setting `micm/L3b` (relative 2e-4, 551.06 s) was not needed. Both reach
  level 3 everywhere; the ranking takes the faster.
- **Why the speedups.** On diurnal-10min at 10,000 cells, MICM's three-stage
  method at relative 1e-2 takes about 9.8 internal steps per 10-minute
  interval, because its loop restarts every interval at one millionth of the
  step; astra-xhigh's fast point takes 1.6 (73,000 steps against 442,000),
  with 9 % of its attempts rejected. The rest of the factor is the Tier 1
  kernels. At level 3 the gap narrows because MICM's time hardly depends on
  its tolerance while the entries' does.
- **No point needed an exact rescore**: no score rounded onto a threshold.
- **Scores at 10,000 cells equal the 1,000-cell verification** within 0.01
  digit for every point, disclosed and held out, so the caveat about the
  last cell group (partial at 1,000 cells, full within the first 1,000 of
  10,000) made no difference.
- **Determinism:** every program produced byte-identical output on all five
  rounds of every scenario (the same SHA-256), so each was scored once per
  scenario.

Every point, with its lowest scores over the 11 scenarios (`ranking.md` has
the per-scenario medians; `timings.csv` per scenario has every run):

| Point | Level 2 | Level 3 | Lowest overall | Lowest species | Time (s) | Disclosed (s) | Held out (s) |
|---|---|---|---:|---|---:|---:|---:|
| sol-xhigh/L2 | yes | no | 2.456 | ISOP 1.445 | 54.52 | 25.79 | 28.74 |
| sol-xhigh/L2_margin | yes | no | 2.800 | ISOP 1.689 | 68.34 | 32.21 | 36.13 |
| sol-xhigh/L3 | yes | yes | 3.343 | ISOP 2.311 | 93.67 | 43.58 | 50.09 |
| sol-xhigh/recommended | yes | yes | 3.707 | ISOP 2.734 | 123.26 | 58.27 | 64.99 |
| sol-xhigh/L4 | yes | yes | 4.659 | ISOP 3.616 | 322.50 | 149.81 | 172.69 |
| sol-max/fast | yes | no | 2.041 | ISOP 1.064 | 39.80 | 18.79 | 21.01 |
| sol-max/balanced | yes | no | 2.290 | ISOP 1.386 | 44.75 | 21.12 | 23.63 |
| sol-max/accurate | yes | yes | 3.251 | ISOP 2.341 | 106.57 | 50.65 | 55.92 |
| sol-max/recommended | yes | yes | 3.571 | ISOP 2.554 | 137.36 | 65.36 | 72.00 |
| sol-max/precise | yes | yes | 4.251 | ISOP 3.069 | 218.38 | 103.30 | 115.08 |
| astra-xhigh/fast | yes | no | 2.041 | ISOP 1.035 | 34.72 | 16.15 | 18.58 |
| astra-xhigh/balanced | yes | no | 2.463 | PHENO 1.569 | 47.60 | 22.38 | 25.22 |
| astra-xhigh/accurate | yes | yes | 3.201 | OCLO 2.390 | 74.47 | 34.42 | 40.05 |
| astra-xhigh/recommended | yes | yes | 3.448 | BRO 2.550 | 86.32 | 40.20 | 46.12 |
| astra-xhigh/strict | yes | yes | 3.861 | ISOP 2.919 | 104.15 | 48.52 | 55.63 |
| astra-max/economy | yes | no | 2.267 | ISOP 1.271 | 62.13 | 29.02 | 33.12 |
| astra-max/fast3 | yes | no | 2.463 | PHENO 1.569 | 86.41 | 39.89 | 46.52 |
| astra-max/balanced | yes | yes | 3.274 | ISOP 2.263 | 109.29 | 48.52 | 60.77 |
| astra-max/conservative | yes | yes | 3.776 | ISOP 2.951 | 163.69 | 71.45 | 92.24 |
| micm/L2 | yes | no | 2.487 | ISOP 1.242 | 465.40 | 221.06 | 244.34 |
| micm/L3 | yes | yes | 3.288 | ISOP 2.053 | 542.66 | 254.96 | 287.70 |
| micm/L3b | yes | yes | 3.425 | ISOP 2.206 | 551.06 | 258.62 | 292.44 |

Digits overall per scenario (identical on every round):

| Point | clean-10min | diurnal-10min | diurnal-30min | diurnal-5min | polluted-10min | tier1-30s | clean-48h-5min | polluted-48h-30min | summer-15min | tropics-20min | winter60-10min |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| sol-xhigh/L2 | 2.81 | 2.84 | 2.54 | 3.30 | 2.84 | 5.74 | 2.47 | 2.46 | 2.61 | 2.53 | 2.67 |
| sol-xhigh/L2_margin | 3.15 | 3.10 | 2.80 | 3.44 | 3.04 | 5.85 | 3.59 | 2.82 | 2.85 | 2.88 | 2.98 |
| sol-xhigh/L3 | 3.53 | 3.50 | 3.62 | 3.70 | 3.34 | 5.93 | 3.83 | 3.43 | 3.39 | 3.62 | 3.71 |
| sol-xhigh/recommended | 3.84 | 3.80 | 4.02 | 3.87 | 3.71 | 5.96 | 4.00 | 3.94 | 3.93 | 4.06 | 4.13 |
| sol-xhigh/L4 | 4.88 | 4.90 | 5.23 | 4.66 | 5.00 | 7.99 | 4.72 | 5.26 | 5.12 | 5.23 | 5.48 |
| sol-max/fast | 2.42 | 2.51 | 2.12 | 2.04 | 2.40 | 5.06 | 2.30 | 2.21 | 2.34 | 2.18 | 2.28 |
| sol-max/balanced | 2.64 | 2.63 | 2.33 | 2.29 | 2.57 | 5.32 | 2.46 | 2.44 | 2.48 | 2.34 | 2.47 |
| sol-max/accurate | 3.30 | 3.27 | 3.70 | 3.41 | 3.34 | 5.82 | 3.54 | 3.74 | 3.25 | 3.61 | 3.67 |
| sol-max/recommended | 3.67 | 3.63 | 3.88 | 3.57 | 3.76 | 5.84 | 3.69 | 3.96 | 3.80 | 3.93 | 3.95 |
| sol-max/precise | 4.42 | 4.46 | 4.30 | 4.25 | 4.43 | 5.86 | 4.31 | 4.39 | 4.42 | 4.37 | 4.44 |
| astra-xhigh/fast | 2.56 | 2.53 | 2.10 | 2.04 | 2.40 | 5.03 | 2.33 | 2.21 | 2.30 | 2.19 | 2.27 |
| astra-xhigh/balanced | 2.77 | 2.73 | 2.74 | 3.16 | 2.65 | 5.49 | 2.46 | 2.78 | 2.53 | 2.54 | 2.67 |
| astra-xhigh/accurate | 3.20 | 3.21 | 3.67 | 3.45 | 3.93 | 6.08 | 3.66 | 3.56 | 3.93 | 3.68 | 3.94 |
| astra-xhigh/recommended | 4.24 | 4.29 | 3.83 | 3.45 | 4.14 | 6.33 | 3.82 | 3.67 | 4.04 | 4.02 | 4.12 |
| astra-xhigh/strict | 4.55 | 4.56 | 4.07 | 4.82 | 4.58 | 6.66 | 3.86 | 4.23 | 4.09 | 4.22 | 4.27 |
| astra-max/economy | 2.62 | 2.59 | 2.27 | 2.79 | 2.46 | 5.36 | 2.46 | 2.40 | 2.40 | 2.28 | 2.42 |
| astra-max/fast3 | 3.22 | 3.17 | 3.34 | 3.14 | 3.09 | 5.46 | 2.46 | 3.36 | 3.31 | 3.26 | 2.94 |
| astra-max/balanced | 3.55 | 3.47 | 3.61 | 3.45 | 3.43 | 5.46 | 3.29 | 3.65 | 3.59 | 3.54 | 3.27 |
| astra-max/conservative | 4.13 | 3.98 | 4.01 | 3.87 | 3.95 | 5.46 | 3.94 | 4.07 | 4.05 | 3.96 | 3.78 |
| micm/L2 | 3.42 | 3.43 | 2.49 | 3.51 | 3.29 | 5.55 | 3.68 | 2.59 | 3.26 | 2.79 | 3.16 |
| micm/L3 | 3.99 | 4.01 | 3.29 | 4.29 | 3.89 | 6.39 | 4.28 | 3.51 | 3.76 | 3.51 | 3.74 |
| micm/L3b | 4.03 | 4.05 | 3.42 | 4.38 | 3.98 | 6.42 | 4.40 | 3.63 | 3.86 | 3.59 | 3.82 |

Isoprene is the worst species for nearly every point on nearly every
scenario; the exceptions are BrO and NO on the short-step and 48-hour clean
scenarios, phenol on clean-48h-5min for the points near level 2, and the
halogen radicals OClO and BrO for astra-xhigh's extrapolation points.

## Procedure

- **Programs** (`programs.txt`): the 19 declared operating points of the four
  entries, built with clang from their collected branches in the
  re-verification of 2026-10-18 19:40 (identical to the committed reports),
  and MICM's three-stage Rosenbrock with 32-cell groups at the tolerances the
  scan below picked. Each scenario's `environment.txt` records every command
  line.
- **MICM's baselines** (`micm-scan/`): 18 tolerance settings on the 11
  scenarios at 1,000 cells, one round. The fastest that reaches level 2 on
  all 11 is relative 1e-2 with absolute 1e-15 (`micm/L2`); level 3, relative
  3e-4 with absolute 1e-15 (`micm/L3`, isoprene 2.05 digits, 0.05 above the
  threshold), with relative 2e-4 (`micm/L3b`) timed as well in case that
  margin did not hold at 10,000 cells. The ranking's MICM baseline at each
  level is the fastest of its programs that reaches the level at 10,000 cells.
- **The session** (`admin/time_tier3.sh`; `time_tier3.log`; `timings.csv`,
  `summary.md`, and `environment.txt` per scenario): the 22 programs on the
  10,000-cell versions of the 11 scenarios, 5 rounds in rotating order, 10
  threads, 5 s between runs, 1,210 runs from 21:19:19 on 2026-10-18 to
  05:46:43 on 2026-10-19, on the SSE's Apple M5 laptop (Mac17,2, 4
  performance and 6 efficiency cores, macOS 27.0.1), on AC power under
  `caffeinate`, with no agent, build, or other work of ours running. Each
  program's first output on a scenario was scored against the scenario's
  1,000-cell reference (`scenarios/score.py`: the first 1,000 cells scored,
  the gate applied to all 10,000); later outputs were scored again only when
  they differed by SHA-256, which never happened. The per-run solver outputs
  were deleted after scoring; the per-run logs stayed in the working
  directory.
- **Load audit** (`load-audit.log`, one sample a minute, 507 samples): AC
  attached throughout, the battery held at 80 % by the system. Processes
  above 5 % of a core other than the benchmark: macOS's media analysis
  daemon in 68 samples (mean 0.6 core), Spotlight's knowledge indexer in 57
  (0.5 core), its metadata store in 24 (0.5 core), XProtect's nightly scans
  in 8 samples (0.1 to 1 core), the PDF service in 7 samples between 21:21
  and 22:12 (indexing documents downloaded for the reference solution's
  literature review, deleted then), the scorer between runs, and this Claude
  Code session in 12 samples at a tenth of a core. The rotation spreads these over
  every program. The round-to-round spread of a point's elapsed time on a
  scenario was 1.3 to 3.9 % of the median on average and 8.4 % at most,
  against gaps of 14 % or more between adjacent ranks, except level 3's
  third place (2.6 %).
- **The ranking** (`ranking.md`, `admin/rank_tier3.py`): for each entry and
  level, its fastest point reaching the level on all 11 scenarios, by the sum
  over the scenarios of the point's median elapsed time; the sums over the
  disclosed and the held-out sets shown beside it. The script flags a score
  that rounds onto a threshold for an exact rescore; none did.
