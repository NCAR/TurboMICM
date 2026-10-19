# Fable's Tier 3 reference solution beside the contest's winners, 2026-10-19

The comparison session deferred on the morning of 2026-10-19 (the run log's
entry of 08:15): Fable's five declared operating points (`tier3/fable/max` at
`cb0f8a3`, not a contestant; wiki
[Fable's Tier 3 reference solution](https://github.com/NCAR/TurboMICM/wiki/Fable-Tier-3-Reference-Solution)),
the contest's two fastest level-2 points and two fastest level-3 points, and
MICM's two baselines, timed in one interleaved session under the ranking's
protocol ([`2026-10-19-tier3-ranking`](../2026-10-19-tier3-ranking/README.md)):
the 10,000-cell versions of the 11 scenarios, six disclosed and five held out,
at 10 threads, with three rounds instead of the ranking's five.

**The contest's ranking stands as published.** This session has only the
winners' points, three rounds, and a non-contestant; it places Fable's
solution beside the winners and changes nothing about the four entries'
order. The contest's points and MICM's baselines reproduced their
ranking-session sums within 2.2 %, and every score is identical to the
ranking's, as the outputs are deterministic.

## Results

Time is the sum over the 11 scenarios of the point's median elapsed time over
3 rounds at 10 threads, in seconds; the sums over the six disclosed and the
five held-out scenarios are shown beside it. Speedup is against MICM's
three-stage Rosenbrock at the ranking's baseline settings (`micm/L2`: relative
1e-2, absolute 1e-15; `micm/L3`: 3e-4 and 1e-15), with 32-cell groups and the
same photolysis update per step. The rows are in order of time; no rank is
assigned, since Fable's solution is not a contestant.

**Level 2** (2 digits overall, 1 per species, on all 11 scenarios):

| Entry | Point | Method | Time (s) | Disclosed (s) | Held out (s) | Speedup over MICM | Lowest overall | Lowest species |
|---|---|---|---:|---:|---:|---:|---:|---|
| Fable (not a contestant) | L2 | RODAS3, MICM's elementary controller, relative 5e-3, absolute 1e-20 | 24.55 | 11.07 | 13.48 | 18.9 | 2.273 | isoprene 1.160 |
| astra-xhigh (the contest's first) | fast | RODAS3, relative 1e-2, absolute 1e-18 | 34.10 | 15.86 | 18.24 | 13.6 | 2.041 | isoprene 1.035 |
| sol-max (the contest's second) | fast | RODAS3, relative 1e-2, 16-cell groups | 39.45 | 18.74 | 20.71 | 11.8 | 2.041 | isoprene 1.064 |
| MICM | L2 | three-stage Rosenbrock, relative 1e-2, absolute 1e-15 | 464.74 | 218.37 | 246.37 | 1 | 2.487 | isoprene 1.242 |

**Level 3** (3 digits overall, 2 per species, on all 11 scenarios):

| Entry | Point | Method | Time (s) | Disclosed (s) | Held out (s) | Speedup over MICM | Lowest overall | Lowest species |
|---|---|---|---:|---:|---:|---:|---:|---|
| Fable (not a contestant) | L3r | ROS3, H211b controller, relative 8e-3, absolute 1e-20 | 49.08 | 21.96 | 27.12 | 11.2 | 3.183 | isoprene 2.326 |
| astra-xhigh (the contest's first) | accurate | extrapolated linearly implicit Euler, order 4, relative 4e-3, absolute 1e-18 | 73.46 | 33.99 | 39.47 | 7.5 | 3.201 | OClO 2.390 |
| sol-xhigh (the contest's second) | L3 | RODAS3.1, relative 3e-4, absolute 1e-16 | 95.69 | 44.70 | 51.00 | 5.7 | 3.343 | isoprene 2.311 |
| MICM | L3 | three-stage Rosenbrock, relative 3e-4, absolute 1e-15 | 549.07 | 256.28 | 292.79 | 1 | 3.288 | isoprene 2.053 |

- **Fable's points are the fastest at both levels.** L2 is 1.39 times faster
  than astra-xhigh's fast point (34.10 against 24.55 s) and 1.61 times faster
  than sol-max's; L3r is 1.50 times faster than astra-xhigh's accurate point
  (73.46 against 49.08 s) and 1.95 times faster than sol-xhigh's L3. The
  disclosed and held-out sums give the same order. The provisional estimate
  from the verification runs (1.4 to 1.6 times the winners' speed at each
  level) held.
- **With more margin.** L2 clears the per-species threshold by 0.16 digit
  (isoprene 1.160 on winter60-10min) and the overall one by 0.27, where the
  contest's level-2 winners clear them by 0.04 to 0.06 (isoprene 1.035 and
  1.064 on diurnal-30min, 2.041 digits overall on diurnal-5min). L3r clears
  level 3 by 0.18 overall (3.183 on tropics-20min) and 0.33 per species
  (isoprene 2.326 on polluted-48h-30min), against 0.20 and 0.39 for
  astra-xhigh's accurate point.
- **Gaps against noise.** The gaps between adjacent places are 39 % and 16 %
  at level 2 and 50 % and 30 % at level 3. The round-to-round spread of a
  point's time on a scenario was 2.7 % of the median on average and 14.3 % at
  most; the five cells above 9 % are single slow rounds, four of them in the
  session's first 14 minutes, when the browser, a photo-analysis daemon, and
  Spotlight were active (see the audit), and the median of three rounds does
  not see a single slow round.
- **Why.** On diurnal-10min at 10,000 cells, Fable's L2 takes 1.9 internal
  steps per 10-minute interval (83,800 steps over 313 cell groups and 144
  intervals, 10 % rejected), more than astra-xhigh's fast point (1.6, 73,000
  steps, 9 % rejected), because its relative tolerance is tighter (5e-3
  against 1e-2); it is faster because its Tier 1 kernels cost less per step,
  23 µs of wall time per group step (one factorization and four solves) at
  10 threads against 36. At level 3, L3r's ROS3 steps and astra-xhigh's
  order-4 extrapolation steps come to nearly the same work on that scenario
  (236,600 factorizations and 710,000 solves against 248,200 and 620,000),
  and the factor of 1.5 is again the kernels' cost per operation. MICM takes
  9.8 steps per interval at either setting, because its loop restarts every
  interval at one millionth of the step.
- **Scores at 10,000 cells equal the 1,000-cell verification** within 0.01
  digit for every one of Fable's points, as for the contest's.
- **Determinism.** Every program produced byte-identical output on all three
  rounds of every scenario (121 program-scenario cells), so each was scored
  once per scenario, and no score rounded onto a threshold (`ranking.md`
  flags none).
- **Fable's other points** reach level 3 everywhere too: L3 (extrapolation,
  order 3) at 54.60 s, recommended (order 4) at 62.14 s, and accurate at
  71.91 s, which takes the time of astra-xhigh's accurate point (73.46 s)
  for 3.87 digits overall and 3.00 per species against 3.20 and 2.39.

Every point, with its lowest scores over the 11 scenarios (`ranking.md` has
the per-scenario medians; `timings.csv` per scenario has every run):

| Point | Level 2 | Level 3 | Lowest overall (scenario) | Lowest species (scenario) | Time (s) | Disclosed (s) | Held out (s) | Speedup over MICM at its level |
|---|---|---|---|---|---:|---:|---:|---:|
| fable-max/L2 | yes | no | 2.273 (winter60-10min) | isoprene 1.160 (winter60-10min) | 24.55 | 11.07 | 13.48 | 18.9 |
| astra-xhigh/fast | yes | no | 2.041 (diurnal-5min) | isoprene 1.035 (diurnal-30min) | 34.10 | 15.86 | 18.24 | 13.6 |
| sol-max/fast | yes | no | 2.041 (diurnal-5min) | isoprene 1.064 (diurnal-30min) | 39.45 | 18.74 | 20.71 | 11.8 |
| fable-max/L3r | yes | yes | 3.183 (tropics-20min) | isoprene 2.326 (polluted-48h-30min) | 49.08 | 21.96 | 27.12 | 11.2 |
| fable-max/L3 | yes | yes | 3.555 (polluted-48h-30min) | isoprene 2.471 (tropics-20min) | 54.60 | 25.09 | 29.51 | 10.1 |
| fable-max/recommended | yes | yes | 3.619 (clean-10min) | OClO 2.692 (diurnal-10min) | 62.14 | 28.88 | 33.27 | 8.8 |
| fable-max/accurate | yes | yes | 3.867 (diurnal-5min) | BrO 2.997 (diurnal-5min) | 71.91 | 32.89 | 39.02 | 7.6 |
| astra-xhigh/accurate | yes | yes | 3.201 (clean-10min) | OClO 2.390 (diurnal-10min) | 73.46 | 33.99 | 39.47 | 7.5 |
| sol-xhigh/L3 | yes | yes | 3.343 (polluted-10min) | isoprene 2.311 (polluted-48h-30min) | 95.69 | 44.70 | 51.00 | 5.7 |
| micm/L2 | yes | no | 2.487 (diurnal-30min) | isoprene 1.242 (diurnal-30min) | 464.74 | 218.37 | 246.37 | 1 |
| micm/L3 | yes | yes | 3.288 (diurnal-30min) | isoprene 2.053 (diurnal-30min) | 549.07 | 256.28 | 292.79 | 1 |

Digits overall per scenario (identical on every round):

| Point | clean-10min | diurnal-10min | diurnal-30min | diurnal-5min | polluted-10min | tier1-30s | clean-48h-5min | polluted-48h-30min | summer-15min | tropics-20min | winter60-10min |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| fable-max/L2 | 2.58 | 2.56 | 2.52 | 2.73 | 2.42 | 5.09 | 2.46 | 2.61 | 2.48 | 2.30 | 2.27 |
| astra-xhigh/fast | 2.56 | 2.53 | 2.10 | 2.04 | 2.40 | 5.03 | 2.33 | 2.21 | 2.30 | 2.19 | 2.27 |
| sol-max/fast | 2.42 | 2.51 | 2.12 | 2.04 | 2.40 | 5.06 | 2.30 | 2.21 | 2.34 | 2.18 | 2.28 |
| fable-max/L3r | 3.56 | 3.51 | 3.27 | 3.47 | 3.43 | 5.77 | 3.59 | 3.52 | 3.31 | 3.18 | 3.42 |
| fable-max/L3 | 3.61 | 3.84 | 3.65 | 4.12 | 3.95 | 5.87 | 3.68 | 3.56 | 3.74 | 3.66 | 3.73 |
| fable-max/recommended | 3.62 | 3.65 | 4.24 | 3.87 | 4.52 | 6.53 | 4.06 | 3.92 | 4.33 | 4.43 | 4.34 |
| fable-max/accurate | 4.81 | 4.72 | 4.57 | 3.87 | 4.54 | 6.84 | 4.27 | 4.48 | 4.46 | 4.55 | 4.67 |
| astra-xhigh/accurate | 3.20 | 3.21 | 3.67 | 3.45 | 3.93 | 6.08 | 3.66 | 3.56 | 3.93 | 3.68 | 3.94 |
| sol-xhigh/L3 | 3.53 | 3.50 | 3.62 | 3.70 | 3.34 | 5.93 | 3.83 | 3.43 | 3.39 | 3.62 | 3.71 |
| micm/L2 | 3.42 | 3.43 | 2.49 | 3.51 | 3.29 | 5.55 | 3.68 | 2.59 | 3.26 | 2.79 | 3.16 |
| micm/L3 | 3.99 | 4.01 | 3.29 | 4.29 | 3.89 | 6.39 | 4.28 | 3.51 | 3.76 | 3.51 | 3.74 |

Isoprene is the worst species for most points on most scenarios; the
exceptions are BrO and NO on the short-step and 48-hour clean scenarios,
phenol on clean-48h-5min for Fable's L2, HYDRALD and MVK for L3r, and the
halogen radicals OClO and BrO for the extrapolation points.

## Procedure

- **Programs** (`programs.txt`): Fable's five declared points, built with
  clang from `cb0f8a3` in the re-verification of 08:04 on 2026-10-19
  (`tier3-work/verify/fable-max-2/`); the contest's four points from the
  re-verification builds of 2026-10-18 19:40, the binaries the ranking
  session ran; and MICM's three-stage Rosenbrock with 32-cell groups at the
  ranking's baseline tolerances. Each scenario's `environment.txt` records
  every command line.
- **The session** (`admin/time_tier3.sh`; `time_tier3.log`; `timings.csv`,
  `summary.md`, and `environment.txt` per scenario): the 11 programs on the
  10,000-cell versions of the 11 scenarios, 3 rounds in rotating order, 10
  threads, 5 s between runs, 363 runs from 12:34:06 to 14:50:42 on
  2026-10-19, on the SSE's Apple M5 laptop (Mac17,2, 4 performance and 6
  efficiency cores, macOS 27.0.1), on AC power under `caffeinate`, with no
  agent, build, or other work of ours running. The laptop was on battery when
  the session was armed at 11:42, so a gate script waited for the adapter,
  attached at 12:33, and started the session a minute later. Each program's
  first output on a scenario was scored against the scenario's 1,000-cell
  reference (`scenarios/score.py`: the first 1,000 cells scored, the gate
  applied to all 10,000); later outputs would have been scored again only if
  their SHA-256 differed, which never happened. The per-run solver outputs
  were deleted after scoring.
- **Load audit** (`load-audit.log`, one sample a minute, 137 samples): AC
  attached in every sample; the battery charged from 75 % to full and was
  then let down to 80 % and held there by the system. Processes above 5 % of
  a core other than the benchmark and the scorer: macOS's media analysis
  daemon in 21 samples (mean 0.5 core, 2.2 at most), Spotlight's metadata
  store in 19 (0.5 core) and its workers, knowledge updater, and
  corespotlightd in 24 process entries (0.1 to 0.7 core), the window server
  in 15 (0.1 core), this Claude Code session and its helper processes in 13
  (0.2 core on average, 0.8 at most, during its two status checks), the
  photo analysis daemon in 2 (1 core, at 12:44 and 12:45), the SSE's browser
  in the first two minutes (0.3 core), and a dozen system services in one or
  two samples each at a tenth of a core. The rotation spreads these over
  every program.
- **The ranking** (`ranking.md`, `admin/rank_tier3.py`): for each entry and
  level, its fastest point reaching the level on all 11 scenarios, by the
  sum over the scenarios of the point's median elapsed time. Its level-2
  table lists sol-xhigh fourth by its L3 point, and its level-3 table has no
  sol-max row, because the session had only each entry's winning points; the
  contest's level-2 point of sol-xhigh and level-3 point of sol-max are in
  the ranking of 2026-10-19.
