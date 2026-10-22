# Tier 2 scenarios

Each directory here is a disclosed day-night scenario: TS1 run through 24 hours,
with its photolysis rates following a factor f(t).

    <name>/scenario.txt             the scenario file
    <name>/cells.csv                the cells: a CELLS_CSV, as in Tier 1
    <name>/photolysis.csv           the photolysis factor f(t)
    <name>/micm_final.csv           MICM's final state at the gate's settings
    <name>/micm_snapshots.csv       MICM's snapshots at the gate's settings
    <name>/reference_final.csv      a high-accuracy reference final state, for scale
    <name>/reference_snapshots.csv  the reference's snapshots
    <name>/micm.log                 how these were made, their scores, and their SHA-256

| Scenario | Cells | Steps | Initial NOx | Snapshots |
|---|---|---|---|---|
| diurnal-30min | 1,000, seed 101 | 48 of 30 min | x 1 | every 6 h |
| diurnal-10min | 1,000, seed 101 | 144 of 10 min | x 1 | every 6 h |
| diurnal-5min | 1,000, seed 101 | 288 of 5 min | x 1 | every 6 h |
| polluted-10min | 1,000, seed 102 | 144 of 10 min | x 10 | every 6 h |
| clean-10min | 1,000, seed 103 | 144 of 10 min | x 0.1 | every 6 h |

- Each runs for 24 hours from 06:00 local solar time, at 40 N on day 80 (the
  March equinox). diurnal-30min, diurnal-10min, and diurnal-5min share their
  cells and f(t).
- The gate also uses scenarios that are not disclosed: other seeds, latitudes,
  seasons, start times, and initial NOx, steps of 5 to 30 minutes, and runs of
  up to 48 hours, all made with make_scenario.py and recording a snapshot every
  6 hours.

## Scenario file

`key = value` lines; `#` starts a comment. Paths are relative to the directory
of the scenario file.

    cells = cells.csv            the cells (CELLS_CSV)
    dt = 600                     the step length, s
    steps = 144                  the number of steps
    photolysis = photolysis.csv  the photolysis factor f(t)
    snapshot_every = 36          optional: record the state every 36 steps

## Running a scenario

- Every cell starts from TS1's InitState, then takes its CELLS_CSV row
  (temperature, pressure, species, and rate parameters), and its air density is
  recomputed from temperature and pressure, as in Tier 1 (see
  ROOT/reference/README.md).
- Step k, for k = 0 to steps - 1, runs from t = k dt to t = (k + 1) dt. At its
  start, every photolysis rate parameter (a label beginning with "j") of every
  cell is set to its CELLS_CSV value times f(k dt), and it stays constant
  through the step. Temperature, pressure, and all other rate parameters stay
  constant throughout. There is no warm-up step.
- photolysis.csv has the header `time_s,factor` and rows of increasing time. f
  is linear between rows, and constant before the first row and after the last.
  The rows are 300 s apart, and f(t) = max(0, cos SZA(t)) / cos SZA(noon), from
  the solar zenith angle SZA; so the photolysis rates in cells.csv are noon
  values, and f = 0 at night.

## Output

Both output files are CSV, with every value written to 17 significant digits
(`%.17g`). Species may come in any order; score.py matches them by name.
- ts1_final.csv: the state after the last step, as in Tier 1. The header is
  `temperature,pressure,` followed by the species names; then one row per cell,
  in CELLS_CSV order.
- ts1_snapshots.csv, for a scenario with snapshot_every = M: the header is
  `step,cell,temperature,pressure,` followed by the species names. Then come
  rows for each of the steps M, 2M, ..., and always the last step, with one row
  per cell for each. `step` counts the steps completed, and `cell` is the
  0-based CELLS_CSV row.

## The gate's settings

MICM's results here are its three-stage Rosenbrock, Tier 1's solver (standard
ordering, in-place Mozart LU, MICM's parameters, error estimate, and step-size
control), at relative tolerance 1e-6 (MICM's default) and absolute tolerance
1e-15 mol m-3 for every species, with 32-cell groups:

    micm_scenario --group-size 32 --atol 1e-15 SCENARIO

The groups are consecutive cells in CELLS_CSV order, cells 0-31, 32-63, and so
on, the last group holding what is left (8 cells of 1,000). Each group is one
MICM state, solved step by step with its own step-size control, steered by its
own error norm; so the grouping changes the result, and the gate's runs use the
same grouping as MICM's.

On every scenario, disclosed and undisclosed, the gate runs

    bin/ts1_scenario_specialized --rtol 1e-6 --atol 1e-15 [--group-size 32] SCENARIO

(`--group-size 32` unless the group size is fixed at build time, `make
GROUP_SIZE=32 bench`) and scores its snapshots against MICM's:

    python3 ROOT/scenarios/score.py ts1_snapshots.csv ROOT/scenarios/<name>/micm_snapshots.csv --level 10

which passes with no NaN, no infinity, and no value below -1e-20 mol m-3, at
least 10 significant digits overall, and at least 9 for every species.

For scale: MICM built with clang and with gcc-15 takes exactly the same internal
steps on every scenario, none rejected, and agrees with itself to 14 digits
overall and at least 13.47 for every species, so roundoff alone leaves more than
three digits of margin. MICM at the gate's settings is itself 4.75 to 5.40
digits from the reference solutions:

| Scenario | MICM's internal steps | Agreement, clang with gcc-15: overall (worst species) | MICM against the reference |
|---|---:|---|---:|
| diurnal-30min | 48,619 | 14.37 (ISOP 13.81) | 4.75 |
| diurnal-10min | 90,980 | 14.28 (ISOP 13.70) | 5.05 |
| diurnal-5min | 145,891 | 14.24 (ISOP 13.69) | 5.21 |
| polluted-10min | 106,393 | 14.07 (ISOP 13.47) | 5.40 |
| clean-10min | 84,273 | 14.35 (ISOP 13.84) | 4.91 |

MICM's default absolute tolerance, 1e-3 mol m-3 (about 24 ppm), is above every
trace gas, and with it the error control does nothing useful over a day. With
32-cell groups, MICM's default stops in step 1 on diurnal-30min and
polluted-10min ("Convergence Exceeded Max Steps"), finishes diurnal-10min with
no correct digit (-7.63), and keeps 1.21 digits on diurnal-5min (worst species
0.21) and 2.96 on clean-10min (worst 2.16).

## Scoring

    python3 ROOT/scenarios/score.py RESULT REFERENCE [--threshold 1e-14] [--level L] [--json]

score.py compares two ts1_snapshots.csv files, or two ts1_final.csv files.
- It prints the overall agreement in significant digits, the five worst
  species, and the gate.
- With --level L, it also checks for L digits overall and L - 1 for every
  species.
- Its exit status is 0 only on a pass.
- Its docstring gives the definitions: following Sandu et al. (1997), a value
  counts when the reference value is at least the threshold (1e-14 mol m-3).

## Reference solutions

Each reference is MICM's six-stage Rosenbrock method (rosenbrock6dae) at
relative tolerance 1e-10 and absolute tolerance 1e-23 mol m-3, checked against
the same method at 1e-9 and 1e-22, with one group of 100 cells per thread on 10
threads. micm.log gives that agreement, the reference's convergence:

| Scenario | Reference convergence: overall (worst species) |
|---|---|
| diurnal-30min | 10.31 (OCLO 8.86) |
| diurnal-10min | 10.22 (CLO 9.34) |
| diurnal-5min | 10.18 (BRO 9.30) |
| polluted-10min | 10.25 (MCO3 9.40) |
| clean-10min | 10.20 (CLO 9.37) |

## Tools

- micm_scenario.cpp runs a scenario with MICM, for the gate's results,
  references, and baselines; its header comment gives the options. Build it
  with:

      clang++ -std=c++20 -O3 -DNDEBUG \
          -I ROOT/micm/include -I ROOT/micm/test/benchmark \
          ROOT/scenarios/micm_scenario.cpp -o bin/micm_scenario

- make_scenario.py writes a scenario, and its docstring gives the options. These
  commands made the disclosed scenarios:

      make_scenario.py diurnal-30min --cells 1000 --seed 101 --dt 1800 --steps 48 --diurnal --snapshot-every 12
      make_scenario.py diurnal-10min --cells 1000 --seed 101 --dt 600 --steps 144 --diurnal --snapshot-every 36
      make_scenario.py diurnal-5min --cells 1000 --seed 101 --dt 300 --steps 288 --diurnal --snapshot-every 72
      make_scenario.py polluted-10min --cells 1000 --seed 102 --dt 600 --steps 144 --diurnal --nox-scale 10 --snapshot-every 36
      make_scenario.py clean-10min --cells 1000 --seed 103 --dt 600 --steps 144 --diurnal --nox-scale 0.1 --snapshot-every 36

  Each also uses the defaults --lat 40 --doy 80 --start 6.
