# Tier 2 scenarios

Each directory here is a disclosed scenario: TS1 run through time, with its
photolysis rates following a factor f(t).

    <name>/scenario.txt             the scenario file
    <name>/cells.csv                the cells: a CELLS_CSV, as in Tier 1
    <name>/photolysis.csv           the photolysis factor f(t)
    <name>/reference_final.csv      the reference final state
    <name>/reference_snapshots.csv  the reference snapshots, for a scenario with snapshots
    <name>/references.log           how the reference was made and checked

| Scenario | Cells | Steps | Photolysis | Initial NOx | Snapshots |
|---|---|---|---|---|---|
| tier1-30s | 1,000, seed 2026 | 101 of 30 s | constant | x 1 | none: final state only |
| diurnal-30min | 1,000, seed 101 | 48 of 30 min | diurnal | x 1 | every 6 h |
| diurnal-10min | 1,000, seed 101 | 144 of 10 min | diurnal | x 1 | every 6 h |
| diurnal-5min | 1,000, seed 101 | 288 of 5 min | diurnal | x 1 | every 6 h |
| polluted-10min | 1,000, seed 102 | 144 of 10 min | diurnal | x 10 | every 6 h |
| clean-10min | 1,000, seed 103 | 144 of 10 min | diurnal | x 0.1 | every 6 h |

- tier1-30s is Tier 1's benchmark on the first 1,000 of its seed-2026 cells. Its
  101 steps are Tier 1's untimed warm-up step and its 100 timed steps, so its
  final state is MICM's Tier 1 final state for those cells.
- The diurnal scenarios run for 24 hours from 06:00 local solar time, at 40 N on
  day 80 (the March equinox). diurnal-30min, diurnal-10min, and diurnal-5min
  share their cells and f(t).
- The ranking times 10,000-cell scenarios, and it also uses scenarios that are
  not disclosed: other seeds, latitudes, seasons, and 48-hour runs, all made
  with make_scenario.py.

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
  - For tier1-30s, f = 1.
  - For the diurnal scenarios, the rows are 300 s apart, and
    f(t) = max(0, cos SZA(t)) / cos SZA(noon), from the solar zenith angle SZA.
    So the photolysis rates in cells.csv are noon values, and f = 0 at night.

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

## Scoring

    python3 ROOT/scenarios/score.py RESULT REFERENCE [--threshold 1e-14] [--level L] [--json]

score.py compares ts1_snapshots.csv with reference_snapshots.csv. For a scenario
without snapshots, it compares ts1_final.csv with reference_final.csv instead.
- It prints the overall accuracy in significant digits, the five worst species,
  and the gate.
- With --level L, it also checks for L digits overall and L - 1 for every
  species.
- Its exit status is 0 only on a pass.
- Its docstring gives the definitions.

For example:

    cd <a work directory>
    ROOT/<run>/bin/ts1_scenario ROOT/scenarios/diurnal-10min/scenario.txt
    python3 ROOT/scenarios/score.py ts1_snapshots.csv \
        ROOT/scenarios/diurnal-10min/reference_snapshots.csv --level 3

## Reference solutions

Each reference is MICM's six-stage Rosenbrock method (rosenbrock6dae) at
relative tolerance 1e-10 and absolute tolerance 1e-23 mol m-3. It is checked
against the same method at 1e-9 and 1e-22. references.log gives that agreement,
and the accuracy of MICM's default three-stage Rosenbrock, both in score.py's
digits.

## Tools

- micm_scenario.cpp runs a scenario with MICM, for references and baselines; its
  header comment gives the options. Build it with:

      SDKROOT=$(xcrun --sdk macosx --show-sdk-path) clang++ -std=c++20 -O3 -DNDEBUG \
          -I ROOT/micm/include -I ROOT/micm/test/benchmark \
          ROOT/scenarios/micm_scenario.cpp -o bin/micm_scenario

- make_scenario.py writes a scenario, and its docstring gives the options. These
  commands made the disclosed scenarios:

      make_scenario.py tier1-30s --cells 1000 --seed 2026 --dt 30 --steps 101 --constant
      make_scenario.py diurnal-30min --cells 1000 --seed 101 --dt 1800 --steps 48 --diurnal --snapshot-every 12
      make_scenario.py diurnal-10min --cells 1000 --seed 101 --dt 600 --steps 144 --diurnal --snapshot-every 36
      make_scenario.py diurnal-5min --cells 1000 --seed 101 --dt 300 --steps 288 --diurnal --snapshot-every 72
      make_scenario.py polluted-10min --cells 1000 --seed 102 --dt 600 --steps 144 --diurnal --nox-scale 10 --snapshot-every 36
      make_scenario.py clean-10min --cells 1000 --seed 103 --dt 600 --steps 144 --diurnal --nox-scale 0.1 --snapshot-every 36

  Each diurnal scenario also uses the defaults --lat 40 --doy 80 --start 6. With
  --cells 10000 and the other arguments unchanged, a command makes the
  10,000-cell version, whose first 1,000 cells are the disclosed ones.
