# Experiment records

Each CSV row is one native run, including failures. `options` is a JSON object of overrides relative to the defaults in `DESIGN.md`; configuration JSON files are the complete scan inputs. `cells` is the number actually integrated. `repeat` is zero-based. `elapsed_ms` is printed by the C or independent MICM C++ executable; Python does not measure solver speed. Group-wide operation counts retain their documented meaning. Accuracy and the worst species come from invoking the unmodified supplied `../scenarios/score.py --json`. Its printed values are rounded to three decimals. Failed runs have no accuracy; missing fields mean unavailable, not zero.

- `families-prefix128.csv`: first 128 cells, all six scenarios, same new-method controls/tolerances; vary Rosenbrock formula or BDF maximum order. Legacy defaults provide separate anchors.
- `precision-prefix128.csv`: vary relative tolerance, absolute tolerance, or BDF matrix age. The final two BDF configurations overlapped a compiler validation job; their compiler identity is marked ambiguous and their timings are excluded from performance conclusions. No declared point depends on them.
- `structure-prefix128.csv`: isolate cell grouping, norm, first-step proposal, controller, Jacobian age, and step holding. Combined holding/age rows test the interaction after the individual choices. ROS2 age rows examine the W-method alternative.
- `controller-prefix128.csv`: controlled growth, safety, H211b, maximum-norm tolerance, and neighboring Rosenbrock work–precision checks. Source and binary hashes are in its metadata file.
- `full-precision.csv`: all 1,000 cells in every disclosed scenario, tolerance sweep and inherited C defaults.
- `operating-points-clang.csv`: earlier 32-cell static controls for the five tolerance points, all disclosed cells and snapshots, five serial native repetitions per point/scenario. The selected final points use the later `final-clang-*` records below.
- `micm-baselines.csv`: preliminary independent C++ comparisons, including failure diagnostics.
- `micm-timings.csv`: five repetitions of both independent three-stage baselines at the same 1,000-cell scale and group size. The accuracy-controlled baseline uses `rtol=3e-4`, `atol=1e-18`.
- `large-structure.csv`: 10,000-cell group-size and scheduling screen, with the supplied reference covering the first 1,000 cells. The gate checks all output cells.
- `throughput-validation.csv`: full disclosed validation for 16-cell dynamic groups, with species diagnostics in the corresponding `-species.csv` file.
- `final-clang-1000.csv` and `final-clang-10000.csv`: five-repeat final operating points with 16-cell dynamic groups. The 1,000-cell study has full reference coverage; the 10,000-cell study has reference coverage for its first 1,000 cells and gate coverage for every output cell.
- `stress-reference-generation.csv`: tight and looser independent MICM6 runs. The tight run scores against itself (16 digits, a setup check); only the looser-vs-tight comparison establishes convergence. All three comparisons exceed 10 overall digits and 9 worst-species digits.
- `stress-final.csv`: final operating points on three extra 64-cell, 48-hour scenarios; `stress-48h.csv` is the earlier 32-cell-group comparison. Both have independent tight references and species diagnostics.
- `micm-large.csv`: independent three-stage baselines on the large extensions, with the same limited reference coverage and complete output gate.

The earliest two scans used the code captured in commit `24499eb`, without projection instrumentation. The next structure scan added step holding and moved endpoint validation outside the clock. Later metadata records source-file hashes, binary hash, compiler/flags, command, and configuration hash. All other prefix timings are screening measurements, not final speed estimates; the 128-cell input does not occupy all ten workers. Empty projection counters in older rows mean that instrumentation did not exist yet. No rows from those versions have been replaced with later results.

Native output files are overwritten in `results/scan-run/` to keep storage bounded. Input, reference, and output hashes are retained in every successful row. Prefix input/reference files are reproducibly constructed in `results/prefix-N/`; species and cells are never selected according to their errors. A study can resume by skipping existing `(configuration, scenario, repeat)` tuples. Use a fresh study name after changing numerical code, compiler, cell subset, or options; mixing them in one CSV is not supported.

Example:

```sh
make bench CC=clang
.conda/bin/python tools/scan.py scans/families.json --study my-families --cells 128
.conda/bin/python tools/scan.py scans/operating-points.json --study my-static-control --repeats 5
.conda/bin/python tools/scan.py scans/operating-points-throughput.json --study my-throughput --repeats 5
.conda/bin/python tools/scan.py scans/micm-baselines.json --study my-micm \
  --backend micm --compiler clang++ --repeats 5
.conda/bin/python tools/plot_scans.py
```

The plotting script uses the committed final study names. All timings in it come from native executables, with unchanged `OPTFLAGS`. Its figures show both overall and worst-species accuracy because the ranking requires both.

Large and additional scenarios can be recreated without writing to the supplied inputs:

```sh
.conda/bin/python tools/large_scenarios.py
.conda/bin/python tools/scan.py scans/operating-points-throughput.json --study my-large \
  --scenario-root results/large --repeats 5
.conda/bin/python tools/stress_references.py
.conda/bin/python tools/scan.py scans/operating-points-throughput.json --study my-stress \
  --scenario-root results/stress \
  --scenarios summer-48h-30min,clean-winter-48h-5min,polluted-48h-10min --metrics
```
