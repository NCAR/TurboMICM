# Experiment records

Each CSV row is one native run, including failures. `options` is a JSON object of overrides relative to the defaults in `DESIGN.md`; configuration JSON files are the complete scan inputs. `cells` is the number actually integrated. `repeat` is zero-based. `elapsed_ms` is printed by the C or independent MICM C++ executable; Python does not measure solver speed. Group-wide operation counts retain their documented meaning. Accuracy and the worst species come from invoking the unmodified supplied `../scenarios/score.py --json`. Its printed values are rounded to three decimals. Failed runs have no accuracy; missing fields mean unavailable, not zero.

- `families-prefix128.csv`: first 128 cells, all six scenarios, same new-method controls/tolerances; vary Rosenbrock formula or BDF maximum order. Legacy defaults provide separate anchors.
- `precision-prefix128.csv`: vary relative tolerance, absolute tolerance, or BDF matrix age. The final two BDF configurations overlapped a compiler validation job; their compiler identity is marked ambiguous and their timings are excluded from performance conclusions. No declared point depends on them.
- `structure-prefix128.csv`: isolate cell grouping, norm, first-step proposal, controller, Jacobian age, and step holding. Combined holding/age rows test the interaction after the individual choices. ROS2 age rows examine the W-method alternative.
- `controller-prefix128.csv`: controlled growth, safety, H211b, maximum-norm tolerance, and neighboring Rosenbrock work–precision checks. Source and binary hashes are in its metadata file.
- `full-precision.csv`: all 1,000 cells in every disclosed scenario, tolerance sweep and inherited C defaults.
- `operating-points-clang.csv`: final five points, all disclosed cells and snapshots, five serial native repetitions per point/scenario. This is the source for reported C timings and work–precision figures.
- `micm-baselines.csv`: preliminary independent C++ comparisons, including failure diagnostics.
- `micm-timings.csv`: five repetitions of both independent three-stage baselines at the same 1,000-cell scale and group size. The accuracy-controlled baseline uses `rtol=3e-4`, `atol=1e-18`.

The earliest two scans used the code captured in commit `24499eb`, without projection instrumentation. The next structure scan added step holding and moved endpoint validation outside the clock. Later metadata records source-file hashes, binary hash, compiler/flags, command, and configuration hash. All other prefix timings are screening measurements, not final speed estimates; the 128-cell input does not occupy all ten workers. Empty projection counters in older rows mean that instrumentation did not exist yet. No rows from those versions have been replaced with later results.

Native output files are overwritten in `results/scan-run/` to keep storage bounded. Input, reference, and output hashes are retained in every successful row. Prefix input/reference files are reproducibly constructed in `results/prefix-N/`; species and cells are never selected according to their errors. A study can resume by skipping existing `(configuration, scenario, repeat)` tuples. Use a fresh study name after changing numerical code, compiler, cell subset, or options; mixing them in one CSV is not supported.

Example:

```sh
make bench CC=clang
.conda/bin/python tools/scan.py scans/families.json --study my-families --cells 128
.conda/bin/python tools/scan.py scans/operating-points.json --study my-final --repeats 5
.conda/bin/python tools/scan.py scans/micm-baselines.json --study my-micm \
  --backend micm --compiler clang++ --repeats 5
.conda/bin/python tools/plot_scans.py
```

The plotting script uses the committed final study names. All timings in it come from native executables, with unchanged `OPTFLAGS`. Its figures show both overall and worst-species accuracy because the ranking requires both.
