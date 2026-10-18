# Run records

Each row records the method/options, scenario, thread count, compiler, C time,
solver operation counts and the unmodified supplied scorer's overall and worst
species digits. Counts sum over independent groups. The `configuration` JSON
specifies deviations from the documented defaults; an absent option has its
DESIGN.md default. Failed runs remain in the tables. Nonlinear iterations are
zero for these linearly implicit methods.

The numerical engine in initial family scans is commit 0388666; later commits
add independent tests, documentation and experiment orchestration. Builds use
`-O3 -DNDEBUG`. The actual build stamp is `build/configuration.txt`. All initial
scans use 1,000-cell disclosed inputs and ten threads, groups of 32 (including
a final group of eight), unless a row states an override.

`smoke.csv` is a diagnostic first run overlapping a test build. Initial
`baseline.csv` and the first four rows of `family.csv` overlapped and shared an
output path due to an orchestration error; their accuracy and time are
provisional and must not be used for conclusions. Isolated corrected trials
supersede them. Other structure-scan times are development measurements;
operating-point repeated timings are the basis for final performance claims.
The scan tool now uses a separate work directory per CSV. Do not launch two
instances with the same output CSV simultaneously.

Reproduce the family scan with:

```
.conda/bin/python tools/scan.py --suite family --output scans/family.csv --resume
```

Available suites are `baseline`, `family`, `loose`, `controller`, and `reuse`.
`--method` and repeated `--option` specify an additional single configuration.
`--repeat N` repeats each scenario. `--scenario NAME` can restrict a scan.
Every run invokes the compiled C executable for time/counters and the supplied
Python scorer only after C execution has ended. Python never runs in a timed
integration. `tools/plot_scans.py` plots the resulting work–precision curves.
