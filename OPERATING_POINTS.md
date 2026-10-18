# Declared operating points

These four configurations solve the complete mechanism in double precision.
All use four-stage RODAS3 with the cell-wise maximum RMS norm, integral step
control, a fresh analytical Jacobian per trial, and a step prediction carried
between external calls. Only the explicit settings below differ from the new
method defaults in [DESIGN.md](DESIGN.md). The program's default method remains
the original MICM three-stage Rosenbrock.

| Point | rtol | atol (mol m-3) | Maximum internal step (s) | Cells per group | Expected level |
| --- | ---: | ---: | ---: | ---: | --- |
| economy | 0.003 | 1e-16 | External interval | 32 | L=2: overall ≥2, every species ≥1 |
| fast3 | 0.003 | 1e-18 | 300 | 32 | L=3: overall ≥3, every species ≥2 |
| **balanced — recommended** | **0.003** | **1e-18** | **200** | **32** | **L=3 with a larger measured margin** |
| conservative | 0.003 | 1e-18 | 100 | 32 | L=3 with additional accuracy |

Use the following commands, replacing `SCENARIO` with a scenario file. Each
writes `ts1_final.csv` and `ts1_snapshots.csv` in its current working directory.
All cores are used by default; `--threads N` can limit the run. The measurements
in [RESULTS.md](RESULTS.md) request ten threads explicitly.

```sh
bin/ts1_scenario --method rodas3 --option rtol=0.003 --option atol=1e-16 --option h_max=0 --option group_size=32 SCENARIO

bin/ts1_scenario --method rodas3 --option rtol=0.003 --option atol=1e-18 --option h_max=300 --option group_size=32 SCENARIO

bin/ts1_scenario --method rodas3 --option rtol=0.003 --option atol=1e-18 --option h_max=200 --option group_size=32 SCENARIO

bin/ts1_scenario --method rodas3 --option rtol=0.003 --option atol=1e-18 --option h_max=100 --option group_size=32 SCENARIO
```

The machine-readable equivalents are in
[scans/operating-points.json](scans/operating-points.json). Point names are report
labels; they are not additional `--method` aliases.

Balanced is recommended for an atmospheric model because the 200-second cap
improves the accumulated error while retaining much of RODAS3's cost advantage.
Fast3 approaches the three-digit acceptance boundary on some disclosed and
independent cases. Balanced leaves more room for changes in conditions, at a
moderate cost. Conservative spends additional work on accuracy but does not
establish a universal four-digit setting. Economy is appropriate only when the
model's chemistry error budget allows the two-digit target.

Across the six disclosed and three independent cases, balanced achieves at
least **3.433 overall digits and 2.490 digits for every species**. Fast3's
corresponding minima are 3.086 and 2.057. These measurements support the choice
of margin rather than treating a local tolerance as a global error guarantee.

The caps restrict internal integration steps; they do not alter the scenario's
photolysis schedule or other chemistry. No species-specific tolerances, reaction
removal, precision changes, Jacobian lagging, or approximate factor reuse are
used in these points. Every accepted state is finite and nonnegative; an
interval that cannot satisfy the solver checks returns failure explicitly.

Accuracy expectations are supported by the six disclosed cases and three
independent 48-hour cases with both compilers. They remain empirical claims
for unseen conditions. [ACCURACY.md](ACCURACY.md) states the metric, projection
limits and scope; [RESULTS.md](RESULTS.md) gives actual digits, projection sizes,
and five-run 10,000-cell timings, with prefix-only accuracy clearly labeled.
