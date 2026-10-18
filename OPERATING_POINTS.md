# Operating points

Five runtime configurations are declared below. `operating_points.json` is the
machine-readable manifest used by the timing and validation tools. All use
32-cell groups and guided OpenMP scheduling with minimum chunk 1. Leave off
`--threads` to use all available cores; `--threads N` limits the run.

| Point | Method | rtol | atol (mol/m3) | Expected disclosed accuracy, overall / worst species |
|---|---|---:|---:|---:|
| L2 | RODAS3.1 | 0.002 | 1e-18 | 2.54 / 1.44 digits |
| L2_margin | RODAS3.1 | 0.001 | 1e-18 | 2.79 / 1.68 digits |
| L3 | RODAS3.1 | 0.0003 | 1e-16 | 3.34 / 2.32 digits |
| **recommended** | **RODAS3.1** | **0.0003** | **1e-18** | **3.70 / 2.73 digits** |
| L4 | six-stage RODAS4 | 0.0001 | 1e-18 | at least 4.5 / 3.5 digits |

The recommended point targets the 0.1% overall / 1% per-species integration
budget in ACCURACY.md, with useful margin. It uses a current analytic Jacobian
at every accepted step, positivity-aware rejection, an I controller and warm
step proposals. Its stricter absolute tolerance protects low concentrations
and their influence on the rest of the mechanism. The faster points trade
that margin for throughput; L4 provides an additional precision option.

All other defaults are those in DESIGN.md: no stale-Jacobian approximation or
step deadband, exact forcing reuse only on retries, start 0.01 s, safety 0.9,
growth limit 6, shrink limit 0.1, max-cell RMS error, and maximum 100000 attempts
per external step. No species or reaction is removed. Accuracy expectations
are measured on the disclosed scenarios, not guarantees about unobserved
conditions; RESULTS.md gives the independent 48-hour checks and timing evidence.

Exact commands (replace `SCENARIO` with its scenario.txt path):

L2:

```sh
bin/ts1_scenario --method rodas31 --option rtol=0.002 --option atol=1e-18 --option group_size=32 --option schedule=2 --option schedule_chunk=1 SCENARIO
```

L2_margin:

```sh
bin/ts1_scenario --method rodas31 --option rtol=0.001 --option atol=1e-18 --option group_size=32 --option schedule=2 --option schedule_chunk=1 SCENARIO
```

L3:

```sh
bin/ts1_scenario --method rodas31 --option rtol=0.0003 --option atol=1e-16 --option group_size=32 --option schedule=2 --option schedule_chunk=1 SCENARIO
```

recommended:

```sh
bin/ts1_scenario --method rodas31 --option rtol=0.0003 --option atol=1e-18 --option group_size=32 --option schedule=2 --option schedule_chunk=1 SCENARIO
```

L4:

```sh
bin/ts1_scenario --method rodas4 --option rtol=0.0001 --option atol=1e-18 --option group_size=32 --option schedule=2 --option schedule_chunk=1 SCENARIO
```

Reproduce five sequential, interleaved timing repetitions with:

```sh
.conda/bin/python tools/operating_points.py --repeat 5
```
