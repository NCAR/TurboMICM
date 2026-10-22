# Validation records

## Current mechanism compiler

`clang-tier2.json` and `gcc-tier2.json` record the clean build times, generated
source hashes and size, six generator parity tests, small MICM gates, full
10,000-cell timings/gates, and final/snapshot scenario scores. Every timed run
has its own gate result. Both compiler runs check all five supplied scenarios
and additional 48-hour/interpolation scenarios. Reference-solution scores are
included for each supplied scenario.

`clang-cases.txt` and `gcc-cases.txt` now contain the 2026-10-22 runs of all
680 unchanged CPU cases. `parameter-checks.json` records alternate generator
settings. `specialization-experiments.json` contains the controlled tuning sweep,
including raw samples, kernel profiles, and every setting's accuracy result.
`scenario-format-checks.json` checks absent and oversized snapshot intervals:
both produce all 33 final cell rows at step eight and match MICM's final state.

Run `.conda/bin/python tools/validate_tier2.py` to reproduce the full main
validation. It creates inputs and run directories under `runs/tier2`, builds
the provided MICM programs in this repository, and never writes to sibling
input directories. The script runs builds, tests, and timing workloads
sequentially. The 18-thread timings have three samples; one-thread timings and
scenario timings have one. Full logs and CSVs are retained in ignored run
directories; the committed records suffice to inspect the results.

## Retained Tier 1 records

The per-case files retain the same format as Tier 1: one result per exact ID in
`tests/cases.txt`, in supplied order. Both current runs also pass the native
controller/linear-solver extension check.

`benchmarks.json` records repeated runs of 10,000 cells, seed 2026, 100 timed
steps of 30 seconds, following one warm-up step. Each record retains the
executable's elapsed time and counters. Native runs also retain the supplied
`compare_final.py` result against MICM for every final-state value.

All programs were built with `-O3 -DNDEBUG`. `clang` and `gcc` labels denote the
C port; `micm` denotes the supplied C++ reference built by `make micm-ref` with
clang++. `*-one` uses `--threads 1`; the other records use the runtime default
(10 cores on this machine). C uses groups of 32, while the reference uses its
default of one group per worker. These grouping differences explain the
different summed statistics. Results and reproduction commands are discussed
in the top-level `RESULTS.md` and `README.md`.

Full CSVs, build products, intermediate sweeps, and raw build logs remain under
the ignored `runs`, `build`, `bin`, and `lib` directories. They are not required
to rebuild or run the C port.
