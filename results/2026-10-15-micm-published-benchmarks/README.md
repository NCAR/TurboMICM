# MICM's published benchmarks, fetched 2026-10-15

MICM's continuous integration benchmarks every commit to `main` and publishes the
history as charts on GitHub Pages (https://ncar.github.io/micm/dev/bench/...).
`data/` holds the charts' data files as fetched from the `gh-pages` branch of
NCAR/micm on 2026-10-15; `latest.txt` is `summarize.py data`, the latest entry of
each series with its history. The latest CPU and GPU entries are for MICM commit
`e15898ce`, the commit TurboMICM uses.

| File | Series | Machine | Cells, steps |
|---|---|---|---|
| `ts1_timing.js` | TS1 wall clock, CPU, best of 3 | GitHub `ubuntu-latest` | 10,000, 30 |
| `ts1_instructions.js` | TS1 callgrind instruction counts, CPU | GitHub `ubuntu-latest` | 2,000, 5 |
| `timing.js`, `instructions.js` | the same for Chapman | GitHub `ubuntu-latest` | as above |
| `gpu_timing.js` | Chapman and TS1 wall clock, CUDA (and early CPU-backend entries) | CIRRUS A10 GPU runner | 10,000, 30 |
| `kokkos_timing.js` | Chapman and TS1 wall clock, Kokkos | CIRRUS A10 GPU runner | 10,000, 30 |
| `gpu_*instructions.js`, `gpu_ts1_timing.js` | early entries from 2026-08-11 | CIRRUS A10 GPU runner | as above |

Every series uses in-place Mozart LU and varies only the matrix layout (standard,
vector1, 2, 4, 8, 128), with all cells identical and one MICM state for all of
them. The wiki's MICM-Configurations page summarizes these results.
