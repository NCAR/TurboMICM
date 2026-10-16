# Validation records

`clang-cases.txt` and `gcc-cases.txt` contain the final successful result for
each exact ID in `tests/cases.txt`, in the supplied order. They were extracted
from `make test-cpu CC=clang` and `make test-cpu CC=gcc-15`, after checking for
missing, duplicated, or failing IDs. Both runs also passed the additional native
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
