# Validation records

## Final Tier 3 evidence

- `tier3-gates.json` records both compilers' 680-case results and the retained
  three-stage TS1 agreement checks at 100 and 10,000 cells, including input,
  reference, output and executable hashes.
- `tier3-clang-cases.txt` and `tier3-gcc-cases.txt` list every supplied case ID.
  These full build/test logs also record the scenario-contract checks.
  `methods-clang.txt` and `methods-gcc.txt` record the new-method analytical,
  order, failure and projection checks. `prototype-clang.csv`
  and `prototype-gcc.csv` record 147 independent comparisons per compiler.
- `stress-references.json` records the independently generated 48-hour cases,
  two tight reference integrations and their bidirectional convergence scores.
- `input-manifest.json` records hashes for 80 input, reference, generator and
  scoring files. It confirms that all six 10,000-cell inputs preserve the
  disclosed 1,000-cell prefix, photolysis and scenario settings exactly.
- `operating-summary.json` is generated from the final scan CSVs by
  `tools/report_results.py`; the same script generates the top-level
  [RESULTS.md](../RESULTS.md).
- `copy-elision.json` verifies unchanged outputs and operation counts after
  removing redundant Jacobian storage. `sdirk-projection-regression.txt` records
  the projected-error defect and regression. `sdirk-fix-operating-points.json`
  verifies all 192 operating-point outputs and counts against earlier runs,
  and reproducibility across the five final timing repeats.

Every native research run, including failures and superseded measurements,
is retained in [scans/](../scans/). The final selected accuracy runs use all
1,000 disclosed cells and all 96 cells of each stress case. The 10,000-cell
timing runs have reference accuracy only for their disclosed 1,000-cell prefix.
The final requirement and provenance audit is in [COMPLETION.md](COMPLETION.md).

## Historical Tier 1 evidence

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
in [TIER1_RESULTS.md](TIER1_RESULTS.md) and the top-level `README.md`.

Full CSVs, build products, intermediate sweeps, and raw build logs remain under
the ignored `runs`, `build`, `bin`, and `lib` directories. They are not required
to rebuild or run the C port.
