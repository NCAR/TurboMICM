# Results

## Build and configuration checks

- `make cpu CC=clang`: passed on the provided macOS machine.
- `make cpu CC=gcc-15`: passed on the provided macOS machine.
- Configuration reader: JSON and YAML tutorial examples parse to 7 species and
  9 and 6 reactions, respectively.
- TS1 importer: parsed the expected 210 species and 547 reactions.
- C code generator: emitted and compiled tables for all 210 species and 547
  reactions under both compilers; `bin/ts1_bench_cpu --inspect` reports those
  counts.
- `make test-cpu CC=clang`: emitted all 680 case IDs; the two Arrhenius
  calculation cases pass and the other 678 currently fail as unported.
- `make test-cpu CC=gcc-15`: emitted the same 680 statuses. The C core checks
  passed before the harness returned failure for unported cases.

## MICM cases and benchmark

Only the two Arrhenius calculation cases have a current C port. The remaining
678 cases fail with an explicit incomplete-port reason. The TS1 integrator, CSV
driver, final-state comparison, and timing benchmark are not implemented, so
there are no MICM comparison or performance results.
