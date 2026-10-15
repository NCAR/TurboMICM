# Portable chemical kinetics in C

Implementation in progress. See DESIGN.md and RESULTS.md for current coverage.

The read-only reference inputs are expected at `../micm` and `../reference`.
`make env` creates the local conda tooling environment; `make cpu` builds the C
library and TS1 benchmark, `make test-cpu` runs the case suite, and `make bench`
builds the benchmark. `make micm-ref` is the only target that builds MICM code.
