# C CPU chemistry port

This repository is the ISO C11 CPU port of MICM for the task described in
[PROMPT.md](PROMPT.md). The library is compiled as C and has no runtime Python
dependency. Python reads JSON/YAML configurations during the build and emits C
mechanism descriptors.

## Build and run

Create the local tooling environment with `make env`, then build with either
`make cpu CC=clang` or `make cpu CC=gcc-15`. The default optimization setting is
`OPTFLAGS=-O3 -DNDEBUG`. The benchmark command is:

```text
bin/ts1_bench_cpu [--threads N] CELLS STEPS DT CELLS_CSV
```

It writes `ts1_final.csv` in the current directory and reports elapsed time and
solver statistics. `make test-cpu` reports the MICM case list in the required
one-line format. `make micm-ref` builds the reference executable when needed.

See [DESIGN.md](DESIGN.md) for the architecture and [DECISIONS.md](DECISIONS.md)
for format and porting decisions. Results are tracked in [RESULTS.md](RESULTS.md).
