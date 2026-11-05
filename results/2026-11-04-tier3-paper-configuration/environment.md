# Environment of the 4 November 2026 session

- Machine: Apple M5 Pro laptop, 18 cores (6 performance, 12 efficiency), 24 GB,
  macOS 26.7.2; the second machine of the Tier 2 verifications, not the
  10-core M5 laptop of the Tier 3 scans and ranking. Timings here compare only
  with each other.
- Compiler: Apple clang 21.0.0 (clang-2100.3.34.2), `-O3 -DNDEBUG`, OpenMP from
  Homebrew's libomp under the home directory (`OMP_CFLAGS` and `OMP_LIBS`
  passed to `make bench`).
- Library: branch `tier3/fable/paper` at `296beb4` (the reference solution
  `tier3/fable/max` at `cb0f8a3` plus the paper's three interval rules as
  options), built in the worktree `~/NSF-NCAR/TurboMICM-paper`.
- MICM baseline: `micm_scenario` built from `tier3/scenarios/micm_scenario.cpp`
  against the MICM sources kept in the ODE repository's `benchmarks/micm`
  (the paper's MICM at commit e15898ce), threads from `std::thread`.
- Scenarios: the six disclosed 10,000-cell scenarios from
  `~/NSF-NCAR/tier2-work/check/disclosed-10000`; the five held-out ones
  regenerated at 10,000 cells with `tier3/scenarios/make_scenario.py` and the
  seeds and arguments in their headers (`--defaults` from the ODE repository's
  `benchmarks/reference/ts1_defaults.csv`); the first 64 cells and the
  photolysis tables of all eleven are byte-identical to the public archive's.
- References for scoring: the archive's 64-cell references
  (`Fast-ODE-Atmospheric-Chem/data/scenarios/*/reference_*.csv`); the scorer
  therefore scores the first 64 cells of each 10,000-cell run and applies the
  gate to every cell.
- Protocol: `admin/time_entries.py` through `time_paper.sh`: three rotating
  rounds at 10 threads, 5 s pauses, one scored run per program and scenario
  unless outputs differ.
