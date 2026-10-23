# Exercise C: luna-max on other mechanisms

Entry: `tier2/luna/max` (3562d5a), worktree `luna-max/`. Route: **the benchmark programs**, built
by the entry's own Makefile with the converted mechanism passed as `TS1_JSON=...`. No driver and
no change to the entry. The tracked `mechanisms/ts1.json` was never swapped; its SHA-256 matches
HEAD and `git status --porcelain` is empty.

## Results

4 cells, one 60 s warm-up step, then 60 steps of 60 s; MICM's default three-stage Rosenbrock
(relative 1e-6, absolute 1e-3); clang 21, `-O3 -DNDEBUG`, 1 thread. Each comparison uses
`compare_final.py` (relative 1e-10). The number after PASS is the largest relative difference.

| Mechanism | Converted | Generator (s) | Specialized C: lines / bytes | `make bench` exit / s | Specialized vs generic | Generic vs MICM | Specialized vs MICM |
|---|---|---:|---:|---:|---|---|---|
| chapman | yes | 0.025 | 304 / 12,126 | 0 / 0.89 | PASS 0 | PASS 6.28e-16 | PASS 6.28e-16 |
| carbon_bond_5 | yes | 0.034 | 9,682 / 478,951 | 0 / 1.38 | PASS 0 | PASS 3.46e-15 | PASS 3.46e-15 |
| carbon_bond_5_no_emissions | yes | 0.034 | 9,542 / 471,952 | 0 / 1.41 | PASS 0 | PASS 4.46e-15 | PASS 4.46e-15 |
| *TS1, for comparison* | (the entry's own) | 0.148 | 28,434 / 1,453,684 | 0 / 3.11 | | | |

- **Generator:** wall time of `python3 tools/generate_c.py IN OUT --specialized`, median of 5
  runs. The generic-table mode takes the same time (0.025, 0.033, 0.032 s). `make bench` runs
  both modes. The generic-table files are 41 / 3,395, 332 / 73,748, and 304 / 69,407 lines /
  bytes.
- **`make bench`:** full build after `make clean`, median of 3 runs (chapman 0.93, 0.89,
  0.88 s; CB5 1.37, 1.39, 1.38 s; CB5 without emissions 1.41, 1.39, 1.41 s). It includes both
  generator runs, compiling the core sources, three links, and `ts1_scenario_specialized`.
  Compiling the specialized file alone takes 0.058, 0.527, and 0.514 s. The builds printed no
  warnings under `-Wall -Wextra -Wpedantic`.
- **Solver statistics:** both programs took 540 internal steps (9 per interval) with none
  rejected. The single-run one-thread elapsed times for the 60 timed steps, generic then
  specialized:
  - chapman: 0.34 and 0.21 ms;
  - CB5: 13.7 and 3.7 ms;
  - CB5 without emissions: 13.4 and 3.7 ms.

  These are single runs, not a benchmark.
- **Emissions:** the 14 emissions in carbon_bond_5 take effect. With them, the final NO is
  6.62e-10 mol m-3; without them it is 2.94e-13, and MICM gives the same values. 61 of 66
  species differ between the two runs.

## Addendum: 64 distinct cells

These runs use the same route and default settings as above. Each mechanism was rebuilt with
`make clean && make bench CC=clang OPTFLAGS="-O3 -DNDEBUG" TS1_JSON=...`, and the binaries are
byte-identical to the main runs'. Each program ran with
`bin/ts1_bench_{cpu,specialized} --threads 1 64 60 60 $C/$M/cells64.csv` in its own empty
directory (`runs64/$M/{generic,specialized}`), and the results were compared with
`$C/$M/micm-run64/ts1_final.csv`. `micm-run64-g32/` holds the identical file and gives the same
results.

| Mechanism | Specialized vs generic | Generic vs MICM | Specialized vs MICM | Generated kernels used (64 cells / 4 cells) | Elapsed, generic / specialized |
|---|---|---|---|---|---:|
| chapman | PASS 0 | PASS 2.46e-15 | PASS 2.46e-15 | true / true | 5.13 / 3.07 ms |
| carbon_bond_5 | PASS 0 | PASS 2.97e-14 | PASS 2.97e-14 | true / true | 219.9 / 62.4 ms |
| carbon_bond_5_no_emissions | PASS 0 | PASS 2.28e-14 | PASS 2.28e-14 | true / true | 219.2 / 60.5 ms |

Every program exited 0 with 540 internal steps and none rejected. One thread with the default
grouping puts all 64 cells in one group.

**How the kernels are chosen (`src/solver.c`).** None of the selection tests involves the cell
count or the group size:
- **Rates:** `chem_solver_create` (lines 862-863) computes each cell's rate constants with
  `kernels->rates` whenever the model has a kernel table.
- **Evaluate:** `evaluate_rhs` (255-258) calls `kernels->evaluate` when three conditions hold:
  there are no constraints, `kernels->species_count == pattern.n`, and
  `kernels->pattern_nonzeros == pattern.nnz`.
- **Factor and solve:** `solve_group` (673-682) installs `kernels->factor` and
  `kernels->solve` under the same size test, provided the factor and solve callbacks are the
  defaults. `chem_rosenbrock_default` sets them to `default_factor_callback` and
  `default_solve_callback`, and `benchmark.c` never changes them.

A cell group only splits cells across threads and shares one step size. Every cell goes through
the same per-cell calls (the loop at line 700), so this entry has no full-group or vector-width
path to fall back from.

The specialized programs link `.kernels=&chem_generated_kernels`. Its `pattern_nonzeros` is 14
for chapman and 743 for both CB5 variants, and the binaries keep the `generated_*` symbols.

**Dynamic check.** I tried two ways to observe the benchmark binary directly, and both hung
waiting for debugger authorization (developer mode is off), so I killed them:
- `lldb --batch` with auto-continue breakpoints on `generated_factor` and
  `default_factor_callback`;
- `sample` on a long run.

Instead, a diagnostic driver checked the selection. `kernel_check/kernel_check.c` is 132 lines
and uses only the public `kinetics.h` API. It links the same solver sources and generated model,
and copies `chem_model` with the generated kernel table. The table keeps its sizes, but each
function is wrapped by a counter that forwards to the generated one. The driver runs the same
problem: same CSV and cell count, group size 0, one thread, a warm-up step, then 60 steps of
60 s.

| Cells | Rates | Evaluate | Factor | Solve | Solver totals incl. warm-up |
|---:|---:|---:|---:|---:|---|
| 4 | 4 | 4,392 | 2,196 | 6,588 | function_calls 1,098, decompositions 549, solves 1,647 |
| 64 | 64 | 70,272 | 35,136 | 105,408 | the same |

The counts are identical for all three mechanisms. They equal the solver's totals times the
number of cells, so every evaluate, factor, and solve call went through the generated kernels;
had the solver fallen back, the counts would be zero. The driver's final states are
byte-identical to `ts1_bench_specialized`'s in all six runs.

The output files cannot show this by themselves: the generated kernels reproduce the generic
arithmetic exactly (specialized vs generic: 0). The elapsed times are consistent with it, though:
the specialized program is 3.5 times faster on CB5 at 64 cells and 3.7 times at 4 cells.

## Commands

```sh
WORK=~/NSF-NCAR/tier2-work/exercises/c/luna-max; C=~/NSF-NCAR/tier2-work/exercises/canonical
# setup as given (links, worktree), then the TS1 sanity build and checks
cd $WORK/luna-max && make cpu bench CC=clang OPTFLAGS="-O3 -DNDEBUG"
#   bin/ts1_bench_{cpu,specialized} --threads 1 100 100 30 ../reference/ts1_cells_100_seed1.csv
#   both PASS against ../reference/micm_ts1_final_100_seed1.csv (3.14e-14) and agree exactly
# convert
python3 $WORK/convert.py $C/$M/canonical.json $WORK/mechanisms/$M.json
# build (timed), per mechanism
make clean && make bench CC=clang OPTFLAGS="-O3 -DNDEBUG" TS1_JSON=$WORK/mechanisms/$M.json
# run, each in its own empty directory ($WORK/runs/$M/{generic,specialized})
bin/ts1_bench_cpu         --threads 1 4 60 60 $C/$M/cells.csv
bin/ts1_bench_specialized --threads 1 4 60 60 $C/$M/cells.csv
# compare
python3 ~/NSF-NCAR/ODE/benchmarks/tier2/reference/compare_final.py specialized/ts1_final.csv generic/ts1_final.csv
python3 .../compare_final.py generic/ts1_final.csv     $C/$M/micm-run/ts1_final.csv
python3 .../compare_final.py specialized/ts1_final.csv $C/$M/micm-run/ts1_final.csv
# finally: make clean && make bench (TS1 again); git status --porcelain is empty
```

Two further checks:
- The generated C files that the builds wrote are byte-identical to the standalone generator's
  output.
- Rerunning each mechanism's programs from a second build reproduced the final states exactly.

## What it took

**How the entry represents a mechanism.** `tools/generate_c.py` reads either a normalized JSON
mechanism or a MICM YAML/JSON configuration with `camp-files`. `tools/translate_ts1.py` writes
`mechanisms/ts1.json` from MICM's TS1 header. That file is the configuration the build reads,
through the Makefile variable `TS1_JSON`.

The generator emits two files from it:
- `build/ts1_model.c`: tables for the generic solver;
- `build/ts1_model_specialized.c` (`--specialized`): literal rate formulas, an unrolled
  reaction and Jacobian update, a no-pivot sparse LU with diagonal Markowitz ordering, and the
  triangular solves.

`src/benchmark.c` is generic over the linked `chem_model`. Its columns are matched by name, and
any CSV column that is not a species or rate parameter is an error. The README documents
`make mechanism CONFIG=...` for other mechanisms, but that target builds only a generic
library. `TS1_JSON` is set with `:=` and is not documented as an override, but a command-line
assignment overrides it in GNU make.

**What the converter maps** (`convert.py`, 124 lines). The canonical format is already close to
the entry's:
- **Rate parameters:** nested under `parameters`: {A, B, C, D, E} for ARRHENIUS; {k0_A, k0_B,
  k0_C, kinf_A, kinf_B, kinf_C, Fc, N} for TROE; {scaling_factor} plus a top-level `label` for
  USER_DEFINED.
- **Products:** `[name, yield]` pairs become `{name, yield}` objects.
- **Third body:** M stays a species with `third_body: true` and `initial: 0.0`. It also stays
  by name in reactant and product lists, as in TS1's `O2 + M + O -> O3 + M`. The generator turns
  each M reactant into a factor of air density and drops M products, as MICM does.
- **Emissions:** they keep an empty reactant list. The generator needs no special case: flux is
  k, with no Jacobian entries.
- **Initial values:** top-level `temperature` and `pressure`, each species' `initial` (0.0 when
  unset), and `custom_parameters` for every rate parameter. The five CB5 photolysis rates
  without an initial value are written as 0.0, as cells.csv and MICM's InitState set them.
- **Species properties:** `molecular_weight` 0.0 and `diffusion_coefficient` null, as ts1.json
  writes for species without them. Only SURFACE reactions read them, and these mechanisms have
  none.

A check confirmed that every converted reaction has the same key structure as ts1.json's
reactions of its type, and that every number survives the conversion exactly.

**Surprises.**
- `make env` did not use an existing environment. It creates a new conda environment (python
  3.12.11, pyyaml 6.0.2) from conda-forge, with a fresh package cache inside the worktree. That
  needs the network, so it was run only with `CONDA_OFFLINE=true`, where it failed:
  `PackagesNotFoundError: ... python=3.12.11`. The empty cache directory it left was removed.
  Without `.conda/bin/python`, the Makefile falls back to `python3` (3.13.1 here). That is
  enough: the generator imports PyYAML only for YAML input.
- Otherwise nothing went wrong. The generator accepted every converted file as is, and the
  builds produced no warnings.
- The specialized path agrees with the generic one bit for bit, as the entry's own
  `test-specialized` claims for TS1.
- Both paths match MICM to 4.5e-15 or better on all three mechanisms.

Files: `convert.py`, `mechanisms/*.json` (converted), `gen/` (standalone generator output and
`generator_timing.json`), `bins/<mech>/` (programs and generated C from each build), `runs/<mech>/`
(final states and logs), and `logs/` (build logs), all in `$WORK`.
