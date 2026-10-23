# Exercise C: luna-xhigh's generator on other mechanisms

Entry: `tier2/luna/xhigh`, commit 1ad13ab, in a detached worktree at `$WORK/luna-xhigh`
(`WORK=~/NSF-NCAR/tier2-work/exercises/c/luna-xhigh`). Compiler: Apple clang, `OPTFLAGS="-O3 -DNDEBUG"`, 1 thread.
2026-10-23.

**Route: the benchmark programs.** No change to the entry was needed. All three mechanisms
convert, generate, compile with no warnings, run, and agree with each other and with MICM
to 1e-10 relative, within 6.3e-16 to 5.0e-15.

## Results

| Mechanism | Species / reactions | Generator | Generated C | `make bench` | Specialized vs generic | Generic vs MICM | Specialized vs MICM |
|---|---|---:|---:|---:|---|---|---|
| chapman | 4 (+M) / 7 | 0.04 s | 602 lines, 27,391 B | exit 0, 0.41 s | PASS 0 (bit-identical) | PASS 6.28e-16 | PASS 6.28e-16 |
| carbon_bond_5 | 66 (+M) / 200 | 0.132 s | 16,456 lines, 815,192 B | exit 0, 7.13 s | PASS 1.51e-15 | PASS 4.18e-15 | PASS 4.33e-15 |
| carbon_bond_5_no_emissions | 66 (+M) / 186 | 0.130 s | 16,190 lines, 801,680 B | exit 0, 6.97 s | PASS 3.9e-15 | PASS 4.81e-15 | PASS 4.98e-15 |
| *TS1, for scale* | 209 (+M) / 547 | 0.56 s | 48,068 lines, 2,404,083 B | exit 0, 66.4 s | | | |

- **Generator:** `tools/generate.py` alone, median of 5 runs.
- **`make bench`:** run after deleting `build/generated/ts1.c`, with the library and driver
  objects already built. It covers generation, compiling the generated C, and linking three
  programs.
- **Clean builds** (`make cpu bench` from an empty build tree): Chapman 1.92 s, CB5 8.28 s,
  CB5 without emissions 8.01 s, TS1 64.5 s and 68.8 s.
- **Solver statistics:** both programs took the same path on every mechanism: 540 internal
  steps, 0 rejected, over the 60 timed steps on 4 cells.
- **Timed loop:** generic against specialized was 0.243 against 0.128 ms for Chapman, 10.3
  against 7.6 ms for CB5, and 10.1 against 7.7 ms for CB5 without emissions.

## Addendum: 64 distinct cells

- **Binaries:** the same ones as the main runs (`$WORK/bin/<mech>/`, built by
  `make bench ... CONFIG=...`).
- **Runs:** `--threads 1 64 60 60 ~/NSF-NCAR/tier2-work/exercises/canonical/<mech>/cells64.csv`,
  each in its own empty directory, `$WORK/runs64/<mech>/{generic,specialized}`.
- **Reference:** `micm-run64/ts1_final.csv`, which is identical to `micm-run64-g32/`.
- **Logs:** `logs/<mech>-64-*`.

| Mechanism | Specialized vs generic | Generic vs MICM | Specialized vs MICM | Generic / specialized timed loop | Generated kernels used |
|---|---|---|---|---:|---|
| chapman | PASS 0 (bit-identical) | PASS 2.46e-15 | PASS 2.46e-15 | 1.26 / 0.84 ms | yes |
| carbon_bond_5 | PASS 2.7e-14 | PASS 1.55e-14 | PASS 1.57e-14 | 60.0 / 46.2 ms | yes |
| carbon_bond_5_no_emissions | PASS 2.24e-14 | PASS 2.31e-14 | PASS 3.44e-14 | 60.1 / 45.4 ms | yes |

On every mechanism, both programs exit 0 with 540 internal steps and none rejected.

**Did the specialized program run the generated kernels?** Yes, in the 4-cell runs and the
64-cell runs alike. This entry has no fallback.

- **Driver:** `src/ts1_bench.c`, `integrate_group()` (lines 262–276). The Makefile compiles
  it into `build/ts1_specialized.o` with `-DCHEM_USE_SPECIALIZED`, and in that build it
  unconditionally returns `chem_generated_integrate_group_steps(cells, ...)`. There is no
  cell-count or group-size test.
- **Grouping:** `main()` makes `min(threads, cells)` groups, so `--threads 1` gives one group
  of 4 or 64 cells.
- **Generated function:** it loops over any `cell_count`, checking only validity and
  overflow. It calls only its own static functions: `generated_rates`, `generated_forcing`,
  `generated_jacobian`, `generated_factor`, `generated_solve`, and `generated_error`.
- **Linkage:** `nm bin/<mech>/ts1_bench_specialized` shows no `_chem_integrate_group` or
  `_chem_integrate_group_steps` (the generic integrators), while `ts1_bench_cpu` has both.
  The specialized program cannot reach the generic solver.
- **Consistent with this:** the specialized timed loop is faster than the generic one at 4
  and 64 cells. CB5's two final states differ in the last digits, as two operation orders
  would; Chapman's are bit-identical.

## How the entry represents and builds a mechanism

- **Generator input:** `tools/generate.py CONFIG OUT.c` reads one of two formats.
  - **The entry's own descriptor JSON:** top-level `species` and `reactions`, the schema of
    `config/ts1.json`, which `tools/translate_ts1.py` writes from MICM's TS1 builder header.
  - **MICM/CAMP-format YAML or JSON** (`camp-files`/`camp-data`): used for the tutorial
    mechanisms in `config/rate_constants_*`, which `make test-specialized` and `make test-cpu`
    exercise.
- **Output:** one C file holding both the `chem_mechanism` descriptor that the generic
  Tier 1 library integrates (`chem_generated_mechanism()`) and the specialized three-stage
  Rosenbrock solver (`chem_generated_integrate_group_steps()`).
- **No C reader:** there is no run-time configuration reader in C, so the generic program
  `bin/ts1_bench_cpu` is also rebuilt for each mechanism.
- **Build input:** the Makefile reads `config/ts1.json` through a variable,
  `CONFIG := config/ts1.json`, and always writes `build/generated/ts1.c`.
- **Benchmark driver:** `src/ts1_bench.c` (both programs) is mechanism-generic. It matches
  CSV columns to the descriptor's species and parameter names, refuses unknown columns,
  takes initial values from the descriptor for anything the CSV lacks, and writes
  `ts1_final.csv` with the integrated species, without M.

## The converter (`$WORK/convert.py`, 109 lines, about 60 of code)

It writes the descriptor schema, mirroring `config/ts1.json` key for key. A check confirmed
the same top-level, species, and reaction keys in the same order.

| Canonical | Entry's descriptor (as in TS1) |
|---|---|
| `ARRHENIUS` A B C D E | `ArrheniusRateConstant`, `parameters` {A, B, C, D, E}, `label` null |
| `TROE` k0_A … N | `TroeRateConstant`, `parameters` {k0_A, k0_B, k0_C, kinf_A, kinf_B, kinf_C, Fc, N} |
| `USER_DEFINED` label, scaling_factor (photolysis, emissions) | `UserDefinedRateConstant`, `label`, `parameters` {scaling_factor}; emissions keep `reactants: []` |
| species M, `third_body: true` | species M, `third_body: true`, `initial` 0.0, still named in reactant and product lists. The generator turns each reactant M into a factor of air density and drops product M. |
| reactants with repeats, products [name, yield] | reactant names with repeats, products {species, yield} |
| `initial` temperature, pressure, concentrations, parameters | top-level `temperature` and `pressure`; each species' `initial` (0.0 if absent); top-level `parameters` {label: value} |
| (no molecular weight or diffusion coefficient) | 0.0 and null: used only by surface reactions, which these mechanisms lack |

- **Units:** no conversion. Both forms use MICM's parameters in mol m-3, and the C-sign and
  Troe conventions are the same as MICM's.
- **Why not the CAMP format:** on that path the generator multiplies Troe `k0_A` and
  `kinf_A` by molecule-cm3 factors, because it assumes the tutorials' units. The
  descriptor format takes mol m-3 values as they are, as TS1 does.

## Commands

```sh
WORK=~/NSF-NCAR/tier2-work/exercises/c/luna-xhigh
mkdir -p $WORK && cd $WORK
ln -sn ~/NSF-NCAR/ODE/benchmarks/micm micm
ln -sn ~/NSF-NCAR/ODE/benchmarks/reference reference
ln -sn ~/NSF-NCAR/ODE/benchmarks/tier2/scenarios scenarios
git -C ~/NSF-NCAR/TurboMICM worktree add --detach $WORK/luna-xhigh tier2/luna/xhigh

# Instead of `make env` (see Deviations): the run root's own pinned environment, read-only
export PYTHONDONTWRITEBYTECODE=1
PY=~/NSF-NCAR/ODE/benchmarks/tier2/luna-xhigh/.conda/bin/python
cd $WORK/luna-xhigh && make cpu bench CC=clang OPTFLAGS="-O3 -DNDEBUG" PYTHON=$PY   # exit 0, 64.5 s

# Sanity: in $WORK/sanity/{generic,specialized}; specialized vs generic: PASS 3.46e-14
../../luna-xhigh/bin/ts1_bench_cpu --threads 1 100 100 30 ../../reference/ts1_cells_100_seed1.csv
../../luna-xhigh/bin/ts1_bench_specialized --threads 1 100 100 30 ../../reference/ts1_cells_100_seed1.csv

# Convert
cd $WORK
for m in chapman carbon_bond_5 carbon_bond_5_no_emissions; do
  python3 convert.py ../../canonical/$m/canonical.json converted/$m.json
done

# Per mechanism, as measure.py does it (logs in $WORK/logs, results in measure.json)
cd $WORK/luna-xhigh
$PY tools/generate.py $WORK/converted/$m.json $WORK/gen-timing/$m.c        # 5 times, median
rm -f build/generated/ts1.c
make bench CC=clang OPTFLAGS="-O3 -DNDEBUG" PYTHON=$PY CONFIG=$WORK/converted/$m.json
cp bin/ts1_bench_cpu bin/ts1_bench_specialized $WORK/bin/$m/; cp build/generated/ts1.c $WORK/generated/$m.c
(cd $WORK/runs/$m/generic     && $WORK/bin/$m/ts1_bench_cpu         --threads 1 4 60 60 ~/NSF-NCAR/tier2-work/exercises/canonical/$m/cells.csv)
(cd $WORK/runs/$m/specialized && $WORK/bin/$m/ts1_bench_specialized --threads 1 4 60 60 ~/NSF-NCAR/tier2-work/exercises/canonical/$m/cells.csv)
C=~/NSF-NCAR/ODE/benchmarks/tier2/reference/compare_final.py; R=~/NSF-NCAR/tier2-work/exercises/canonical/$m/micm-run/ts1_final.csv
python3 $C $WORK/runs/$m/specialized/ts1_final.csv $WORK/runs/$m/generic/ts1_final.csv
python3 $C $WORK/runs/$m/generic/ts1_final.csv $R
python3 $C $WORK/runs/$m/specialized/ts1_final.csv $R

# Addendum: 64 distinct cells, with the same saved binaries
(cd $WORK/runs64/$m/generic     && $WORK/bin/$m/ts1_bench_cpu         --threads 1 64 60 60 ~/NSF-NCAR/tier2-work/exercises/canonical/$m/cells64.csv)
(cd $WORK/runs64/$m/specialized && $WORK/bin/$m/ts1_bench_specialized --threads 1 64 60 60 ~/NSF-NCAR/tier2-work/exercises/canonical/$m/cells64.csv)
#   compared as above, against canonical/$m/micm-run64/ts1_final.csv
nm $WORK/bin/$m/ts1_bench_specialized | grep -E ' T _chem_(integrate_group|integrate_group_steps)$'   # nothing

# Clean-build times, then restore TS1 (a clean build without CONFIG)
rm -rf build bin; make cpu bench CC=clang OPTFLAGS="-O3 -DNDEBUG" PYTHON=$PY CONFIG=$WORK/converted/$m.json
rm -rf build bin; make cpu bench CC=clang OPTFLAGS="-O3 -DNDEBUG" PYTHON=$PY
cmp build/generated/ts1.c $WORK/gen-timing/ts1.c        # identical to the original TS1 C
# The sanity runs repeated in $WORK/sanity-after: ts1_final.csv byte-identical to before
git -C $WORK/luna-xhigh status --porcelain               # empty
```

## Deviations and notes

- **`make env` not run.** In this entry it runs `conda env create --prefix .conda` with
  `CONDA_PKGS_DIRS` set to a new `.conda_pkgs` in the worktree. That downloads Python and
  PyYAML from conda-forge, which the rules forbid (no network access, no installs).
  - **What I used instead:** the Makefile's own knob (`PYTHON ?= ./.conda/bin/python`),
    pointed at the environment the agent built during its run,
    `ODE/benchmarks/tier2/luna-xhigh/.conda`. It has Python 3.12.8 and PyYAML 6.0.2,
    exactly the pins in `environment.yml`.
  - **No writes:** I ran it with `PYTHONDONTWRITEBYTECODE=1`, and a `find -newer` showed
    nothing written under ODE.
  - **Same output:** the TS1 C it generated matched the build's byte for byte.
- **No swap of `config/ts1.json`.** The build takes the configuration as the make variable
  `CONFIG`, so each converted mechanism was passed as `CONFIG=...`. The tracked file was
  never touched, so it needed no backup or restore.
  - **Forcing regeneration:** make decides by mtime and always writes
    `build/generated/ts1.c`, so I deleted that file before each mechanism's build.
- **Surprises.**
  - **Nothing failed.** Emissions (reactions with no reactants) are handled: each one only
    adds to its product's forcing, and the entry's own tutorial tests include a CAMP
    emission. A reactant M becomes `value *= air_density` in the rate code.
  - **Fast generator, large output:** CB5 generates in 0.13 s but makes 0.8 MB of C,
    because the generator emits a separate `#pragma omp simd` loop for every LU update
    inside per-pivot functions. Its `make bench` takes about 7 s; TS1's (2.4 MB) takes 66.4 s.
  - **TS1 names stay:** the generated file keeps the name `ts1.c`, and the outputs keep the
    names `ts1_final.csv` and the "TS1 CSV" messages, whatever the mechanism. This is
    cosmetic only.
- **Files** in `$WORK`:
  - `convert.py`, the converter;
  - `converted/*.json`, its output;
  - `measure.py`, orchestration only, calling no entry API;
  - `measure.json`, the raw measurements;
  - `generated/*.c`, `bin/<mech>/`, `runs/<mech>/{generic,specialized}/ts1_final.csv`, and
    `logs/`.
