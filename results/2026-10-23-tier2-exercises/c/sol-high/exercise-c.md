# Exercise C, sol-high: the generator on another mechanism

Entry: `tier2/sol/high` (07848a8), worktree `$WORK/sol-high`, where
`WORK=~/NSF-NCAR/tier2-work/exercises/c/sol-high`. Route: **the benchmark programs**
(swap the TS1 configuration and rebuild). No change to the entry was needed.
At the end, `git status --porcelain` is empty and `configs/ts1.json` matches its
backup (SHA-256 `e21000a7…f5b7`).

## Results

Protocol: the default `make bench CC=clang OPTFLAGS="-O3 -DNDEBUG"`, then
`bin/ts1_bench_cpu` and `bin/ts1_bench_specialized` with `--threads 1 4 60 60 cells.csv`,
each in its own empty directory. Comparisons use `compare_final.py` (rtol 1e-10) and
report the largest relative difference.

| | chapman | carbon_bond_5 | carbon_bond_5_no_emissions |
|---|---|---|---|
| converted / generator ran | yes / yes | yes / yes | yes / yes |
| generator, CLI (in-process) | 0.032 s (0.5 ms) | 0.041 s (8.5 ms) | 0.040 s (8.4 ms) |
| generated C (width 32) | 452 lines, 17,090 B | 9,685 lines, 352,001 B | 9,419 lines, 342,035 B |
| generated C compile alone | 0.23 s | 2.80 s | 2.48 s |
| `make bench` after swap (clean) | 3.63 s (4.38 s), exit 0 | 6.20 s (6.94 s), exit 0 | 5.72 s (6.82 s), exit 0 |
| generic / specialized ran | yes / yes | yes / yes | yes / yes |
| specialized vs generic | PASS 0 † | PASS 0 † | PASS 0 † |
| generic vs MICM | PASS 6.28e-16 | PASS 4.58e-15 | PASS 4.12e-15 |
| specialized vs MICM | PASS 6.28e-16 † | PASS 4.58e-15 † | PASS 4.12e-15 † |
| **Supplementary: `SPECIALIZE_WIDTH=4`** | | | |
| build | 3.38 s, exit 0 | 5.48 s, exit 0 | 5.37 s, exit 0 |
| specialized vs generic | PASS 0 | PASS 4.38e-15 | PASS 3.61e-15 |
| specialized vs MICM | PASS 6.28e-16 | PASS 1.77e-15 | PASS 5.66e-15 |
| elapsed ms: generic / default spec. / width-4 spec. | 0.228 / 0.221 / 0.109 | 10.68 / 10.62 / 3.48 | 10.69 / 10.19 / 3.32 |

For reference, TS1 in the same session: generator 0.074 s CLI (33 ms in-process),
26,300 lines and 990,603 bytes of generated C, which compiles alone in 13.6 s; the
TS1 rebuild after the restore took 16.7 s.

† **Caveat: the default build's specialized program does not run the generated
kernels on 4 cells.** `make bench` generates the kernels for a fixed width of 32
cells (`--fixed-width $(SPECIALIZE_WIDTH)`, which defaults to `GROUP_SIZE` = 32), and
each generated kernel starts with `if (s->cells != 32) { chem_forcing(s, y, f); return; }`
(likewise for the Jacobian, alpha shift, factorization, solve, and rates). The
4-cell problem forms one 4-cell group, so the specialized program runs the generic
chemistry and linear-algebra kernels. Only the generated Rosenbrock loop is its own
(fused stages and a SIMD error norm, both of which keep the generic order of
arithmetic). Its output is bit-identical to the generic program's, as expected. The
supplementary rows rebuild with `SPECIALIZE_WIDTH=4`, the Makefile's documented
control (DESIGN.md: "positive values specialize that width"). That build does run
the generated kernels (`cells != 4` guards; CB5 runs about 3 times faster than
generic), and all of them pass. Chapman's
width-4 result is still bit-identical to generic: its operation order is unchanged,
and the generator says it keeps accumulation order. For CB5, the reciprocal pivots
change the last bits. Elapsed times are single 4-cell runs, a rough guide only. The
addendum below runs 64 distinct cells with the default build, so the generated
kernels run as shipped.

## Addendum: 64 distinct cells with the default build

The build and route are the same as the main runs: the Makefile's defaults (width 32)
and a swap of `configs/ts1.json`. Each program ran with
`--threads 1 64 60 60 canonical/MECH/cells64.csv` in its own empty directory and was
compared with `micm-run64/ts1_final.csv`. The comparisons against
`micm-run64-g32/` gave the same results; that file is byte-identical.

| | chapman | carbon_bond_5 | carbon_bond_5_no_emissions |
|---|---|---|---|
| `make bench` after swap | 3.41 s, exit 0 | 6.29 s, exit 0 | 5.89 s, exit 0 |
| specialized vs generic | PASS 0 | PASS 2.56e-14 | PASS 2.18e-14 |
| generic vs MICM | PASS 2.46e-15 | PASS 2.29e-14 | PASS 2.28e-14 |
| specialized vs MICM | PASS 2.46e-15 | PASS 1.82e-14 | PASS 2.69e-14 |
| generated kernels used, 64 cells | **yes** (0 fallback calls) | **yes** (0) | **yes** (0) |
| generated kernels used, 4 cells (main runs) | **no** (all calls fall back) | **no** | **no** |
| elapsed ms, generic / specialized (64 cells) | 1.21 / 0.86 | 56.5 / 42.5 | 56.0 / 42.4 |

**The code path that chooses the kernels.** `src/bench.c:26` starts the group size at
`CHEM_GROUP_SIZE`, which is 32 (`GROUP_SIZE`). Lines 42-46 cap it at cells per thread
unless `--group-size` is given, and line 47 caps it at the cell count. With
`--threads 1`, 4 cells form one 4-cell group, and 64 cells form two 32-cell groups.
The default build generates with `--fixed-width 32`, and each of the six generated
kernels branches on the group width:

- `sp_forcing`: `if (s->cells != 32) { chem_forcing(s, y, f); return; }`;
- `sp_jacobian` falls back to `chem_negative_jacobian`;
- `sp_alpha` to `chem_alpha_minus_jacobian`;
- `sp_factor` to `chem_factor`;
- `sp_backsolve` to `chem_linear_solve`;
- `ts1_update_rates` to `chem_update_rates`.

Only the generated Rosenbrock loop (fused stages, SIMD error norm) runs in both
cases.

**Runtime check.** I did not use a debugger, because macOS developer mode is off and
lldb would raise an authorization dialog. Instead, `probe/` compiles the unmodified
`build/ts1_specialized.c` of each mechanism with a forced-include header
(`probe/probe.h`, 6 lines of `#define`). The header renames those six fallback calls
to counting wrappers (`probe/probe.c`, 17 lines), and the result links like
`bin/ts1_bench_specialized`. In every run the probe's `ts1_final.csv` is
byte-identical to the real program's. Counts, the same for all three mechanisms:

- at 4 cells, every call falls back: forcing 1,098, Jacobian 549, alpha 549,
  factor 549, solve 1,647, rate update 1 (549 = 540 timed + 9 warm-up internal steps);
- at 64 cells, no call falls back: all six counts are 0.

The numbers agree. In the 4-cell default runs, specialized output is bit-identical
to generic. At 64 cells, CB5's specialized output differs from generic in the last
bits (the reciprocal-pivot LU) and takes 25 % less time. Chapman stays bit-identical
because its generated operation order equals the generic order. At 32-cell width
the generic kernels vectorize well, so the gain is 1.33 times, against 3.1 times for
the width-4 build at 4 cells. All timings are single runs.

## What it took

- **Format.** The build reads `configs/ts1.json`, a path hard-coded in the Makefile
  (`build/ts1.c` and `build/ts1_specialized.c` both depend on it). The file is extended
  CAMP JSON: a `camp-data` list with one `CHEM_SPEC` per species and one `MECHANISM`.
  The entry's `tools/translate_ts1.py` produced it from MICM's header. Both
  `tools/configure.py` (generic descriptor data) and `tools/specialize.py` (generated
  solver) read it with the same reader, `normalize` plus `reorder`.
- **The converter**, `convert.py` (84 lines), mirrors the TS1 file:
  - species as `{"name", "type": "CHEM_SPEC", "initial concentration"}`, with 0 for
    species absent from `initial.concentrations`;
  - M as `{"initial concentration": 0, "third body": true}`, kept in reactant and
    product lists as TS1 keeps it; the reader drops it from the dynamic species and
    multiplies the rate by air density once per occurrence;
  - reactants as name lists with repeats; products as `{"species", "yield"}` lists;
  - ARRHENIUS and TROE parameters under their own keys;
  - USER_DEFINED with `scaling_factor` and `"MUSICA name"` = label. Emissions are
    the same type with an empty reactant list, and the reader accepts them as it
    accepts its own EMISSION type;
  - `"rate units": "mol m-3"`, `"initial conditions"` (temperature and pressure), and
    `"initial rate parameters"` (label to value).

  Every key set matches TS1's file, except that species have no molecular weight:
  canonical.json has none, and the reader needs it only for SURFACE reactions. The
  entry's reader accepted all three files on the first try (4/66/66 dynamic species,
  7/200/186 reactions, 3/39/25 labels), with no warnings in any build.
- **The drivers are mechanism-generic.** `src/bench.c` and `src/bench_io.h` refer only
  to the symbol `ts1_data`. CSV columns are matched against the descriptor's
  species and label names, and the output columns come from the same descriptor, so
  no driver of my own was needed. Negative yields (PAR in CB5), emissions,
  third-body reactions, and TROE all worked unchanged.
- **The generator** is fast. Most of its CLI time is Python start-up; the in-process
  generation takes at most 8.5 ms for CB5. Of CB5's 5.7 to 6.2 s rebuild, compiling
  the generated C takes 2.5 to 2.8 s. Most of the rest goes to the library sources,
  which the Makefile compiles again into each of the three programs.
- **Setup note:** `make env` downloaded about 25 MB of conda-forge packages (Python 3.12,
  PyYAML) because the Makefile points `CONDA_PKGS_DIRS` at an empty cache inside the
  worktree. That was the only network use.

## Commands

```sh
WORK=~/NSF-NCAR/tier2-work/exercises/c/sol-high
# setup as prescribed: three links, git worktree add --detach $WORK/sol-high tier2/sol/high,
# make env; make cpu bench CC=clang OPTFLAGS="-O3 -DNDEBUG"; TS1 sanity runs in $WORK/sanity
cd $WORK && python3 -I convert.py ../../canonical/MECH/canonical.json converted/MECH/ts1.json
# generator timing (also run 5x under a perf_counter wrapper), output kept in $WORK/gen/MECH
cd $WORK/sol-high && .conda/bin/python tools/specialize.py $WORK/converted/MECH/ts1.json \
    $WORK/gen/MECH/ts1_specialized.c --prefix ts1 --batch 1 --fixed-width 32
$WORK/run_exercise.sh protocol chapman carbon_bond_5 carbon_bond_5_no_emissions
$WORK/run_exercise.sh width4   chapman carbon_bond_5 carbon_bond_5_no_emissions
$WORK/run_exercise.sh clean    chapman carbon_bond_5 carbon_bond_5_no_emissions
$WORK/run_exercise.sh restore
# addendum: default build, 64 distinct cells, the fallback probe; then restore again
$WORK/run_exercise.sh cells64  chapman carbon_bond_5 carbon_bond_5_no_emissions
$WORK/run_exercise.sh restore
```

`cells64` builds with the defaults and checks that `build/ts1_specialized.c` matches
`gen/MECH` (width 32, six `cells != 32` guards). It runs both programs with
`--threads 1 64 60 60 canonical/MECH/cells64.csv` in `runs/MECH/{generic,specialized}-64/`
and compares each pair and both against `micm-run64` and `micm-run64-g32`. It then
builds and runs the probe in `probe/MECH/` at 4 and at 64 cells and checks it against
the real program's output with `cmp`.

`run_exercise.sh` copies the backup of `configs/ts1.json` to `$WORK/backup/` once,
then for each mechanism:

1. it copies `converted/MECH/ts1.json` over `configs/ts1.json`;
2. it runs `/usr/bin/time -p make bench CC=clang OPTFLAGS="-O3 -DNDEBUG" [SPECIALIZE_WIDTH=4]`,
   logged in `logs/`, and checks that `build/ts1_specialized.c` is identical to `gen/MECH`;
3. it runs `bin/ts1_bench_cpu` and `bin/ts1_bench_specialized` with
   `--threads 1 4 60 60 canonical/MECH/cells.csv`, each in a fresh `runs/MECH/<label>/`;
4. it runs `python3 .../tier2/reference/compare_final.py` on each pair.

`restore` copies the backup back, rebuilds, and re-checks TS1 against MICM on the
100-cell reference: generic PASS 4.25e-14, specialized PASS 3.23e-14. It then lists
tracked changes, and there are none. This held after the main runs and again after
the addendum.
