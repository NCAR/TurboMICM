# Exercise C, sol-xhigh: the generator on another mechanism

Entry: `tier2/sol/xhigh` (7af0738), detached worktree at `sol-xhigh/`. Route: **the benchmark
programs**: the converted mechanism replaced `mechanisms/ts1.json`, and `make bench` rebuilt
both `bin/ts1_bench_cpu` (generic) and `bin/ts1_bench_specialized` (generated) for it. No
change to the entry was needed, and no driver was written.

## Results

4 cells, `--threads 1 4 60 60 cells.csv` (one warm-up step of 60 s, then 60 steps of 60 s),
default settings (three-stage Rosenbrock, relative 1e-6, absolute 1e-3), clang
`-O3 -DNDEBUG`. Comparisons use `compare_final.py` (rtol 1e-10) and give the largest relative
difference.

| Mechanism | Converted | Generator (s) | Generated C (lines / bytes) | `make bench` (exit, s) | Generic ran | Specialized ran | Specialized vs generic | Generic vs MICM | Specialized vs MICM |
|---|---|---|---|---|---|---|---|---|---|
| chapman | yes | 0.049 | 670 / 22,471 | 0, 0.65 | yes | yes | PASS 0 | PASS 0 | PASS 0 |
| carbon_bond_5 | yes | 0.075 | 36,218 / 1,105,036 | 0, 4.65 | yes | yes | PASS 0 | PASS 0 | PASS 0 |
| carbon_bond_5_no_emissions | yes | 0.068 | 35,944 / 1,091,562 | 0, 4.84 | yes | yes | PASS 0 | PASS 0 | PASS 0 |
| TS1 (for scale) | n/a | 0.192 | 112,387 / 3,431,254 | 0, 18.14 | yes | yes | PASS 0 | PASS 0 (100-cell check) | PASS 0 (100-cell check) |

Every final value in every comparison is **bitwise equal** to MICM's: 24 of 24 values for
chapman, 272 of 272 for each carbon_bond_5 variant. The entry also writes its output columns
in MICM's order. Both paths take 540 internal steps with 0 rejections on every mechanism. The
emissions take effect: NO ends at 6.62e-10 mol m-3 with them and 2.94e-13 without them.

What was measured:
- **Generator time:** the median of 3 runs of `tools/specialize.py` with the Makefile's
  arguments, Python start-up included. The first runs took 0.151 and 0.161 s, cold. The output
  is byte-identical across runs and to the build's own `generated/sp_ts1.c`.
- **Generated C:** `generated/sp_ts1.c`, the specialized kernels with the descriptors embedded.
  The generic path's descriptor file `generated/ts1.c` (`tools/configure.py`, 0.05 s) is a
  further 34, 661, and 633 lines (4,832, 119,135, and 112,130 B).
- **`make bench`:** the incremental build after the swap. It regenerates both files, compiles
  them, and links the three programs, including the scenario driver; `libchem.a` does not depend
  on the configuration and is not rebuilt. It gave 0 warnings. TS1's row is the rebuild after
  restoring the original configuration.
- **Timed region (not a benchmark: one run each, 4 cells, 1 thread):** generic against
  specialized was 0.226 against 0.100 ms for chapman, 10.20 against 3.41 ms for carbon_bond_5,
  and 11.35 against 3.40 ms for carbon_bond_5_no_emissions.

## Addendum: 64 distinct cells

For each mechanism I rebuilt the entry exactly as in the main runs: the same swap and
`make bench` with the Makefile's defaults. Every mechanism-dependent product was regenerated or
rebuilt, and both generated files were byte-identical to the main run's. Each program then ran
with `--threads 1 64 60 60 canonical/MECH/cells64.csv` in its own empty directory under
`runs64/MECH/`. MICM's file is `micm-run64/ts1_final.csv`; the `micm-run64-g32` file is
identical, and every comparison against it gives the same result.

| Mechanism | Specialized vs generic | Generic vs MICM | Specialized vs MICM | Bitwise equal to MICM | Generated kernels used (64 cells) | Generated kernels used (4 cells) |
|---|---|---|---|---|---|---|
| chapman | PASS 0 | PASS 0 | PASS 0 | 256 / 256 | yes, fixed-stride (`*_fixed`) | yes, runtime-stride (`*_dynamic`) |
| carbon_bond_5 | PASS 0 | PASS 0 | PASS 0 | 4,224 / 4,224 | yes, fixed-stride | yes, runtime-stride |
| carbon_bond_5_no_emissions | PASS 0 | PASS 0 | PASS 0 | 4,224 / 4,224 | yes, fixed-stride | yes, runtime-stride |

Both paths take 1,080 internal steps with 0 rejections, which is 2 groups × 540.

**The evidence: the code path that chooses the kernels.**
- **The driver.** With `-DCK_SPECIALIZED`, the benchmark calls only
  `ts1_specialized.update_rates` and `ts1_specialized.solve` (`src/bench_common.h:7-8`); there
  is no call to the generic `ck_update_rates` or `ck_solve`.
- **The solve.** `ts1_specialized.solve` is the generated `ts1_specialized_solve`, the
  three-stage Rosenbrock template `src/specialized_solver.h` included into `generated/sp_ts1.c`.
  It calls the generated `_forcing`, `_jacobian`, `_factor`, and `_linear` functions, and the
  rate update is the generated `_update`.
- **The choice.** Each of those four generated functions dispatches like this:
  ```c
  if (s->cells == 32)  /* --width 32, from GROUP_SIZE */
    { static const ck_chemistry_kernel kernels[] = {ts1_specialized_forcing_0_fixed, ...}; ck_run_chemistry_kernels(...); }
  else
    { static const ck_chemistry_kernel kernels[] = {ts1_specialized_forcing_0_dynamic, ...}; ck_run_chemistry_kernels(...); }
  ```
  Both helper sets are generated for the mechanism; carbon_bond_5_no_emissions' binary has 118
  of each. The `_fixed` set makes the stride a constant of 32, and the `_dynamic` set takes it
  at run time. `ck_run_*_kernels` (`src/specialized.c`) only loops over the helper array; with
  tile 0, that is the whole group.
- **Group size.** Without `--group-size`, the default width of 32 shrinks to
  ceil(cells / threads) (`src/bench.c:41-52`):
  - 4 cells on 1 thread make one group of 4, which takes the `_dynamic` branch. The summed
    statistics confirm one group: 540 internal steps.
  - 64 cells on 1 thread make 2 groups of 32, which take the `_fixed` branch. The summed
    statistics confirm two groups: 1,080 = 2 × 540.

So the specialized program ran generated kernels in both cases. For this entry the 4 identical
cells do not cause a fallback to generic code. They only select the runtime-stride versions of
the same generated kernels.

**Supporting timing.** These are interleaved runs with the carbon_bond_5_no_emissions build,
the minimum of 3 rounds, 300 steps of 60 s, timed region only. The specialized program is
faster in every configuration, and the fixed-stride branch is the fastest:

| Configuration | Generic (ms) | Specialized (ms) | Ratio | Generated helpers |
|---|---|---|---|---|
| 4 cells | 47.5 | 15.9 | 3.0 | `_dynamic` |
| 64 cells, default 32-cell groups | 255.2 | 183.8 | 1.39 | `_fixed` |
| 64 cells, `--group-size 16` | 301.5 | 227.1 | 1.33 | `_dynamic` |

The single 60-step runs of `run_cells64.py` are too short to time and are not used; for example,
carbon_bond_5 read 54.9 against 53.2 ms. The specialized program's advantage shrinks with larger
groups because the generic library amortizes its indexing over the cells of a group.

**What didn't work:**
- **No runtime profile.** Developer mode is off. I did not use lldb, because it might raise an
  authorization prompt. `/usr/bin/sample` hung, and I killed it after about 2 minutes.
- **A hazard in my harness, now fixed.** macOS's GNU Make 3.81 compares modification times to
  the whole second. My first 64-cell pass copied each configuration in within the second of the
  previous build's outputs, so make kept the previous mechanism's products:
  - carbon_bond_5's specialized program was still chapman's ("Unknown CSV column: AACD");
  - carbon_bond_5_no_emissions' generic program was still carbon_bond_5's, emissions included,
    and failed against MICM.

  `run_cells64.py` now waits for the next second before each copy and checks that every
  mechanism-dependent product is newer than the configuration. The first pass is kept in
  `attempt1-make-timestamps/`. The main 4-cell runs were not affected: their make logs show
  every step regenerating both files, and their generated C matched the standalone generator's.

## Commands

```
WORK=~/NSF-NCAR/tier2-work/exercises/c/sol-xhigh
mkdir -p $WORK && cd $WORK
ln -sn ~/NSF-NCAR/ODE/benchmarks/micm micm
ln -sn ~/NSF-NCAR/ODE/benchmarks/reference reference
ln -sn ~/NSF-NCAR/ODE/benchmarks/tier2/scenarios scenarios
git -C ~/NSF-NCAR/TurboMICM worktree add --detach $WORK/sol-xhigh tier2/sol/xhigh
# .conda/bin/python wrapper (see "What it took"), then:
cd $WORK/sol-xhigh && make env && make cpu bench CC=clang OPTFLAGS="-O3 -DNDEBUG"    # exit 0, 20.5 s
# sanity check: both programs, --threads 1 100 100 30 ../reference/ts1_cells_100_seed1.csv,
# each in its own directory (runs/sanity-ts1/): PASS 0 against MICM

cd $WORK
for m in chapman carbon_bond_5 carbon_bond_5_no_emissions; do
  python3 convert.py ../../canonical/$m/canonical.json converted/$m.json; done
cp -p sol-xhigh/mechanisms/ts1.json backup/ts1.json        # identical to HEAD's (cmp)
python3 run_mech.py chapman
python3 run_mech.py carbon_bond_5 carbon_bond_5_no_emissions
```

For each mechanism, `run_mech.py` (157 lines; it only automates these commands) does the
following:
```
cp converted/MECH.json sol-xhigh/mechanisms/ts1.json
(cd sol-xhigh && .conda/bin/python tools/specialize.py mechanisms/ts1.json $WORK/gen/MECH/sp_ts1_runK.c \
   --symbol ts1_specialized --batch 8 --split 32 --algebra-batch 8 --width 32 --tile 0 --fuse 0 \
   --joint-batch 1 --lu-order column)                                   # K = 0, 1, 2, timed
(cd sol-xhigh && make bench CC=clang OPTFLAGS="-O3 -DNDEBUG")           # timed
(cd runs/MECH/generic     && $WORK/sol-xhigh/bin/ts1_bench_cpu         --threads 1 4 60 60 CELLS)
(cd runs/MECH/specialized && $WORK/sol-xhigh/bin/ts1_bench_specialized --threads 1 4 60 60 CELLS)
python3 ~/NSF-NCAR/ODE/benchmarks/tier2/reference/compare_final.py A B   # S-G, G-MICM, S-MICM
# at the end, always: cp backup/ts1.json sol-xhigh/mechanisms/ts1.json
```
CELLS is `~/NSF-NCAR/tier2-work/exercises/canonical/MECH/cells.csv`. The logs are in `logs/MECH/`,
the numbers in `results/MECH.json`, and the runs in `runs/MECH/`.

After the restore, the TS1 configuration's SHA-256 matched HEAD's (`fef0ecdc…705a`).
`make bench` then rebuilt TS1 (exit 0, 18.1 s), and both programs passed the 100-cell TS1 check
against MICM again (difference 0). `git -C sol-xhigh status --porcelain` was empty.

The addendum used two more commands:
```
python3 -B run_cells64.py chapman carbon_bond_5 carbon_bond_5_no_emissions
# per mechanism: wait past the second of the newest build output; cp converted/MECH.json
# sol-xhigh/mechanisms/ts1.json; make bench CC=clang OPTFLAGS="-O3 -DNDEBUG"; check the rebuild;
# (cd runs64/MECH/generic     && $WORK/sol-xhigh/bin/ts1_bench_cpu         --threads 1 64 60 60 CELLS64)
# (cd runs64/MECH/specialized && $WORK/sol-xhigh/bin/ts1_bench_specialized --threads 1 64 60 60 CELLS64)
# compare_final.py: S-G, G-micm-run64, S-micm-run64 (and against micm-run64-g32); restore ts1.json
# timing, with the no_emissions build: 3 interleaved rounds of the 6 configurations above,
# 300 steps (diag/timing/results.txt)
```
`run_cells64.py` is 125 lines; the results are in `results/MECH-cells64.json`. TS1 was then
rebuilt again (exit 0, 16.8 s, every product rebuilt) and passed its check with difference 0.
The worktree was clean again.

## What it took

**The format.** `mechanisms/ts1.json` is CAMP JSON, with three top-level keys:
- `"concentration units": "mol m-3"`;
- `"initial state"`: temperature, pressure, `concentrations`, and `custom`, the label values;
- `"camp-data"`: CHEM_SPEC entries, then one MECHANISM entry with `reactions`.

The generic descriptors (`tools/configure.py`) and the generator (`tools/specialize.py`, which
imports configure's reader) both read this same file. The benchmark programs refer only to the
generated symbols `ts1` and `ts1_specialized`. They map CSV columns by name to temperature,
pressure, species, or rate-parameter labels, so any mechanism's cells.csv works unchanged.

**The converter** (`convert.py`, 118 lines, standard library only) maps each canonical feature
the way ts1.json expresses it:
- **Third body:** `third_body: true` becomes `"parameterization": "air-density"`. M stays in
  the reactant and product lists, as in ts1.json. The reader multiplies k by the air density
  once for each M reactant and drops M from the state, which matches MICM's `SetThirdBody()`.
- **Products:** `[name, yield]` becomes `{"species", "yield"}`.
- **Reactants:** the list is kept, repeats included, for stoichiometry.
- **ARRHENIUS and TROE:** the parameters are copied under the same keys, since TS1 uses the
  canonical names.
- **USER_DEFINED:** `label` becomes `"MUSICA name"`, and `scaling_factor` is kept.
- **Initial values:** concentrations and parameter values go into `"initial state"`. Values the
  canonical file omits are left out, and the reader defaults them to 0.
- **Units:** values pass through unchanged, since TS1's file also declares `mol m-3`.

**Nothing in the mechanisms needed a workaround.** The 14 emissions are USER_DEFINED reactions
with an empty reactant list. `configure.read` accepts them, and the C library and the generator
treat them as zero-order: forcing only, no Jacobian entries. Negative yields (PAR) and species
that are both reactant and product (PAR, ROR) pass through as well. carbon_bond_5_no_emissions
was therefore not needed as a fallback.

**The one obstacle was setup, not the entry's design.** This entry's `make env` creates `.conda`
by downloading its pinned packages from conda-forge, about 70 MB; the Tier 2 verification and
the A/B exercise logs show the download. The no-network rule rules that out. I satisfied
`make env` without network: a 4-line shell wrapper at the worktree's git-ignored
`.conda/bin/python` runs the existing environment that `make env` built from the identical
`environment.yml` during this entry's Tier 2 verification
(`~/NSF-NCAR/tier2-work/verify/sol-xhigh/clang/sol-xhigh/.conda`: Python 3.12.12, PyYAML 6.0.3).
The wrapper sets `PYTHONDONTWRITEBYTECODE=1`, so nothing is written there. The tools call no
subprocesses and no clang-format, so that environment is all they use.

**Surprises:**
- **Bitwise agreement.** The entry agrees with MICM bit for bit on all three new mechanisms, as
  it does on TS1. It reproduces MICM's species ordering (diagonal Markowitz with MICM's tie
  breaks), reaction ordering, and operation order closely enough that the LU rounding is the
  same.
- **Generation is cheap; compiling dominates.** The generator takes 0.05 to 0.08 s; a
  carbon_bond_5 build takes 4.7 s, almost all of it compiling 36,000 lines of generated C.
- **The scenario format assumes TS1's label naming.** This is outside the measured runs. The
  specialized scenario driver (`src/scenario.c:301`) scales the rate parameters whose labels
  begin with `j`, which is how the disclosed scenarios' README (line 49) defines photolysis
  parameters, following TS1's naming. A day-night scenario for these mechanisms, whose
  photolysis labels are `PHOTO_*`, would need the labels renamed or the convention generalized.
  This is the scenario format's convention, not a defect of the entry.
