# Exercise C, luna-high: the generator on another mechanism

Entry: `tier2/luna/high` at 2142282, detached worktree `$WORK/luna-high`
(`WORK=~/NSF-NCAR/tier2-work/exercises/c/luna-high`). Route: **the benchmark
programs**. A converted mechanism replaced `mechanisms/ts1.json` and was built with
`make bench`. No driver of my own, and no change to the entry.

## Results

clang, `-O3 -DNDEBUG`, 1 thread, the mechanism's 4 cells: one warm-up step and then
60 steps of 60 s, with MICM's default Rosenbrock settings (rtol 1e-6, atol 1e-3).
Comparisons use `compare_final.py` (rtol 1e-10). "PASS 0" means the largest relative
difference is zero.

| Mechanism | Converted | Generator ran | Generator (s) | Generated C (lines / bytes) | `make bench` exit / s | Generic ran | Specialized ran | Specialized vs generic | Generic vs MICM | Specialized vs MICM |
| --- | --- | --- | ---: | ---: | --- | --- | --- | --- | --- | --- |
| chapman | yes | yes | 0.027 | 301 / 8,993 | 0 / 2.94 | yes | yes | PASS 0 | PASS 0 | PASS 0 |
| carbon_bond_5 | yes | yes | 0.033 | 10,610 / 362,820 | 0 / 5.19 | yes | yes | PASS 0 | PASS 0 | PASS 0 |
| carbon_bond_5_no_emissions | yes | yes | 0.032 | 10,498 / 358,647 | 0 / 5.10 | yes | yes | PASS 0 | PASS 0 | PASS 0 |
| *TS1, for scale* | (as committed) | yes | 0.060 | 31,335 / 1,124,511 | 0 / 13.2–13.5 | yes | yes | PASS 0 | 1.32e-14 (100-cell input) | 1.32e-14 |

All three mechanisms reproduce MICM exactly:
- every value of every species in every cell is text-identical to MICM's `ts1_final.csv`;
- the generic and specialized outputs are byte-identical.

The runs are not trivial:
- 61 of CB5's 66 state species change during the run;
- with and without emissions, CB5 ends in different states (O3 2.10e-7 against
  1.71e-6 mol m-3), so the 14 emissions take effect.

Solver counters agree on both paths: 540 internal steps, 0 rejected. Timed region
(4 cells, so only indicative):

| Mechanism | Generic (ms) | Specialized (ms) |
| --- | ---: | ---: |
| chapman | 0.535 | 0.209 |
| carbon_bond_5 | 20.06 | 4.52 |
| carbon_bond_5_no_emissions | 20.29 | 4.25 |

How the measurements were taken:
- **Generator time:** the median of 5 runs of
  `python3 tools/generate_c.py IN OUT`, including Python startup (about 0.015 s).
- **Generated C lines:** counted with `wc -l`. The generator writes no final newline,
  so each file has one more line than the count.
- **Build time:** the wall time of `make bench`, which does three things in sequence:
  - it generates the C;
  - it runs three full clang compiles: the generic benchmark, the specialized
    benchmark, and the specialized scenario runner;
  - there are no object files, so each compile rebuilds every source.

The build gave no compiler warnings under `-Wall -Wextra -Wpedantic`. The generated
file in the build was byte-identical to a separate run of the generator outside the
tree.

## Addendum: 64 distinct cells

Setup:
- The build was the same as in the main runs: the converted mechanism swapped in
  for `mechanisms/ts1.json`, then `make bench CC=clang OPTFLAGS="-O3 -DNDEBUG"`.
- Each program ran
  `--threads 1 64 60 60 ~/NSF-NCAR/tier2-work/exercises/canonical/MECH/cells64.csv`
  in its own empty directory under `$WORK/runs64/MECH/`.
- The reference is `micm-run64/ts1_final.csv`. MICM's `micm-run64-g32/` file is
  identical, and comparing against it gives the same results.

| Mechanism | Specialized vs generic | Generic vs MICM | Specialized vs MICM | Generated kernels used | Generic / specialized (ms) |
| --- | --- | --- | --- | --- | ---: |
| chapman | PASS 0 | PASS 0 | PASS 0 | yes | 8.52 / 3.12 |
| carbon_bond_5 | PASS 0 | PASS 2.33e-14 (CLO, cell 58) | PASS 2.33e-14 | yes | 330.35 / 69.07 |
| carbon_bond_5_no_emissions | PASS 0 | PASS 2.06e-16 (SULF, cell 41) | PASS 2.06e-16 | yes | 314.36 / 68.30 |

- On every mechanism, the generic and specialized files are byte-identical.
- The solver counters agree on both paths: 540 internal steps, 0 rejected.
- `make bench` took 3.04, 5.20, and 4.91 s.

**The specialized program ran the generated kernels in both the 4-cell and the
64-cell runs.** No code path depends on the number of cells or on the group size:

1. **`src/main.c`** is shared by both programs. It solves the warm-up (line 584) and
   every timed step (line 605) with `chem_rosenbrock3_solve_group_with_workspace`.
   With `--threads 1` it passes one group of all the cells: 4, or 64.
2. **The group solver** (`src/solver.c` 804–983) works through the group one cell at
   a time. For each cell it calls `solver_forcing` (877, 900),
   `solver_negative_jacobian` (880), `solver_factor` (891), and
   `solver_triangular_solve` (914).
3. **Each of these wrappers** (`src/solver.c` 49–66, 365–381) calls
   `m->specialized_ops->…` whenever that pointer is non-NULL, and the generic
   routine otherwise. That pointer is the only switch.
4. **The pointer is set at build time.** `generated/ts1.c` sets
   `chem_generated_mechanism.specialized_ops = &generated_ops` under
   `#ifdef CHEM_SPECIALIZED`, and NULL otherwise. The Makefile passes
   `-DCHEM_SPECIALIZED` only to `ts1_bench_specialized` and
   `ts1_scenario_specialized`.
5. **The built binaries confirm it.** `$WORK/check_ops.py` reads each mechanism's
   build with `nm` and `dyld_info -fixups`, without running it:
   - in `ts1_bench_specialized`, `_chem_generated_mechanism+0x48` (the
     `specialized_ops` field) rebases to `_generated_ops`;
   - that table's four slots rebase to `_generated_forcing`,
     `_generated_negative_jacobian`, `_generated_factor`, and
     `_generated_solve_factored`;
   - `ts1_bench_cpu` has no `_generated_*` symbols and no fixup at that field, so
     its pointer is NULL.

   The output is saved in `$WORK/runs64/MECH/ops_*.txt`.
6. **The timings agree.** On carbon_bond_5 the specialized program is 4.4 times
   faster with 4 cells and 4.8 times faster with 64, with identical solver counters.

The "generated kernels" are the forcing, the Jacobian, the LU factorization, and
the triangular solve. Rate constants come from the shared generic evaluator in both
programs, by design (DESIGN.md).

This entry's benchmark has no 32-cell grouping:
- with `--threads 1`, the 64 cells form one group that shares one step controller,
  as in MICM's default run;
- only the scenario runner takes `--group-size`.

So no fallback to generic code is possible here.

## Commands

```sh
WORK=~/NSF-NCAR/tier2-work/exercises/c/luna-high
mkdir -p $WORK && cd $WORK
ln -sn ~/NSF-NCAR/ODE/benchmarks/micm micm
ln -sn ~/NSF-NCAR/ODE/benchmarks/reference reference
ln -sn ~/NSF-NCAR/ODE/benchmarks/tier2/scenarios scenarios
git -C ~/NSF-NCAR/TurboMICM worktree add --detach $WORK/luna-high tier2/luna/high
cd $WORK/luna-high
CONDA_OFFLINE=true make env    # fails: PackagesNotFoundError python=3.12.11 (see below)
make cpu bench CC=clang OPTFLAGS="-O3 -DNDEBUG"           # TS1: 13.5 s
# sanity, each in its own empty directory: both pass against
# micm_ts1_final_100_seed1.csv (1.32e-14) and are identical to each other
bin/ts1_bench_cpu         --threads 1 100 100 30 ../reference/ts1_cells_100_seed1.csv
bin/ts1_bench_specialized --threads 1 100 100 30 ../reference/ts1_cells_100_seed1.csv

cd $WORK
C=~/NSF-NCAR/tier2-work/exercises/canonical
python3 convert.py $C/MECH/canonical.json converted/MECH.json          # for each MECH
python3 luna-high/tools/generate_c.py converted/MECH.json gen/MECH.c   # timed, 5 runs
cp -p luna-high/mechanisms/ts1.json backup/ts1.json                    # = HEAD's blob
bash run_exercise.sh chapman
bash run_exercise.sh carbon_bond_5 carbon_bond_5_no_emissions
#   for each MECH:
#     cp converted/MECH.json luna-high/mechanisms/ts1.json
#     (cd luna-high && /usr/bin/time -p make bench CC=clang OPTFLAGS="-O3 -DNDEBUG")
#     (cd runs/MECH/generic     && $WORK/luna-high/bin/ts1_bench_cpu         --threads 1 4 60 60 $C/MECH/cells.csv)
#     (cd runs/MECH/specialized && $WORK/luna-high/bin/ts1_bench_specialized --threads 1 4 60 60 $C/MECH/cells.csv)
#     compare_final.py specialized generic; generic MICM; specialized MICM
#   on exit: cp backup/ts1.json luna-high/mechanisms/ts1.json, and
#   `git diff --quiet` confirms it matches HEAD
cd luna-high && make bench CC=clang OPTFLAGS="-O3 -DNDEBUG"   # TS1 again (13.2 s); sanity passes again
git -C $WORK/luna-high status --porcelain                      # empty

# addendum: 64 distinct cells
cd $WORK && bash run_cells64.sh
#   for each MECH: the same swap and `make bench` as above, then
#     python3 check_ops.py luna-high/bin/ts1_bench_specialized   (and bin/ts1_bench_cpu)
#     (cd runs64/MECH/generic     && $WORK/luna-high/bin/ts1_bench_cpu         --threads 1 64 60 60 $C/MECH/cells64.csv)
#     (cd runs64/MECH/specialized && $WORK/luna-high/bin/ts1_bench_specialized --threads 1 64 60 60 $C/MECH/cells64.csv)
#     compare_final.py specialized generic; generic micm-run64; specialized micm-run64; specialized micm-run64-g32
#   on exit: TS1 restored (matches HEAD)
cd luna-high && make bench CC=clang OPTFLAGS="-O3 -DNDEBUG"   # TS1 rebuilt; sanity passes again
git -C $WORK/luna-high status --porcelain                      # empty
```

Logs and outputs are in `$WORK/runs/MECH/` (4 cells) and `$WORK/runs64/MECH/` (64 cells, plus
`ops_*.txt`): `build.log`, `generic.log`, `specialized.log`, and
`{generic,specialized}/ts1_final.csv`. The converted
mechanisms are in `$WORK/converted/` and the generated C in `$WORK/gen/`.

## What it took

**Finding the input.**
- The Makefile's `generated/ts1.c` rule hard-codes `mechanisms/ts1.json`, with no
  variable, so the exercise swapped that file and then restored it.
- The format is the entry's normalized `"c-mechanism-1"` JSON. Two of its tools
  write it:
  - `tools/import_ts1.py` parses MICM's `ts1_mechanism.hpp` but hard-codes 210
    species and 547 reactions;
  - `tools/read_config.py` reads CAMP configurations. The tests use it to run the
    generator on MICM's tutorial mechanisms.
- canonical.json maps directly onto the normalized JSON, so neither tool was needed.

**The converter** (`convert.py`, 102 lines, 70 of them code) maps:
- **Species:** `{name, third_body, molecular_weight: 0.0, diffusion: 0.0}`, in the
  source's order. Only surface reactions use molecular weight and diffusion, and many
  TS1 species have 0.0 for both. M stays a species with `"third_body": true`, as in
  TS1.
- **Types:**
  - `ARRHENIUS` becomes `arrhenius` (A, B, C, D, E);
  - `TROE` becomes `troe` (k0_A, k0_B, k0_C, kinf_A, kinf_B, kinf_C, Fc, N);
  - `USER_DEFINED` becomes `user_defined` (label, scaling_factor).

  Types are written in lower case, as in TS1; the generator accepts either case.
  Parameters go in a nested `parameters` object, with `scaling_factor` on every
  reaction, as `import_ts1.py` writes them.
- **Terms:**
  - each listed reactant becomes one `{name, coefficient: 1.0}`, with repeats kept
    and M included, as in TS1;
  - each product `[name, yield]` becomes `{name, coefficient: yield}`;
  - an emission becomes `"reactants": []`.
- **Initial values:** `defaults.species` lists every non-third-body species, with 0.0
  where none is given, as TS1's InitState does. `defaults.rate_parameters` lists every
  label, followed by temperature and pressure.
- No unit conversion was needed: both formats use mol m-3.

**What was surprising or worth noting:**
- **Nothing in the run path is TS1-specific.**
  - `src/main.c` maps the CSV's columns by name to the generated species and custom
    labels; its error text still says "unknown TS1 CSV column".
  - The output columns come from the generated species table, with the third body
    left out.
- **"Generic" still means generated tables.**
  - `bin/ts1_bench_cpu` compiles the same `generated/ts1.c` without
    `CHEM_SPECIALIZED`. That is the table-driven Tier 1 solver over generated
    tables, not a mechanism read at run time.
  - So both programs need regeneration and a rebuild for each new mechanism.
- **The generator needed no special cases.**
  - Emissions (no reactants), negative yields (CB5's PAR), repeated reactants, and
    reactions with M on both sides all went through as written.
  - CB5's five labels with no initial value came through as 0.0 from
    canonical.json and cells.csv.
- **The results are bit-identical to MICM.** The generator:
  - reorders species by its replica of MICM's Markowitz ordering, which keeps
    MICM's `size_t` underflow for product-only species (per its comment);
  - orders reactions by rate-constant type;
  - keeps the generic library's accumulation order.

  That is the likely reason for the exact agreement. TS1 itself differs from MICM
  by 1.32e-14, so exact agreement is not guaranteed.
- **`make env` needs network access.**
  - It creates `.conda/` from conda-forge with a fresh package cache in
    `.conda-pkgs/`. Forced offline, it fails at once (PackagesNotFoundError,
    python=3.12.11), and I did not retry it with network access.
  - The generator uses only the standard library. The Makefile falls back to
    `python3` (3.13.1 here) when `.conda/bin/python` is absent, so `make bench`
    worked unchanged.
  - Only `tests/test_specialized.py`, which calls `.conda/bin/python` directly, would
    need the environment.
  - I removed the empty `.conda-pkgs/` that the attempt left behind.
- **Side check: a mechanism with no user-defined rates** (`$WORK/scratch-no-labels/`,
  not one of the three).
  - The generator emits empty `custom_labels[]` and `default_custom_parameters[]`
    arrays, which give 4 `-Wpedantic` warnings: C23 empty initializers and
    zero-length arrays.
  - `src/main.c`'s size check, `cells > (size_t)-1 / nr`, then divides by zero,
    which is undefined behavior.
  - The clang -O3 build still ran and matched the exact solution, exp(-0.66).

**The entry is unchanged.**
- `mechanisms/ts1.json` was restored from the backup and matches HEAD.
- The TS1 binaries were rebuilt, and the 100-cell sanity check passes again.
- `git status --porcelain` is empty. Only the ignored `bin/` and `generated/` remain.
