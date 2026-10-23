# Exercise C, astra-max: the generator on other mechanisms

Entry: `tier2/astra/max` (cf16213), detached worktree at `$WORK/astra-max`
(`WORK=~/NSF-NCAR/tier2-work/exercises/c/astra-max`). Route: **benchmark programs**.
The configuration the build reads, `configs/ts1.json`, was swapped for each converted
mechanism, then `make bench` was run, then `bin/ts1_bench_cpu` and
`bin/ts1_bench_specialized`. The original file was restored afterwards (SHA-256 `73843ae0…e65e`,
identical to the backup), and `git status --porcelain` shows no tracked change. Every
mechanism converted, generated, compiled without warnings, solved, and passed every
comparison at 1e-10 relative. The specialized output is bitwise identical to the generic
output in every run. No change to the entry was needed.

## Results: 4 identical cells (`cells.csv`), 1 thread, 1 warm-up step, then 60 × 60 s

| Mechanism | Evolving species / reactions | Converted | Generator (s) | Generated C: lines / bytes | `make bench`: exit / s | Generic ran | Specialized ran | Specialized vs generic | Generic vs MICM | Specialized vs MICM | Generated kernels used |
|---|---|---|---:|---:|---|---|---|---|---|---|---|
| chapman | 4 / 7 | yes | 0.035 | 863 / 25,381 | 0 / 1.12 | yes | yes | PASS 0 | PASS 6.28e-16 | PASS 6.28e-16 | yes (runtime-stride form) |
| carbon_bond_5 | 66 / 200 | yes | 0.048 | 36,844 / 915,508 | 0 / 5.16 | yes | yes | PASS 0 | PASS 1.77e-15 | PASS 1.77e-15 | yes (runtime-stride form) |
| carbon_bond_5_no_emissions | 66 / 186 | yes | 0.048 | 36,454 / 900,047 | 0 / 5.15 | yes | yes | PASS 0 | PASS 5.66e-15 | PASS 5.66e-15 | yes (runtime-stride form) |
| *TS1, for scale* | *209 / 547* | – | *0.169* | *114,209 / 2,863,675* | *0 / 14.9–15.4* | | | | | | |

How each column was measured:

- **Comparisons:** `compare_final.py`, giving the largest relative difference. They cover 24 values (4 cells × 6 columns) for chapman and 272 (4 × 68) for each CB5 variant. M is not a column: like MICM, the entry holds no concentration for M.
- **Solver statistics:** both programs report the same counts, 540 internal steps with no rejections.
- **Elapsed time, one run each at 1 thread (not a timing session):**
  - chapman: generic 0.224 ms, specialized 0.091 ms;
  - carbon_bond_5: 10.61 and 3.14 ms;
  - carbon_bond_5_no_emissions: 9.99 and 3.19 ms.
- **Generator time:** the median of 5 runs of the exact command the Makefile runs (`tools/specialize.py … --block 32 --fixed 1 --lu mozart --fusion 4`). It includes Python start-up and imports, about 27 ms.
- **Generated C:** `tools/specialize.py`'s three outputs, `ts1_specialized.c` + `ts1_specialized_kernels.c` + `ts1_specialized.h`. The build also generates the generic constructor `build/generated/ts1.c` with `tools/configure.py`: 39 lines (2,670 bytes) for chapman, 392 (51,900) for carbon_bond_5, 350 (47,627) for no_emissions, and 1,543 (185,815) for TS1. Per-file sizes:
  - chapman: 164 + 689 + 10 lines;
  - carbon_bond_5: 937 + 35,897 + 10;
  - no_emissions: 887 + 35,557 + 10.

  The `_kernels.c` file dominates because every kernel is emitted twice, in a runtime-stride form and a fixed-width form. The generator's own metadata reports 751 Jacobian nonzeros and 2,501 LU operations for both CB5 variants, and 14 and 18 for chapman.
- **`make bench` time:** the median of 3 builds after the swap. Only the mechanism-dependent targets rebuild: `ts1.c`, `ts1.o`, `ts1_bench_cpu`, `ts1_specialized.c`, its two objects, and the two specialized programs. Individual builds:
  - chapman: 1.118, 1.159, 0.918 s;
  - carbon_bond_5: 5.213, 5.155, 4.998 s;
  - no_emissions: 5.263, 5.151, 5.091 s;
  - TS1, the same incremental rebuild: 15.40, 14.89, 15.17 s;
  - the first clean `make cpu bench` of the fresh worktree, library included: 16.75 s.

## Addendum: 64 distinct cells (`cells64.csv`), two full 32-cell groups

The entry was built exactly as in the main runs: the swap, then `make bench` with the Makefile's defaults (`GROUP_SIZE=32`, `KERNEL_WIDTH=32`). Each log shows all 8 targets rebuilt. MICM's result is `micm-run64/ts1_final.csv`, byte-identical to `micm-run64-g32/`.

| Mechanism | Specialized vs generic | Generic vs MICM | Specialized vs MICM | Generated kernels used | Elapsed, generic / specialized (ms) |
|---|---|---|---|---|---|
| chapman | PASS 0 | PASS 2.46e-15 | PASS 2.46e-15 | yes (fixed-width `*_fixed` form) | 1.125 / 0.593 |
| carbon_bond_5 | PASS 0 | PASS 1.82e-14 | PASS 1.82e-14 | yes (fixed-width `*_fixed` form) | 55.03 / 40.81 |
| carbon_bond_5_no_emissions | PASS 0 | PASS 2.69e-14 | PASS 2.69e-14 | yes (fixed-width `*_fixed` form) | 54.26 / 40.66 |

The comparisons cover 384 values for chapman and 4,352 for each CB5 variant. Both programs report the same counts: 1,080 internal steps (2 groups × 540) with no rejections.

### Did the specialized program run the generated kernels? Yes, in both cases

This is established from the code path, in the source and in the built binaries. There is no runtime trace: lldb needs developer mode, which is disabled on this machine, and enabling it would prompt for admin credentials.

1. **The program calls the generated solver.** `bin/ts1_bench_specialized` is compiled with `-DCH_SPECIALIZED`, so `tools/ts1_io.h:12-18` maps `ch_update_rates` to `ch_rates_ts1_specialized` and `ch_solve` to `ch_solve_ts1_specialized`.
2. **The generated solver uses generated kernels.** The generated `ts1_specialized.c` ends with `#define CH_FORCE ch_force_ts1_specialized` (and likewise for `CH_JACOBIAN`, `CH_FACTOR`, `CH_LINEAR`) and then `#include "specialized_solver.h"`; for chapman these are lines 158–164. So the Rosenbrock loop calls only generated functions, apart from the library's `ch_normalized_error`, `ch_post_solve`, and the step-size controller.
3. **Each dispatcher picks one of two generated forms.** There are six dispatchers (rates, forcing, Jacobian, shift, factor, linear), emitted by `Emitter.kernel` in `tools/specialize.py`. Each tests `if (n == CH_DEFAULT_GROUP_SIZE && n % CH_KERNEL_WIDTH == 0)`, where n is the group's cell count. If the test holds, it calls the `*_fixed` blocks, compiled with the constant stride 32. Otherwise it calls the runtime-stride blocks `ts1_specialized_<kind>_kernel_<k>(n, …)`. Both forms are generated and specific to the mechanism, with literal indices and coefficients. There is no branch back to the generic library.
4. **4 cells take the runtime-stride form.** `tools/ts1_bench.c:23` sets `if (group_size > cells) group_size = cells;`, which gives one group of n = 4, so every dispatcher takes the runtime-stride branch.
5. **64 cells take the fixed-width form.** They form two groups of n = 32, which equals both `CH_DEFAULT_GROUP_SIZE` and `CH_KERNEL_WIDTH` (the Makefile passes `-D…=32`), so every dispatcher takes the fixed-width branch.
6. **The binaries agree.** In the built programs (`otool -tV`), `_ch_force_kernel_ts1_specialized` compiles to `cmp x22, #0x20; b.ne <runtime-stride blocks>`, and otherwise calls the `*_fixed` blocks: 1 block for chapman, 7 for carbon_bond_5, 6 for no_emissions. `_ch_solve_ts1_specialized` calls only the generated dispatchers and kernels, plus `_ch_normalized_error` and `_ch_post_solve`. It never calls `_ch_force`, `_ch_jacobian`, or `_ch_solve`, although the generic library is linked in for the model's setup.
7. **Timing is consistent with this.** At 4 cells the specialized program takes 3.1 ms against 10.6 ms for generic, so it is not running the generic code.

## Extra checks (not required)

- **64 identical cells**, built from `cells.csv` repeated 16 times and compared with MICM's 4-cell file repeated 16 times: everything passes. Specialized vs generic is 0 for all three; against MICM the largest differences are 6.28e-16, 1.77e-15, and 5.66e-15.
- **The entry's own generator test**, `tests/specialized.c`, built for the TS1 slot by `make bin/test_specialized_ts1`, run after each swap.
  - **What it checks:** it compares generic and generated rates, forcing, Jacobian, LU, linear solve, and adaptive trajectories at 1e-10, at 1, 3, 7, 32, and 33 cells, with default and tight tolerances.
  - **chapman:** PASS.
  - **carbon_bond_5:** 1 cell passed with both tolerance settings. At 3 cells with default tolerances, both the generic and the generated solver returned "Convergence exceeded maximum steps", with identical counts (1002 steps, 227 rejections each).
  - **carbon_bond_5_no_emissions:** at 1 cell with default tolerances, both failed the same way with identical counts (1001 steps, 194 rejections each).
  - **Interpretation:** no kernel mismatch was reported before the failing solves, so this is a limit of the test's synthesized inputs, not a disagreement between generic and generated code. The test perturbs temperature, pressure, and density, scales parameters by 0.93 per step, and, for mechanisms of at most 100 species, sets every rate parameter that is 0 or ≥ 1 to 0.001 s⁻¹. That last rule turns on PHOTO_O2, PHOTO_CL2, PHOTO_N2O5, PHOTO_PAN, and PHOTO_PANX, which is the likely cause (not verified).
- **TS1 sanity check**, before the swaps and after the final TS1 rebuild: 100 cells of seed 1. Specialized equals generic bitwise, and both are within 3.23e-14 of MICM.

## Commands

```sh
WORK=~/NSF-NCAR/tier2-work/exercises/c/astra-max
# setup, as specified
mkdir -p $WORK && cd $WORK && ln -sn …/benchmarks/micm micm && ln -sn …/benchmarks/reference reference && ln -sn …/tier2/scenarios scenarios
git -C ~/NSF-NCAR/TurboMICM worktree add --detach $WORK/astra-max tier2/astra/max
cd $WORK/astra-max && make env && make cpu bench CC=clang OPTFLAGS="-O3 -DNDEBUG"   # 11.4 s and 16.75 s
# conversion
.conda/bin/python $WORK/convert.py ~/NSF-NCAR/tier2-work/exercises/canonical/MECH/canonical.json $WORK/converted/MECH.json
# generator timing and size, outside the tree (identical to what make generated: cmp)
.conda/bin/python tools/specialize.py $WORK/converted/MECH.json $WORK/gen/MECH/ts1_specialized.c --name ts1_specialized --block 32 --fixed 1 --lu mozart --fusion 4
# main runs: $WORK/run_exercise.sh. It backs up configs/ts1.json to $WORK/backup/ts1.json (cp -p, SHA-256 recorded),
# then for each MECH:
cp $WORK/converted/MECH.json configs/ts1.json
make bench CC=clang OPTFLAGS="-O3 -DNDEBUG"          # timed
(cd $WORK/runs/MECH/generic     && …/bin/ts1_bench_cpu         --threads 1 4 60 60 …/canonical/MECH/cells.csv)
(cd $WORK/runs/MECH/specialized && …/bin/ts1_bench_specialized --threads 1 4 60 60 …/canonical/MECH/cells.csv)
python3 …/tier2/reference/compare_final.py A B       # specialized/generic, generic/MICM, specialized/MICM
# and restores configs/ts1.json from the backup on exit (trap), checked with cmp
# addendum: $WORK/run_cells64.sh (same, with --threads 1 64 60 60 …/MECH/cells64.csv, against micm-run64)
# extras: run_fixed.sh (64 identical cells), run_harness.sh (make bin/test_specialized_ts1), time_builds.sh (repeat timings, TS1 rebuild)
```

## What it took

- **The converter** (`convert.py`, 124 lines) is a near one-to-one mapping, because the entry reads MICM's own camp-data JSON, the format of its `configs/ts1.json`:
  - **Species:** `{"name", "type": "CHEM_SPEC"}`; the third body M gets `"parameterization": "air density"`, as TS1's M does. Species keep the canonical order, which is MICM's header order, just as TS1's file keeps its header's order.
  - **Products:** `[name, yield]` becomes `{"species", "yield"}`.
  - **USER_DEFINED:** `label` becomes `"MUSICA name"`, and `scaling_factor` becomes `"scaling factor"`.
  - **ARRHENIUS and TROE:** the parameters keep their names.
  - **M:** stays in reactant and product lists, as in TS1. The reader multiplies the rate by the air density for each M reactant and drops M from the products.
  - **Emissions:** stay USER_DEFINED with `"reactants": []`, which the reader and generator accept as they are. The reader also knows an `EMISSION` type, which TS1 does not use.
  - **Rate units:** `"rate units": "mol m-3"`, so the reader does no molecule cm⁻³ scaling.
  - **Initial state:** temperature, pressure, concentrations, and parameters. The canonical `environment` is dropped: the programs recompute air density as P/(RT), and the time step is an argument.
  - **Dropped fields:** the canonical `source` field has no counterpart.
  - **Error handling:** any unknown reaction type or field raises an error. None occurred.
- **The generator is mechanism-generic, as claimed.** `tools/specialize.py` contains no TS1 names. It handled a 4-species mechanism, the CB5 variants with zero-reactant reactions, and negative yields (PAR −0.66), with no special handling. It is fast (35–48 ms), but the C is large: about 37,000 lines and 0.9 MB for CB5.
- **The configuration path is fixed.** The Makefile hard-codes `configs/ts1.json` in the rules for `build/generated/ts1.c` and `ts1_specialized.c`; there is no variable. The only way in was the swap.
- **Pitfall: one-second timestamps.** macOS's GNU Make 3.81 compares mtimes to the second, so a configuration swap in the same second as the previous build's outputs is treated as up to date. This happened twice in my extra checks. In one, a test build did nothing, so the first "PASS" for carbon_bond_5 came from chapman's binary. In the other, a timing build regenerated both C files and relinked the generic program for the new mechanism, but did not recompile `ts1_specialized_kernels.o` or relink the specialized programs, which stayed on the previous mechanism. Both results were discarded (kept in `logs/invalid-same-second/`). All scripts now pause 2 s before each swap, and every reported result comes from a build whose log shows all 8 mechanism-dependent targets rebuilt. The main 4-cell runs were complete builds, which was checked from the logs and by `cmp` of the generated C against separate generator runs. Anyone switching mechanisms quickly with this Makefile can get stale or mixed binaries without any message.
- **Other observations:**
  - `make env` downloaded the pinned packages (python 3.12.12, numpy, pyyaml) from conda-forge into the worktree; it was part of the prescribed setup.
  - The retained `tools/translate_ts1.py` (C++ header to configuration) is TS1-only: it asserts 210 species and 547 reactions. The reader and generator do not depend on it.
  - The scenario driver treats only parameters whose names start with `j` as photolysis, following TS1's naming, so MUSICA's `PHOTO_*` labels would not follow a day–night curve. The scenario driver was not used here.

## Files in `$WORK`

- **Converter and report:** `convert.py` (converter), `make_report.py` (builds `exercise-c.json` from the run files).
- **Scripts:** `run_exercise.sh`, `run_cells64.sh`, `run_fixed.sh`, `run_harness.sh`, `time_builds.sh`. The 2 s pause before each swap was added after the first runs.
- **Inputs and generated code:**
  - `converted/MECH.json`: the converted mechanisms;
  - `backup/ts1.json` (+ `.sha256`): the backup of the original configuration;
  - `gen/MECH/`: the generator's output, with the build's own copy in `from-make/`.
- **Runs:**
  - `runs/MECH/` (4 cells), `runs/MECH-cells64/` (the addendum, with the binaries, generated C, and disassembly), `runs/MECH-64cells/` (identical cells), and `runs/sanity-ts1*/`;
  - in each run directory, `*/ts1_final.csv`, `stdout.txt`, and `compare-*.txt`.
- **Logs:** `logs/`, one file per build and test, plus `logs-*.txt` with each script's console output.
