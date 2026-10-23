# Exercise C (another mechanism): sol-max

Entry: `tier2/sol/max` (commit ca0119e), detached worktree `sol-max/`. clang, `-O3 -DNDEBUG`, one thread.

**Route: benchmark programs.** `configs/ts1.json` is the file the build reads (`CONFIG_ts1` in the Makefile). Each converted mechanism took its place, then `make bench`, then `bin/ts1_bench_cpu` and `bin/ts1_bench_specialized` ran with the mechanism's `cells.csv`. No entry change was needed and no driver was written. The original file was restored at the end.

## Results

| | chapman | carbon_bond_5 | carbon_bond_5_no_emissions |
|---|---|---|---|
| Species (variables) / reactions / LU nonzeros | 5 (4) / 7 / 14 | 67 (66) / 200 / 742 | 67 (66) / 186 / 742 |
| Converted | yes | yes | yes |
| Generator ran: `specialize.py` time (s) | 0.029 | 0.071 | 0.070 |
| Generic-loader generator, `configure.py` (s) | 0.029 | 0.044 | 0.044 |
| `generated/ts1_specialized.c` lines / bytes | 902 / 29,827 | 22,867 / 732,601 | 22,488 / 722,785 |
| `generated/ts1.c` (generic loader) lines / bytes | 250 / 7,198 | 7,783 / 215,193 | 7,527 / 207,794 |
| `make bench`: exit / wall time (s) | 0 / 0.95 | 0 / 6.94 | 0 / 7.46 |
| Compile of `ts1_specialized.c` alone (s) | 0.17 | 4.87 | 4.84 |
| Generic ran / specialized ran | yes / yes | yes / yes | yes / yes |
| Specialized vs generic | PASS 0 | PASS 0 | PASS 0 |
| Generic vs MICM | PASS 0 | PASS 4.38e-15 (CLO) | PASS 1.4e-15 (HONO) |
| Specialized vs MICM | PASS 0 | PASS 4.38e-15 (CLO) | PASS 1.4e-15 (HONO) |
| Timed solve, generic / specialized (ms) | 0.148 / 0.073 | 6.360 / 1.643 | 6.677 / 1.666 |

TS1 itself, for scale: `specialize.py` takes 0.834 s and `configure.py` 0.308 s. They write 66,075 lines (2,222,161 bytes) and 21,455 lines (611,548 bytes). After the restore, `make bench` took 29.6 s, 23.7 s of it compiling `ts1_specialized.c`; the first full build (`make cpu bench`) took 33.7 s.

Notes on the measurements:
- Generator times are medians of 5 runs, measured with `perf_counter` around the Makefile's own command lines and written outside the tree. The files they wrote are byte-identical to the ones `make bench` then wrote in `generated/`.
- Comparisons used `compare_final.py` (rtol 1e-10, columns matched by name, all 6 or 68 columns present). "PASS 0" means bit-identical.
- Every run had 540 internal steps and no rejections, in both programs.
- The gate discriminates: carbon_bond_5's result compared with MICM's no-emission state fails on 244 of 272 values. So the 14 emissions are really exercised, and the solution moves well away from its initial state (59 of 66 species change).

## Addendum: 64 distinct cells

The build was exactly as in the main runs: the same swap of `configs/ts1.json`, `make bench` with the Makefile's defaults (0.94, 6.93, and 7.13 s), and generated C byte-identical to the main runs'. Both programs ran with `--threads 1 64 60 60 cells64.csv`, each in its own empty directory.

| 64 cells (two 32-cell groups) | chapman | carbon_bond_5 | carbon_bond_5_no_emissions |
|---|---|---|---|
| Specialized vs generic | PASS 0 | PASS 0 | PASS 0 |
| Generic vs MICM (`micm-run64`) | PASS 0 | PASS 8.46e-15 (XO2N, cell 53) | PASS 2.85e-14 (ETH, cell 27) |
| Specialized vs MICM (`micm-run64`) | PASS 0 | PASS 8.46e-15 | PASS 2.85e-14 |
| Generated kernels used: 4 cells / 64 cells | yes / yes | yes / yes | yes / yes |
| Timed solve, generic / specialized (ms) | 0.874 / 0.929 | 45.382 / 34.368 | 45.171 / 33.877 |

- `micm-run64-g32/ts1_final.csv` is byte-identical to `micm-run64/ts1_final.csv` and gives the same results.
- Every run had 1080 internal steps (540 per group) and no rejections, in both programs.

**Did the specialized program run the generated kernels?** Yes, at both sizes. sol-max's specialized program has no fallback to generic code at all.
- `src/bench.c`, built with `-DCHEM_SPECIALIZED` (lines 11-18), maps the loader, `chem_calculate_rates`, and `BENCH_SOLVE` to `chem_load_ts1_specialized`, `chem_rates_ts1_specialized`, and `chem_solve_ts1_specialized`.
- All three are in `generated/ts1_specialized.c`, which calls none of the library's numerical routines: grep finds no `chem_solve`, `chem_lu_mozart`, `chem_linear_solve`, or `chem_calculate_rates`.
- Each kernel exists in two generated variants, chosen by `if (nc == 4 && s == 4)`:
  - a literal-stride "fixed" 4-cell kernel;
  - a runtime-stride "dynamic" kernel.

  Both are emitted by `tools/specialize.py` from the mechanism. In carbon_bond_5's file (kept as `$WORK/gen/carbon_bond_5/ts1_specialized.c`), the choices are made at:
  - line 6782 (force);
  - line 9810 (Jacobian);
  - line 22724 (factorization, in `ts1_specialized_attempt`);
  - line 22634 (triangular solves).
- `chem_solve_ts1_specialized` (from `tools/rosenbrock3_blocked.c.in`) does two things:
  - once per internal step, it computes the group's forcing with `ts1_specialized_force(model, cells, stride)`;
  - it cuts the group into 4-cell tiles (`tile = 4`, `q.stride = 4`), and on each tile it runs the Jacobian, factorization, three solves, and stage forcing.
- **4 cells:** one group of 4 is exactly one tile. The group forcing call has nc = s = 4, and the single tile has q.cells = q.stride = 4, so every call takes the fixed generated branch.
- **64 cells:** two groups of 32. The group forcing (nc = s = 32) takes the dynamic generated branch (`ts1_specialized_force_0..N`). Each group's eight full tiles take the fixed generated branch for everything else.
- **Not generated:** only the shared step-size controller and Rosenbrock coefficients, state allocation, and CSV input and output.
- **No runtime trace:** macOS developer mode is off, so lldb would raise an authorization dialog. Identical outputs cannot tell the paths apart, because the generator keeps the generic code's arithmetic order.
- **Timing (CB5 without emissions, 5 repeats, Tier 1 program against specialized):**
  - 4 cells: 6.42-6.51 ms against 1.66-1.90 ms (3.8 times);
  - 64 cells: 44.9-45.3 ms against 33.7-34.2 ms (1.33 times).

  Per cell-step, the Tier 1 program gets cheaper with 32-cell groups (2.98 to 1.31 µs), because its group-wide loops amortize per-call overhead. The specialized program gets slightly dearer (0.78 to 0.98 µs). The 64-cell ratio is close to the entry's own 1.2 to 1.3 times on TS1 with 100 cells. Chapman's 4 species leave both programs overhead-bound, and there the two are about even.

## Commands

```sh
WORK=~/NSF-NCAR/tier2-work/exercises/c/sol-max; CAN=~/NSF-NCAR/tier2-work/exercises/canonical
CMP=~/NSF-NCAR/ODE/benchmarks/tier2/reference/compare_final.py
# setup as specified (links, worktree), then
cd $WORK/sol-max && make env && make cpu bench CC=clang OPTFLAGS="-O3 -DNDEBUG"      # 10.1 s, 33.7 s
# sanity: both programs on ../reference/ts1_cells_100_seed1.csv: PASS 3.48e-14 against MICM

python3 $WORK/convert.py $CAN/$m/canonical.json $WORK/converted/$m.json            # each m

# generator timing (outside the tree; the Makefile's command lines)
.conda/bin/python tools/configure.py  $WORK/converted/$m.json $WORK/gen/$m/ts1.c --name ts1 --kernel-chunk 32
.conda/bin/python tools/specialize.py $WORK/converted/$m.json $WORK/gen/$m/ts1_specialized.c --name ts1_specialized \
  --kernel-chunk 32 --linear-chunk 1 --fuse-kernels 1 --fuse-stages 1 --specialize-rates 1 --linear unrolled \
  --linear-loop rows --vector-error 1 --cell-stride 4 --initialize-kernels 1 --solver-tile 4

cp configs/ts1.json $WORK/backup/ts1.json        # SHA-256 9e14c904...48b6; git blob 88a9da89
# for m in chapman carbon_bond_5 carbon_bond_5_no_emissions:
cp $WORK/converted/$m.json configs/ts1.json
/usr/bin/time -p make bench CC=clang OPTFLAGS="-O3 -DNDEBUG"
(cd $WORK/runs/$m/generic     && $WORK/sol-max/bin/ts1_bench_cpu         --threads 1 4 60 60 $CAN/$m/cells.csv)
(cd $WORK/runs/$m/specialized && $WORK/sol-max/bin/ts1_bench_specialized --threads 1 4 60 60 $CAN/$m/cells.csv)
python3 $CMP runs/$m/specialized/ts1_final.csv runs/$m/generic/ts1_final.csv
python3 $CMP runs/$m/generic/ts1_final.csv     $CAN/$m/micm-run/ts1_final.csv
python3 $CMP runs/$m/specialized/ts1_final.csv $CAN/$m/micm-run/ts1_final.csv

cp $WORK/backup/ts1.json configs/ts1.json; git status --porcelain     # empty
make bench CC=clang OPTFLAGS="-O3 -DNDEBUG"      # 29.6 s; TS1 gate PASS 3.48e-14 again, both programs

# addendum, for each m: the same swap and make bench, then
(cd $WORK/runs/$m/cells64-generic     && $WORK/sol-max/bin/ts1_bench_cpu         --threads 1 64 60 60 $CAN/$m/cells64.csv)
(cd $WORK/runs/$m/cells64-specialized && $WORK/sol-max/bin/ts1_bench_specialized --threads 1 64 60 60 $CAN/$m/cells64.csv)
python3 $CMP runs/$m/cells64-specialized/ts1_final.csv runs/$m/cells64-generic/ts1_final.csv
python3 $CMP runs/$m/cells64-generic/ts1_final.csv     $CAN/$m/micm-run64/ts1_final.csv
python3 $CMP runs/$m/cells64-specialized/ts1_final.csv $CAN/$m/micm-run64/ts1_final.csv   # and micm-run64-g32
# then restored again: SHA-256 unchanged, git status empty, make bench 30.9 s, TS1 gate PASS 3.48e-14
```

Logs and outputs are in `$WORK`:
- `logs/`: the make logs, which hold only timings because the Makefile is silent, and before/after mtimes;
- `gen/`: the timed generator outputs, `generator_times.json`, and compile-only objects;
- `runs/<mechanism>/{generic,specialized}/`: each 4-cell run's `ts1_final.csv` and stdout;
- `runs/<mechanism>/{cells64-generic,cells64-specialized}/`: the same for each 64-cell run (`runs/carbon_bond_5_no_emissions/repeat-timing/` holds the repeat timings);
- `converted/`: the converted mechanisms;
- `backup/ts1.json`: the original configuration.

## What it took

**Format.** sol-max reads CAMP-style YAML or JSON through `tools/configure.py`, and `tools/specialize.py` shares that reader and the symbolic layout. `configs/ts1.json` is one JSON file. It holds `"concentration units": "mol m-3"`, the `camp-data` records (210 `CHEM_SPEC` and one `MECHANISM` with 547 reactions), and an `initial state`. `convert.py` (90 lines) maps each canonical field to TS1's own key and rejects any field it does not map:

| Canonical | sol-max (as in `configs/ts1.json`) |
|---|---|
| species `{name, third_body}` | `{"name", "type": "CHEM_SPEC"}`; M gets `"parameterization": "third body"`, exactly TS1's M record. TS1's molecular weights are only used by SURFACE reactions, and the canonical has none, so they are omitted (optional, as in the entry's Robertson example) |
| `reactants` (repeats) | `"reactant sequence"`, M kept as in TS1's four M reactions: the library multiplies k by air density once per third-body reactant, and ignores M as a product |
| `products` `[name, yield]` | `"product sequence"` `[{"species", "yield"}]` |
| ARRHENIUS A–E; TROE k0_A–C, kinf_A–C, Fc, N | the same keys |
| USER_DEFINED `label`, `scaling_factor` | `"MUSICA name"`, `"scaling factor"` |
| emission (no reactants) | USER_DEFINED with `"reactant sequence": []`. TS1 has none; the reader and the C loader accept zero reactants. The reader's own `EMISSION` type was not used, so as to mirror TS1 |
| units already mol m-3 | top-level `"concentration units": "mol m-3"`. Without it the reader assumes molecule cm-3 and rescales A, k0_A, and kinf_A |
| `initial` | `"initial state"`: temperature, pressure, concentrations, rate parameters. `environment` is not carried: the programs compute air density from T and p (ideal gas), which gives the canonical value |

**Check.** Passing the converted files back through the entry's own `read_configuration` gives the canonical species, third-body flags, reaction multiset, and initial values exactly.

**No driver needed.** Both programs read the CSV by name:
- a USER_DEFINED reaction's MUSICA name becomes its custom parameter's name verbatim, so the `PHOTO_*` and `EMIS_*` columns of `cells.csv` work without renaming;
- species are looked up the same way.

The Makefile always builds whatever `configs/ts1.json` holds under the C name `ts1`. Another mechanism therefore replaces TS1 (one mechanism per build), and the programs keep their `ts1_*` names and still write `ts1_final.csv`. `make bench` regenerated and recompiled only `generated/ts1.c` and `generated/ts1_specialized.c`, then relinked the library, the three programs, and `tutorials`.

**What "generic" means here.** `bin/ts1_bench_cpu` is the Tier 1 program. It calls `chem_load_ts1(1)`, which installs the configure.py-generated force and Jacobian kernels into the library's generic Rosenbrock solver and sparse LU. `bin/ts1_bench_specialized` is the fully generated solver: unrolled LU, four-cell tiles, and its own loader generated from the same file.

**Surprises:**
- Chapman is bit-identical to MICM, even in column order, because the species are reordered the same way as MICM's.
- The specialized output is bit-identical to the generic one for all three mechanisms; the generator is designed to keep the arithmetic order.

**Hazards for other mechanisms, which did not bite here:**
- `chem_model_finalize` rejects any species that no reaction uses, and the program then reports "cannot build mechanism". Every canonical species, M included, is used.
- `ts1_scenario_specialized` treats as photolysis only the custom parameters whose names begin with `j`, TS1's naming. `make bench` builds it, though this exercise does not use it, and it would scale none of the `PHOTO_*` labels.

**Rules.** `configs/ts1.json` was restored (SHA-256 unchanged, equal to git's blob), and `git status --porcelain` is empty. Nothing was written outside `$WORK` except the ignored build products of `make` inside the worktree.
