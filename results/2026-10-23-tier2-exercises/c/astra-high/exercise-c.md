# Exercise C, astra-high: the generator on other mechanisms

Entry: `tier2/astra/high` (f478080), detached worktree `astra-high/`. Route: **benchmark programs**.
The build reads `mechanisms/ts1.json`. The Makefile hard-codes that path and `--symbol ts1`, so each
converted mechanism was swapped in temporarily. The drivers (`src/ts1_bench.c`, `src/ts1_driver.h`)
resolve CSV columns and output species by name, so they run any mechanism unchanged. Converter: 76
lines. Driver: none. Entry changes needed: none. At the end `mechanisms/ts1.json` is byte-identical
to HEAD, TS1 is rebuilt and passes its sanity check (difference 0 against MICM), and
`git status --porcelain` is empty.

## Results: 4 identical cells (cells.csv), 1 warm-up + 60 steps of 60 s, default ROS3 settings

Clang, `OPTFLAGS="-O3 -DNDEBUG"`, one thread, default `make bench` (GROUP_SIZE 32).
"PASS 0" means `compare_final.py` passed with a largest relative difference of 0: every value is
bit-identical.

| Mechanism | Converted | Generator (s) | Generated C (lines / bytes) | make bench (exit, s) | Generic ran | Specialized ran | Spec. vs generic | Generic vs MICM | Spec. vs MICM | Generated kernels used |
|---|---|---|---|---|---|---|---|---|---|---|
| chapman | yes | 0.030 | 569 / 20,202 | 0, 1.01 | yes | yes | PASS 0 | PASS 0 | PASS 0 | **no** (generic fallback) |
| carbon_bond_5 | yes | 0.050 | 15,649 / 682,538 | 0, 7.03 | yes | yes | PASS 0 | PASS 0 | PASS 0 | **no** (generic fallback) |
| carbon_bond_5_no_emissions | yes | 0.052 | 15,303 / 669,240 | 0, 5.62 | yes | yes | PASS 0 | PASS 0 | PASS 0 | **no** (generic fallback) |
| *TS1, for scale* | | 0.660 | 44,974 / 2,048,850 | 0, 23.4 | | | | | | |

- Generator time: `tools/specialize.py`, run as make runs it, median of 7 runs, including Python
  start-up. The generic builder (`tools/configure.py`) takes 0.03 s. The make-generated
  `build/specialized/ts1.c` is byte-identical to a standalone run.
- make bench: median of 3 full rebuilds after the swap (chapman 1.02/1.01/0.99, CB5 7.03/6.98/7.35,
  CB5 without emissions 6.09/5.62/5.41 s). Each rebuild regenerates and compiles the builder and the
  specialized C and links the three drivers; `libchem.a` is unchanged. Other exercises' compiles
  were running on the machine (load about 3 on 18 cores).
- No compiler warnings with `-Wall -Wextra -Wpedantic`.
- Emissions do enter the result. CB5 and CB5 without emissions differ in 244 of 272 values (CO
  1.3e-5 against 6.0e-6 mol m-3), and the entry matches MICM exactly for both.

### Same 4-cell problem, rebuilt with `make bench ... GROUP_SIZE=4`

This uses a documented Makefile variable, not a code change. The generated stride becomes 4, equal
to the 4-cell group, so the generated kernels run.

| Mechanism | Spec. vs generic | Spec. vs MICM | elapsed_ms: generic / specialized, default / specialized, GROUP_SIZE=4 |
|---|---|---|---|
| chapman | PASS 0 | PASS 0 | 0.23 / 0.22 / 0.07 |
| carbon_bond_5 | PASS 0 | PASS 0 | 10.8 / 10.6 / 2.1 |
| carbon_bond_5_no_emissions | PASS 0 | PASS 0 | 10.6 / 10.0 / 2.2 |

Each figure is a single run, so the timings are indicative only.

## Addendum: 64 distinct cells (cells64.csv), same default build

| Mechanism | Spec. vs generic | Generic vs MICM (micm-run64) | Spec. vs MICM | Generated kernels used | elapsed_ms, median of 5: generic / specialized |
|---|---|---|---|---|---|
| chapman | PASS 0 | PASS 0 | PASS 0 | **yes** | 1.25 / 1.54 |
| carbon_bond_5 | PASS 0 | PASS 0 | PASS 0 | **yes** | 59.0 / 48.0 |
| carbon_bond_5_no_emissions | PASS 0 | PASS 0 | PASS 0 | **yes** | 58.7 / 48.0 |

`micm-run64-g32/ts1_final.csv` is the same file as `micm-run64/ts1_final.csv`, and the specialized
results also pass against it with difference 0. The generated C was identical to the main runs'.

### Did the specialized program run the generated kernels?

**Code path.** For the specialized binary, `src/ts1_bench.c:16-20` sets
`group_size = CHEM_GROUP_SIZE` (`-DCHEM_GROUP_SIZE=32`) and then
`if (group_size>cells) group_size=cells;`.

`src/ts1_driver.h` maps `chem_update_rates` and `chem_solve` to `ts1_update_rates` and `ts1_solve`
(lines 19-20). It also calls `ts1_specialized_attach` (line 35), which installs the generated
`factor` and `solve` as `s->linear`; `src/rosenbrock_step.inc:44,63` calls them with `s->cells`.

Each generated kernel then chooses by group size. Line numbers are from `runs/carbon_bond_5/specialized.c`:

- `ts1_update_rates`: `if (s->cells!=32) chem_update_rates(s); else rates(s);` (4205)
- `forcing`: `if (s->cells!=32) { chem_forcing(s,y,f); return; }` (6859)
- Jacobian macro: `((s)->cells==32 ? jacobian(...) : chem_solver_negative_jacobian(...))` (15641)
- `factor`: `if (n!=32) return chem_lu_factor(a,n,v,CHEM_LU_MOZART);` (13939)
- `solve`: `if (n!=32) { chem_lu_solve(a,n,v,x); return; }` (15361)

With 4 cells there is one 4-cell group, so every kernel takes the generic branch. With 64 cells
there are two 32-cell groups, so every kernel takes the generated branch. The README documents the
fallback. In the 4-cell run the generated integrator instantiation (`ts1_solve3`), the
variable-stride `diagonal` and the rate cache still run; only the mechanism-specific kernels fall
back.

**Run-time confirmation.** I built an instrumented copy in `coverage/`, outside the entry, from the
entry's unchanged sources and the same generated C, adding `-fprofile-instr-generate
-fcoverage-mapping`. Its outputs are identical to the real runs.

| Mechanism | 4 cells: generated kernel functions entered | 4 cells: generic kernels entered | 64 cells: generated entered | 64 cells: generic entered |
|---|---|---|---|---|
| chapman | 0 / 5 | all 5 | 5 / 5 | none |
| carbon_bond_5 | 0 / 32 | all 5 | 32 / 32 | none |
| carbon_bond_5_no_emissions | 0 / 30 | all 5 | 30 / 30 | none |

The generated kernel functions are `rates_*`, `forcing_add_*`, `jacobian_add_*`, `factor_*` and
`triangular_*`. The generic kernels are `chem_update_rates`, `chem_forcing`,
`chem_solver_negative_jacobian`, `chem_lu_factor` and `chem_lu_solve`.

In these OpenMP builds clang also charges every `omp simd` loop entry and iteration to a function's
entry counter. For example, `diagonal` shows 549 × 331 = 181,719. The counts are therefore inflated,
but 0 still means the function never ran.

Timing is weak evidence. Setting the driver's `CHEM_GROUP_SIZE=16` forces the fallback in the same
binary: CB5 then takes 68.4 ms generic and 63.6 ms specialized, against 59.0 and 48.0 ms with groups
of 32. The specialized integrator and rate cache stay in use after the fallback, and for chapman the
generated kernels at stride 32 are slower than the generic ones.

## Commands

```
WORK=~/NSF-NCAR/tier2-work/exercises/c/astra-high; mkdir -p $WORK && cd $WORK
ln -sn .../ODE/benchmarks/micm micm; ln -sn .../reference reference; ln -sn .../tier2/scenarios scenarios
git -C ~/NSF-NCAR/TurboMICM worktree add --detach $WORK/astra-high tier2/astra/high
# make env would download from conda-forge into a repo-local cache, so it ran offline from a copy of
# the package cache this entry's earlier exercise (ab/astra-high) had already downloaded:
cp -Rp ../ab/astra-high/astra-high/.conda_pkgs astra-high/.conda_pkgs
cd astra-high && CONDA_OFFLINE=true make env          # exit 0, 8 s, nothing downloaded
make cpu bench CC=clang OPTFLAGS="-O3 -DNDEBUG"        # exit 0, 24.4 s; TS1 sanity: difference 0
python3 $WORK/convert.py canonical/MECH/canonical.json $WORK/configs/MECH/ts1.json
$WORK/run_exercise.sh   # backup ts1.json; per MECH: cp config -> mechanisms/ts1.json;
                        # make bench CC=clang OPTFLAGS="-O3 -DNDEBUG"; in empty dirs:
                        # bin/ts1_bench_cpu --threads 1 4 60 60 cells.csv (and ts1_bench_specialized);
                        # compare_final.py; then make bench ... GROUP_SIZE=4 and rerun specialized; restore
python3 $WORK/measure.py [MECH]   # generator x7 and make bench x3 per MECH (1.1 s pause per swap)
$WORK/run_cells64.sh    # per MECH: swap, default make bench, both programs on 64 cells, compare
                        # with each other and micm-run64(-g32); CHEM_GROUP_SIZE=16 timing; restore, rebuild
coverage/: instrumented copy (clang -fprofile-instr-generate -fcoverage-mapping), llvm-profdata, summarize.py
```

Logs and outputs: `runs/MECH/` (4 cells), `runs64/MECH/` (64 cells), `measure/` (timings), `gen/`
(standalone generator outputs), `coverage/` (run-time confirmation).

## What it took

- **Converter (76 lines),** mapping each feature as `mechanisms/ts1.json` does:
  - Species become `CHEM_SPEC` records.
  - The third body M becomes `"parameterization": "air density"`. It stays in the reactant lists,
    where the reader multiplies the rate by air density once per occurrence, and in the product
    lists, where the reader drops it, as TS1 does.
  - One `MECHANISM` record carries `"rate units": "mol m-3"`. This line is essential: the reader's
    default is molecule cm-3, which would silently rescale A, k0_A and kinf_A.
  - Reactants are ordered lists with repeats. Products are `{species, yield}` lists.
  - ARRHENIUS keeps A–E and TROE keeps k0_*, kinf_*, Fc and N.
  - USER_DEFINED becomes `"scaling factor"` plus `"MUSICA name"`, the label, which is also the
    cells.csv column. Emissions are USER_DEFINED reactions with an empty reactant list, and the
    reader and generator accept them as is.
  - `"initial state"` holds the given concentrations (never M, which the reader rejects), then the
    rate parameters, then temperature and pressure. Air density is omitted because the entry
    computes it from the ideal gas law, as MICM's driver does.
  - Molecular weights are not needed, since there are no SURFACE reactions.
  - No file was hand-edited.
- **Generator:** accepted all three mechanisms on the first try. It took 0.03–0.05 s and produced
  15.6k lines of C for CB5's 200 reactions and 66 species (TS1: 45k for 547 reactions and 209
  species). The C compiled without warnings and was bit-identical to the entry's generic path and
  to MICM on 4 and 64 cells. The entry reproduces MICM's operation order exactly
  (rate-type sort, Markowitz tie-breaking, Mozart LU and accumulation order), as it does for TS1.
- **Surprises:**
  - With the exercise's 4 cells, the specialized program does not run its generated kernels (see
    above), so the 4-cell "specialized vs generic" comparison mostly compares the generic kernels
    with themselves. The 64-cell runs and the GROUP_SIZE=4 build close that gap.
  - `/usr/bin/make` is GNU Make 3.81, which compares timestamps in whole seconds. A swap that lands
    in the same second as the previous build is not seen as newer. Chapman builds in about 1 s, and
    two of my first timed rebuilds were partial, so I discarded and redid them with a pause.
  - `make env` needs network for a fresh worktree, since it pins Python 3.12.12 with a repo-local
    package cache.
