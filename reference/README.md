# Reference inputs and tools

Read-only inputs for the MICM-to-C runs. ROOT is the directory above this one.
MICM is ROOT/micm at commit e15898ce (v3.7.0-259), built in ROOT/micm/build with
Apple clang 21: Release, double precision, OpenMP on. All 67 ctest tests pass.

## MICM's test cases

`micm_cases.txt` lists the 680 cases of MICM's test suite: 672 GoogleTest cases in
59 executables, and 8 plain programs (MICM's tutorials). Each line is
`<id> <kind>`:

    <ctest name>:<Suite.Test> gtest      a GoogleTest case
    <ctest name> program                 a plain program

The executable for ctest name `X` is `ROOT/micm/build/test_X`. MICM runs every
test from `ROOT/micm/build`, and some tests read or write files there, so copy
the build directory into your own directory and run the executables from the
copy, for example:

    cp -R ROOT/micm/build my_micm_build
    cd my_micm_build && ./test_rosenbrock --gtest_filter='RosenbrockSolver.*'

(`ctest` in the copy would still run the originals in `ROOT/micm/build`.)
`ROOT/admin/make_case_list.py` wrote the list from the build.

## TS1

TS1 is `ROOT/micm/test/benchmark/ts1_mechanism.hpp`: 210 species (209 in the
solver state) and 547 reactions.

`micm_ts1_ref.cpp` is MICM's run of the TS1 benchmark, the reference for the gate
and for timing. Build it with clang++ and your OPTFLAGS; on this machine clang
links only with SDKROOT set:

    SDKROOT=$(xcrun --sdk macosx --show-sdk-path) clang++ -std=c++20 -O3 -DNDEBUG \
        -I ROOT/micm/include -I ROOT/micm/test/benchmark \
        ROOT/reference/micm_ts1_ref.cpp -o bin/micm_ts1_ref

    bin/micm_ts1_ref [--threads N] [--group-size G] CELLS STEPS DT CELLS_CSV
    bin/micm_ts1_ref --write-defaults FILE

- Solver: MICM's three-stage Rosenbrock with its default parameters and
  tolerances, in micm_bench's default configuration (standard ordering, in-place
  Mozart LU).
- Each cell starts from TS1's `InitState`, then takes the values in its
  CELLS_CSV row (temperature, pressure, species, and rate parameters); air
  density is then recomputed from temperature and pressure. CELLS uses the first
  CELLS rows.
- One untimed warm-up step, then STEPS timed steps of DT seconds. Prints
  `elapsed_ms=` and the solver statistics summed over the timed steps, and writes
  the final state to `ts1_final.csv` in the working directory.
- Cells are solved in groups, each group a MICM state with its own step-size
  control, and the statistics are summed over the groups. By default there is
  one group per thread: with N threads (default: all cores) the cells are split
  into N contiguous groups. With `--group-size G` (added 2026-10-15) the cells are
  split into groups of G, and each thread takes a contiguous block of groups, as
  an OpenMP static schedule does, solving each group's STEPS steps in turn with
  one solver per thread. On the 100-cell input below, 1 and 4 threads, and groups
  of 7 or 32 cells, all give identical final states.

`ts1_defaults.csv` is TS1's `InitState` for one cell (`--write-defaults`).

`make_cells_csv.py CELLS SEED` writes a CELLS_CSV in which every cell differs;
its docstring gives the perturbations. It needs only the Python standard library.

`ts1_cells_100_seed1.csv` is `make_cells_csv.py 100 1`, and
`micm_ts1_final_100_seed1.csv` is MICM's final state for it:

    bin/micm_ts1_ref --threads 1 100 100 30 ts1_cells_100_seed1.csv

`compare_final.py YOURS MICM [--rtol 1e-10] [--atol 1e-30]` applies the gate to
two `ts1_final.csv` files and exits 1 on failure.

The ranking run is 10,000 cells, 100 steps, DT = 30 s, with a CELLS_CSV from
`make_cells_csv.py` and an undisclosed seed. With seed 2026, MICM built as above
takes 137,383 ms on one thread of this machine (median of 3 runs: 136,052,
137,383, and 138,060 ms; 9 internal steps per 30 s step, none rejected).
