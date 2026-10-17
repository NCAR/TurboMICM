# Porting MICM's test cases

Every case in `../reference/micm_cases.txt` (680 lines, `<id> gtest|program`)
becomes one C function here, registered under MICM's id. `make test-cpu`
prints one line per case in the list's order: `<id> PASS|FAIL|N/A <reason>`.
A case with no port prints `FAIL not ported`.

## Layout

- `harness.h`: the `TMC_TEST`, `TMC_CASE`, `TMC_EXPECT_*`, `TMC_ASSERT_TRUE`,
  and `TMC_NA` macros, and the path helpers `tmc_test_micm_path`,
  `tmc_test_reference_path`, `tmc_test_data_path`.
- `main.c`: runs the cases; `tools/test_registry.py` writes
  `generated/tests/registry.h` (the groups) and `case_list.inc` (the ids) at
  build time from the files here.
- `test_<ctest name>.c`: one file per MICM test executable (or a related set),
  exporting `const tmc_test_case *tmc_tests_<name>(int *n)` whose array ends
  with `{ NULL, NULL }`. The registry finds that signature by text, so keep it
  on one line, exactly `const tmc_test_case *tmc_tests_<name>(int *n)`.
- `data/`: expected outputs and inputs copied or converted from MICM's tests.

## Rules

- Port every case: check the same quantities with the same or tighter
  tolerances. Where MICM throws, the C API returns an error category and code
  (`tmc_error`); check those. A case that tests C++ mechanics with no C
  equivalent (template deduction, SFINAE, overload resolution) is `TMC_NA`
  with a one-line reason; a case about behavior (copying, errors, layouts,
  iteration order) is ported to the equivalent C behavior.
- MICM's test sources are read-only, under `tmc_test_micm_path("test/...")`
  (`../micm` by default; the Makefile sets `TMC_TEST_MICM`).
- Numerical comparisons: `TMC_EXPECT_DOUBLE_EQ` is GoogleTest's 4-ULP
  comparison; `TMC_EXPECT_NEAR` is absolute; `TMC_EXPECT_REL` is relative.
- Randomized MICM tests (random sparse matrices, lognormal values) use fixed
  seeds in C (`srand(1)` or a small LCG in the test); the comparison is still
  against an independent computation inside the test, so the seed does not
  change what is verified.
- Layouts: MICM's `Matrix` is `tmc_dense` with `L = 1`; `VectorMatrix<L>` is
  `tmc_dense` with that `L`. `SparseMatrix<T, StandardOrdering>` is
  `tmc_sparse` with `L = 1`; `VectorOrdering<L>` is `L`. Compressed sparse row
  and column are `TMC_SPARSE_CSR` and `TMC_SPARSE_CSC`. A MICM test that is
  instantiated for several policies (for example `Matrix` and
  `VectorMatrix<1..4>`) becomes one C case that runs the same check for each
  layout, since the case id is the same.
- Library gaps: a test may need a function the library lacks. Add it to
  `include/tmc.h` and the matching `src/*.c` file when it is a self-contained
  accessor or utility (a getter, a copy, an element-wise operation) that does
  not change existing behavior, keep the style of the surrounding code, and
  list the addition in your report. Anything touching the solver's numerics
  or the builder is reported, not changed.
- Build and run one group while porting:
  `make test-cpu TEST_FILTER=<ctest name>` (the filter is a substring of the id).
  The library's own build is `make cpu`.
