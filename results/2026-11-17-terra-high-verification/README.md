# Verification of terra-high, 2026-11-17

terra-high, the third Terra entry (gpt-6.1-terra at high effort, Codex CLI
0.169.0), was rebuilt from its branch, retested, and audited with both
compilers on 2026-11-17 at 08:46, eleven minutes after the run ended. It was
timed with every other entry that evening, in one session (see [Timing](#timing)).

| File | Contents |
|---|---|
| `tests.md` | `make test-cpu` results for terra-high with each compiler |

## Builds, tests, and audit

terra-high (`a8a2ea0`, 79 commits) was built with `make env`, `make cpu bench
OPTFLAGS="-O3 -DNDEBUG"`, and `make test-cpu`, once with `CC=clang` (Apple
clang 21.0.0) and once with `CC=gcc-15` (GCC 15.2.0), by
`admin/verify_tier3.py --stages build,tests,ts1,audit`, the method of
[`2026-11-17-terra-max-verification`](../2026-11-17-terra-max-verification/README.md).

- **Tests** (`tests.md`): 675 PASS, 0 FAIL, 5 N/A with each compiler, every id
  exactly once, as the agent reported; the N/A cases are terra-max's four C++
  template-concept cases and `GroupedView.CategoryTagsAreDistinct`, a check of
  C++ tag types.
- **How the PASS lines are made.** Unlike terra-max and terra-xhigh, whose test
  programs print one result per case, `tools/test_cpu.py` maps each case id to
  one of 47 test programs (`tools/test_cpu.py:9-293`) and reports the id PASS
  when that program exits 0 (`tools/test_cpu.py:316-324`): the 675 PASS ids
  share 47 exit statuses, 579 of them in programs that carry ten or more ids
  (86 from `build/test_sparse_functions`, 42 from `build/test_sparse_views`).
  The file states that "each listed executable must check every assertion in
  each mapped case" (`tools/test_cpu.py:8`). A partial, read-only audit on
  2026-11-17: the asserts are live in the `-O3 -DNDEBUG` build, since every
  test source undefines `NDEBUG` before `<assert.h>`
  (`tests/test_doolittle.c:3-6`). But `build/test_doolittle` and
  `build/test_linear_solve` factor and solve one cell's matrix in standard
  storage only (`tests/test_doolittle.c:13-50`; `tests/test_linear_solve.c:13-42`),
  and carry 24 of MICM's vector-ordering ids
  (`lu_decomposition_{doolittle,mozart}[_in_place]:…VectorOrdering`,
  `linear_solver[_in_place]:…VectorOrdering`): those 24 PASS lines rest on no
  check of a vector-ordered layout. terra-max and terra-xhigh test these
  cases at grouped widths. Three more programs carry 5 vector-ordered ids
  (`ConstraintSet.VectorMatrix1-3`, `ProcessSet.VectorMatrix`,
  `Process.VectorMatrix`) with no grouped-layout code in their sources (not
  read further); the other 157 vector-ordered ids come from programs that do
  build grouped layouts, and the rest of the mapping was not traced. The
  judges recorded the same reporting pattern for luna-max and luna-xhigh.
- **TS1 gate:** passes with both compilers; largest relative difference from
  MICM 1.25e-13 with clang (BCARY, cell 6848) and 1.43e-13 with gcc-15
  (BCARYO2VBS, cell 4054).
- **Flags, libraries, sources:** only `-std=c11`, warnings, `-O3 -DNDEBUG`,
  include paths, and OpenMP flags. The benchmark, the one program in `bin/`,
  links only the OpenMP runtime and `libSystem`; the 56 test and tutorial
  programs in `build/`, checked by hand, link only `libSystem` (three also the
  OpenMP runtime). The scan of 114 C sources for banned constructs matched
  nothing. The TS1 generator runs with `--reference-order` (`Makefile:4`),
  which sorts a configuration's reactions by rate type and applies MICM's
  Markowitz ordering, unsigned-cost tie-break included
  (`tools/generate.py:41-63,66-77`); it applies to any mechanism.
- **Warnings:** none, with either compiler.
- **Same work, same clock:** `bench/ts1_bench.c` solves each cell on its own,
  with its own step-size control, an `omp parallel for` over the 10,000 cells,
  static schedule unless `OMP_SCHEDULE` is set (`bench/ts1_bench.c:173-189`).
  It computes the rate constants in one untimed warm-up interval per cell
  (`:174-180`) and times the 100 intervals with `omp_get_wtime` (`:181-190`).
  The solver allocates its working memory on each call, inside the timer:
  `chem_rosenbrock_solve` `calloc`s n² + 10n doubles (about 366 KB for TS1's
  209 state species) and frees them before returning (`src/rosenbrock.c:167-175`),
  once per cell and 30 s step. Its solver work at ten threads: 9,000,000
  internal steps, all accepted, 18,000,000 function calls, 9,000,000 Jacobian
  updates and decompositions, and 27,000,000 solves, MICM's 9 steps per
  interval in every cell.
- **Self-report** (one line, not a timing): its RESULTS.md gives 25,689 ms
  against 26,324 ms for MICM's default, one run each at 10 threads.

## The case mapping, continued (2026-11-17, evening)

A second read-only pass covered the 19 programs whose ids name a layout, an
ordering, or a solver variant; they carry 470 of the 675 PASS ids. For each, the
question was whether the program runs its checks in every variant its ids name.

- **Layouts and orderings.** Eleven programs run their checks in every layout and
  ordering their ids name: `test_matrix`, `test_matrix_function`, and
  `test_matrix_views` at group widths 0 (row-major) to 4; `test_sparse_matrix`,
  `test_sparse_functions`, `test_sparse_views`, and `test_typed_views` at widths
  0–4 in CSR and CSC order; `test_grouped_views` at widths 0, 1, 2, 4, and 8 in
  both orders; `test_alpha_jac` at widths 0–4, both orders, and 1–4 cells;
  `test_reorder` at widths 0–4; and `test_chapman_regression`, whose grouped
  solve runs at widths 0–4 with the five parameter sets
  (`tests/test_chapman_regression.c:203-220`). Each runs one sequence of checks
  per variant, not one check per id: the 86 ids of `test_sparse_functions` rest
  on one 80-line sequence run in ten layouts (`tests/test_sparse_functions.c:35-114`).
- **Solver variants.** `test_chain_rates` and `test_nonlinear` run the five
  Rosenbrock parameter sets and backward Euler where their ids name them;
  `test_doolittle` and `test_linear_solve` run Mozart and Doolittle, sparse and
  dense, in place and copied, from five initial values
  (`tests/test_doolittle.c:13-50`).
- **36 vector-ordered ids rest on standard storage.** Besides the 24 LU and
  linear-solver ids above: `RosenbrockSolver.VectorNormalizedError` and
  `VectorNormalizedErrorWithConstraints`, checked on cell-major arrays of 1–5
  cells (`tests/test_normalized_error.c:12-39`);
  `RegressionRosenbrock.VectorRateConstants`, `VectorForcing`, and
  `VectorJacobian`, and `RosenbrockSolver.VectorTerminator`, checked one cell at
  a time (`tests/test_nonlinear.c:224-249,271-337`); `ConstraintSet.VectorMatrix1`
  to `3` and `ConstraintSet.VectorizedMatricesRespectGridCellIndexing`, on
  cell-major arrays (`tests/test_constraint_set.c:31-78`); and
  `ProcessSet.VectorMatrix` and `Process.VectorMatrix`, in programs with no
  grouped-layout code. MICM's versions of these cases run on vector-ordered
  matrices.
- **Two constraint ids.** The two `NormalizedErrorWithConstraints` ids rest on a
  block with no algebraic rows that computes one norm five times
  (`tests/test_normalized_error.c:26-39`); MICM's cases check that the norm
  includes every variable under each of the five parameter sets (MICM
  `test/unit/solver/test_rosenbrock.cpp:213-243`).
- **Not traced.** The other 28 programs carry 205 ids that name no layout or
  solver variant; whether each of their checks covers every assertion of each
  mapped case was not compared with MICM's tests.

For the SSE, with the Luna gate findings: terra-high's 675 PASS lines are 47
programs' exit statuses, and 36 of them name a vector-ordered layout that no
check of theirs builds.

## Timing

Timed on 2026-11-17, 20:49 to 23:39, with every other entry and MICM's three
programs, in
[`2026-11-17-tier1-all-twelve-timing`](../2026-11-17-tier1-all-twelve-timing/README.md):
5 rounds at 10 and 1 threads, every run passing the gate.

| | 10 threads, median ms (range) | vs MICM's default | vs MICM's fastest | 1 thread, median ms (range) | vs MICM's default | vs MICM's fastest |
|---|---:|---:|---:|---:|---:|---:|
| terra-high | 25,343 (25,173–26,021) | 1.08× | 0.38× | 144,516 (143,201–144,696) | 1.01× | 0.35× |

It is ninth of the twelve entries at both thread counts, 1.15× faster than
terra-xhigh at ten threads and 1.17× on one, and gains 5.70× from one thread to
ten; its times include the work arrays it allocates and frees on every call. Its
self-report, 25,689 ms against 26,324 ms for MICM's default in one ten-thread
run, compares with 25,343 and 27,265 ms here.

## Compliance

`admin/scan_commands.py` on the session record: 952 commands. Two path
findings, both benign: a mistyped `../micmicm` in a command that otherwise
reads `../micm/` (13:30 UTC) and `/Users` as a search pattern for `rg` inside
the run's own files (15:11 UTC). None naming another run or the experiment's
private folders, no installs or downloads, no connector, MCP, or web-search
calls, and no subagents. All 110 file edits went to 71 files inside the run
directory.
