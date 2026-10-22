# Completion audit

The Tier 2 requirements were checked against the final files, actual build/test
logs and numerical output. RESULTS.md records the measurements; ignored local
artifacts remain under `results/`. The original Tier 1 audit is retained below.

| Tier 2 requirement | Current evidence |
|---|---|
| Work and commits stay in this repository; reference trees read-only | Git top level equals this workspace. Environments, caches, builds, reference executables and numerical outputs are local. |
| Python generator reads the Tier 1 MICM/CAMP format | `tools/specialize.py` shares `configure.read`; no mechanism-name branch in binding or emission. |
| Any supported mechanism, including TS1 and the Tier 1 configuration files | `make test-specialized`: TS1, three original configurations, all-rate and zero-order fixtures; 70 numerical cases pass with both compilers. |
| Reaction/parameter edits regenerate without generator edits | All schedules derive from parsed data; structural edit test changes a TS1 reaction and its parameter. Tests derive dimensions from the input rather than freezing TS1 counts. |
| All runtime quantities remain data | Generated rates read conditions/custom arrays; kernels read concentrations; the solve reads dt and per-state tolerances. Input changes and day-night transitions are exercised numerically. |
| Preserve MICM's three-stage method and controller | `src/specialized_solver.h`; numerical kernel/trajectory/statistics comparisons and stringent day-night agreement. No changed error norm or factorization reuse. |
| Generated C is a build product | Ignored `generated/`; Make dependencies include configurations, generator, template and tuning stamp; both clean builds regenerate it. |
| Keep Tier 1 library, CPU benchmark and gate | Exact 680-case list equals the reference; both compilers report 675 PASS / 0 FAIL / 5 original C++-only N/A. Generic TS1 passes the small and full-size MICM comparisons. |
| Pinned local conda environment; build-time Python only | `environment.yml`, local `make env` creation, `CONDA_PKGS_DIRS` inside workspace. Both benchmark dependency inspections show only the system C runtime and OpenMP. |
| ISO C11, doubles, portable SIMD and permitted libraries | Both builds use C11/pedantic flags without diagnostics. Source audit finds no forbidden vector extensions, intrinsics, attributes, builtins, assembly, long double or OS APIs. |
| CC/OPTFLAGS and compiler-specific OpenMP | Actual clang 21 and GCC 15.3 clean builds; direct GCC-to-clang switch also passes both suites. Make adds only language/warning/include/OpenMP/link flags. |
| Tunable defaults, correctness at other settings | DESIGN.md lists all controls. Alternate generated fixture uses fusion, row updates, odd stride, uneven tiles and different batches/split; numerical checks pass on both compilers. |
| Specialized benchmark contract | Shared Tier 1 CLI/CSV/group/output code, generated kernels, one untimed warm-up, timed steps, elapsed_ms/statistics; 100 and 10,000-cell MICM gates pass. |
| Scenario runtime contract | `src/scenario.c`: relative paths, original custom values, piecewise-linear factors, every j-label updated per step, no chemistry warm-up, correct snapshots/partial groups/tolerance overrides. |
| Every disclosed scenario passes with both compilers | Ten final-state and ten snapshot gates pass with `score.py --level 10`, zero nonfinite/below-floor values; roughly 14 overall digits and at least 13.424 for the worst species. |
| Evidence beyond disclosed scenarios | Independent MICM runs for 48-hour northern-winter and 40-hour southern-winter inputs; both compilers pass, including final snapshots off the interval and partial groups. |
| Performance, build/source sizes, attempted specializations | RESULTS.md gives 18/1-thread C/MICM times, both clean build times, generated C size and day-night scores; SPECIALIZATION.md includes retained and dropped experiments. |
| Required documents | SPECIALIZATION.md, README.md, DESIGN.md, DECISIONS.md, RESULTS.md and environment.yml are present and updated. |

## Retained Tier 1 audit

This section records the original Tier 1 requirements. Numerical results are recorded in
RESULTS.md; implementation and extension details are in DESIGN.md.

| Requirement | Evidence |
|---|---|
| Own repository and workspace; reference inputs read-only | All source, generated products, environments, caches, copied reference executables and run outputs are inside this repository. |
| Replace CPU functionality with C | `include/chem*.h`, `src/*.c`; standard/vector dense and sparse operations, rate types, state/system operations, both solvers, constraints, external models and utilities. |
| Every reference case, exact IDs and bounds | `tests/cases.txt` equals the supplied 680-case list; runner covers 675 PASS and five C++-only N/A. Source assertions, numerical vectors, formulas and tolerances are retained in C tests. |
| Plain-program output | Eight tutorial fixture comparisons cover deterministic printed states/statistics. Timing and interleaved thread-start logging are excluded as documented nondeterminism. |
| C++ mechanics get C equivalents | Error statuses replace exceptions; descriptor assignment and clone/move replace object mechanics; runtime layout/shape/writability metadata replaces policies. Only pure compile-time concepts/categories are N/A. |
| Generic mechanisms in C and YAML/JSON | Common descriptor/binding path; three example configurations compare against C-built mechanisms for rates, forcing, Jacobians and integration. |
| Own configuration reader and retained TS1 translation | `tools/configure.py`, `tools/import_ts1.py`, committed `mechanisms/ts1.json`; regenerated JSON compared byte-for-byte. Ten reader tests cover invalid inputs and extra rate schemas. |
| No runtime Python or MICM dependency | Static C library and C benchmark; linked-library inspection; reference C++ appears only in `micm-ref` and reference runs. |
| Preserve numerical methods | All five Rosenbrock tables, RMS control, algebraic step-change error, initialization/backtracking, backward Euler Newton/reductions, Mozart and Doolittle LU are tested. |
| Extensible solvers | Parameter tables, controller callbacks, linear-operation callbacks and explicit model/state interfaces; extension guide in DESIGN.md. |
| ISO C11 / permitted libraries and SIMD | Builds with C11 and pedantic diagnostics, standard C/libm/OpenMP only; contiguous cell loops use restrict/omp simd; source audit finds no target-specific code or long double. |
| No tuned compiler flags | Make passes `CC` and `OPTFLAGS` (default `-O3 -DNDEBUG`) and adds only language, warnings, includes, OpenMP and linking. Sanitizers were a separate diagnostic build. |
| Cross-compiler builds | `make test-cpu cpu` passes with clang and gcc-15; configuration stamp rebuilds C on compiler/flag changes. macOS SDK/OpenMP flags selected conditionally. |
| Local, reproducible Python tooling | Pinned environment.yml; `.conda` and `.conda-pkgs`, ignored by Git. Every Python run uses the prefix interpreter. |
| Parameterized portable performance choices | Group size, thread count and OpenMP schedule configurable; explicit defaults and valid ranges in DESIGN.md. |
| Benchmark contract | CSV names/defaults, actual final partial group, untimed warm-up, timed steps, summed statistics, final named CSV; default all OpenMP threads and optional thread limit. |
| TS1 gate in every cell/species | Supplied fixture and full-size 10,000-cell comparison use atol 1e-30 + rtol 1e-10. Both compilers pass; several group sizes/schedules and another input seed checked. |
| Required documents and commit | README.md, DESIGN.md, DECISIONS.md, RESULTS.md and environment.yml; completed work committed. |

The five exclusions are `ViewCategorySFINAE.VectorTypes`, `.MatrixViews`,
`.SparseMatrixViews`, `.MatrixTypesExcluded`, and
`GroupedView.CategoryTagsAreDistinct`. Their source bodies assert only C++ type
concepts and overload/category identity. The unified C view interface has no
such type dispatch; all corresponding runtime view operations are tested.
