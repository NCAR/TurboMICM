# Completion audit

This audit follows PROMPT.md's requirements. Numerical results are recorded in
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
