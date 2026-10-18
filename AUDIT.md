# Retained Tier 1 completion audit

This section records the inherited Tier 1 audit; Tier 3 is audited separately below. Numerical results are recorded in
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

# Tier 3 completion audit

The following evidence was inspected after implementation. `tools/report_results.py`
rejects incomplete matrices before writing RESULTS.md; the underlying CSVs,
compiler/test logs and live command results were also inspected directly.

| Requirement | Current evidence |
|---|---|
| Own workspace and read-only inputs | All authored source, environments, builds, large run outputs and records are in this repository. Inputs used are only the supplied micm, reference and scenarios directories; their source files were not modified. |
| Preserve complete chemistry | mechanisms/ts1.json remains unchanged: 210 species, 209 integrated variables, 547 reactions. New stage loops use the same generic forcing/Jacobian/LU APIs. No species or reaction is removed, lumped or replaced. |
| Preserve Tier 1 methods and gate | Original solver and backward Euler APIs remain. Exact 680-ID manifest check; both compilers give 675 PASS, zero FAIL, the same five pure C++ N/A. Fresh 100-cell and full 10000-cell comparisons check all values at atol 1e-30 + rtol 1e-10. Clang is bit-identical; GCC max relative differences are 3.27e-14 and 1.43e-13. |
| Runtime numerical additions | include/chem_integrator.h, src/integrator.c and parameters.c: eight experimental method choices, I/PI/H211b control, error norms, positivity rejection, warm starts and exact/approximate reuse options. bin/ts1_scenario selects every new method at runtime; all have scenario scan rows. |
| Mechanism-generic testing | Independent two-cell nonlinear dimer prototype; separate exact linear/Riccati convergence and invariant checks for all eight methods. No species-name specialization in an integrator. Experimental APIs explicitly reject DAEs; the retained complete DAE solvers remain available. |
| C/prototype agreement | 48 strict comparisons per compiler have identical counts and solutions within 3e-13 (observed maxima 4.61e-15 clang, 3.77e-15 GCC). 64 cache comparisons have identical C trajectories and exactly one call saved per rejection. Aggressive b=1 controller differences are explicitly recorded/explained in prototypes/README.md and roundoff_sensitivity.csv. |
| Systematic research and sources | REFERENCES.md records every source used and what was taken. SCANS.md describes controlled changes and negative results; scans/*.csv preserves configurations, scenarios, supplied scores, operation counts and C times, including provisional/failing trials. Splitting/precision/BDF directions are identified as untested, not claimed improvements. |
| Accuracy decision | ACCURACY.md connects the literature to an explicit 0.1% overall / 1% species RMS budget, distinguishes the supplied pooled metric from the paper's worst-species metric, and reports measured margins and limitations. |
| Scenario contract | Shared CSV defaults; relative scenario paths; piecewise-linear photolysis sampled at each external-step start; all other conditions fixed; no warm-up. C timer includes rate updates and chemistry, excludes reading/setup/snapshot copies/output. Output uses 17-digit values and correct step/cell order. Independent MICM parity and partial-group/snapshot-remainder/parser tests pass with both compilers. |
| Five configurations and recommendation | OPERATING_POINTS.md and operating_points.json agree. All five points pass their declared levels on all six disclosed scenarios: 150 clang trials, plus 30 GCC checks. Three additional 48-hour scenarios add 15 explicit level checks, all PASS. |
| High-accuracy independent conditions | Three 128-cell cases, NOx x0.02 to x50, latitudes 0 and 60, two days. Tight six-stage MICM references converge by 10.076–10.349 overall digits and at least 8.992 species digits. Arguments and work/convergence records are versioned. |
| 10000-cell behavior | Nine grouping/scheduling trials, then all five declared points run diurnal-10min with 10000 cells and pass the full finite/nonnegative gate. Accuracy is honestly restricted to the disclosed first 1000 reference cells. This is distinct from the full 10000-cell Tier 1 roundoff gate. |
| Timing evidence and plots | RESULTS.md derives medians from five shuffled sequential C repetitions and compares against independently compiled MICM3 baselines (three repetitions at each qualified level). All speed data come from C/C++. Overall and species work–precision PNGs were visually inspected; matching PDF exports are included. |
| Portable implementation and flags | Both requested compilers build with make cpu, make test-cpu and make bench. ISO C11/pedantic builds, source audit, correct conditional SDK/OpenMP flags, and no extra optimization/target/FP/LTO/PGO flags. Release linkage contains only libSystem and libomp. No runtime Python, OS-specific API, intrinsic, vector extension or long double. |
| Memory/undefined behavior | Address/undefined-behavior sanitizer runs pass all eight reused-method probes, exact-order tests and scenario/partial-group/parser checks. Release clang build restored afterward. |
| Pinned local tooling | environment.yml pins Python, PyYAML, clang-format, matplotlib and numpy without build strings. make env uses the workspace prefix/cache. Python is tooling/prototype/scoring only, never invoked by the library or benchmark. |
| Required documents and commit | README.md, DESIGN.md, DECISIONS.md, RESULTS.md, ACCURACY.md, SCANS.md, REFERENCES.md, OPERATING_POINTS.md, environment.yml, prototypes/ and scans/ are present and updated. All requested work is committed; large reproducible inputs/builds stay ignored. |

Undisclosed inputs remain the experimenter's future evaluation. The implementation
uses generic adaptivity and nonfinite/negativity rejection, and the declared
points have broad current validation; no claim of a uniform mathematical error
bound on unavailable conditions is substituted for that evaluation.
