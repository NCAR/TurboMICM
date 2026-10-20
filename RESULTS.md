# Results

## MICM case list

`make test-cpu` prints one line for each of the 680 IDs in
`../reference/micm_cases.txt`, followed by aggregate and TS1 reference-check
lines. The clean results for the current source are:

| Compiler | PASS | FAIL | N/A | TS1 serial / parallel |
| --- | ---: | ---: | ---: | --- |
| Apple clang 21, macOS arm64 | 662 | 0 | 18 | PASS / PASS |
| GCC 15.2, macOS arm64 | 662 | 0 | 18 | PASS / PASS |

The 662 passing cases cover implemented rate-constant formulas and parameter
references, configuration examples, tutorials, state and object operations,
process and constraint construction and Jacobians, dense and sparse matrices,
grouped views, solver utilities, LU operations, analytical integrations,
FlowTube and Chapman regressions, DAE behavior, and OpenMP. Eighteen aerosol
cases now run against C mechanism equivalents for state layout and updates,
forcing, Jacobians, and single-cell and multi-cell solves with backward Euler
and Rosenbrock. Sixteen external-model constraint cases have focused C
regressions for algebraic metadata, process-only and constraint-only behavior,
DAE conservation, algebraic-row forcing/Jacobians, combined constraints,
multi-cell solves, equilibrium and conservation behavior, temperature-dependent
parameters, and finite-difference Jacobians.

The 18 N/A cases describe C++-specific behavior with no meaningful C counterpart:
one rate-store lambda callback entry; four compile-time view-category checks and
one grouped-view category-tag check; five version metadata checks; three
solver-builder type/matrix-policy checks; two sparse algebraic-row pattern
checks (the C DAE solver uses dense Jacobian storage and has no sparse callback
element-access path); one reusable solver-builder lifetime check (C passes
mechanism and constraint descriptors directly); and one external-model
insertion-order check (C has no model registry or order-dependent mask). Every
other case in the supplied list passes.

## TS1 comparisons

All results below use `-O3 -DNDEBUG`, double precision, and the same TS1 input
and 10-thread setting on the provided Apple arm64 machine. For the 10,000-cell
seed 2026 comparison, C used **21,199.606 ms** and MICM used **28,453.001 ms**
in an earlier local run. A fresh comparison used **21,397.132 ms** for C and
**27,850.738 ms** for MICM. Both runs reported 9,000 internal steps, with matching
accepted/rejected counts, function calls, Jacobian updates, decompositions,
and solves. All 2,110,000 final cell/species values passed
`../reference/compare_final.py`; maximum relative error was `1.77e-13`
(BCARY, cell 6848). The task's hidden seed and five-run median were not measured
here, so these timings are indicative local comparisons, not the ranked result.

The fixed 100-cell, 100-step, 30-second seed 1 input passes for one and four
threads. A 1,000-cell seed 2026 run passed all 211,000 values with maximum
relative error `7.74e-14`. Earlier exploratory 100-cell seed 1 comparison had
maximum relative error `9.7e-13`.

## Build verification

Clean `make test-cpu` builds completed with Apple clang 21 and GCC 15.2. Both
reported 662 PASS, 0 FAIL, and 18 N/A, and passed the serial and four-thread
TS1 checks. `make cpu` and `make bench` also built the benchmark with both
compilers. The linked benchmark depends on the C runtime, libm, and OpenMP only;
it has no Python linkage. All cases in the supplied list are accounted for as
PASS or an API-specific N/A.
