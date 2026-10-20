# Results

## MICM case list

`make test-cpu` prints one line for each of the 680 IDs in
`../reference/micm_cases.txt`, followed by aggregate and TS1 reference-check
lines. Results for the current source are:

| Compiler | PASS | FAIL | N/A | TS1 serial / parallel |
| --- | ---: | ---: | ---: | --- |
| Apple clang 21, macOS arm64 | 55 | 613 | 12 | PASS / PASS |
| GCC 15.2, macOS arm64 | 55 | 613 | 12 | PASS / PASS |

The 55 passing case lines cover implemented rate-constant formulas, rate
parameter descriptor references, the two configuration-versus-hand-built
examples, the README chemistry example, multiple grid cells, the OpenMP tutorial,
the vectorized storage example, the solve_results tutorial, seven Jacobian
verification cases, and the OpenMP TS1 path. The 12 N/A cases cover C++ callable
storage, template type traits, and release metadata with no meaningful
equivalent in the current C API. The remaining 613 cases are explicitly
reported as FAIL; they cover
unimplemented MICM builder/state APIs, generic matrix and view utilities,
standalone linear algebra cases, constraint/DAE and alternate integration
paths, external/aerosol model APIs, and a small output mismatch in the four-
and six-stage DA tutorial at 1400 s. The five parameter sets match MICM's step
statistics; the two DA trajectories
differ in one printed species value. The per-case reasons are printed by the
runner. These results do not satisfy the prompt's all-cases correctness gate.

## TS1 comparisons

All results below use `-O3 -DNDEBUG`, double precision, and the same TS1 input
and 10-thread setting on the provided Apple arm64 machine. For the 10,000-cell
seed 2026 comparison, C used **21,388.226 ms** and MICM used **27,943.482 ms**
in one local run.
Both reported 9,000 internal steps, with matching accepted/rejected counts,
function calls, Jacobian updates, decompositions, and solves. All 2,110,000
final cell/species values passed `../reference/compare_final.py`; maximum
relative error was `1.77e-13` (BCARY, cell 6848). The task's hidden seed and
five-run median were not measured here, so this timing is an indicative local
comparison, not the ranked result.

The fixed 100-cell, 100-step, 30-second seed 1 input also passes for one and
four threads. A 1,000-cell seed 2026 run passed all 211,000 values with maximum
relative error `7.74e-14`. Earlier exploratory 100-cell seed 1 comparison had
maximum relative error `9.7e-13`.

## Build verification

Clean `make test-cpu` builds completed with Apple clang 21 and GCC 15.2. Both
compiler runs reported 55 PASS, 613 FAIL, and 12 N/A, and passed the serial and
four-thread TS1 checks. The target exits nonzero while failed cases remain. The
solver_configuration tutorial's DA trajectories still differ from MICM in
one printed species value at 1400 s, although its step statistics match. The
benchmark and library contain no Python linkage. The gate remains incomplete
because many MICM cases have not been ported.
