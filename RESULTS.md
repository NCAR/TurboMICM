# Results

## MICM case list

`make test-cpu` prints one line for each of the 680 IDs in
`../reference/micm_cases.txt`, followed by aggregate and TS1 reference-check
lines. Results for the current source are:

| Compiler | PASS | FAIL | N/A | TS1 serial / parallel |
| --- | ---: | ---: | ---: | --- |
| Apple clang 21, macOS arm64 | 47 | 621 | 12 | PASS / PASS |
| GCC 15.2, macOS arm64 | 47 | 621 | 12 | PASS / PASS |

The 47 passing case lines cover implemented rate-constant formulas, rate
parameter descriptor references, the two configuration-versus-hand-built
examples, the README chemistry example, multiple grid cells, the OpenMP tutorial,
the vectorized storage example, and the OpenMP TS1 path. The 12 N/A cases cover
C++ callable storage, template type traits, and release metadata with no
meaningful equivalent in the current C API. The remaining 621 cases are
explicitly reported as FAIL; they cover unimplemented MICM builder/state APIs,
generic matrix and view utilities, standalone linear algebra cases,
constraint/DAE and alternate integration paths, external/aerosol model APIs,
and one tutorial statistics mismatch. The per-case reasons are
printed by the runner. These results do not satisfy the prompt's all-cases
correctness gate.

## TS1 comparisons

All results below use `-O3 -DNDEBUG`, double precision, and the same TS1 input
and 10-thread setting on the provided Apple arm64 machine. For the 10,000-cell
seed 2026 comparison, C used **20,465.211 ms** and MICM used **27,195.973 ms**.
Both reported 9,000 internal steps, with matching accepted/rejected counts,
function calls, Jacobian updates, decompositions, and solves. All 2,110,000
final cell/species values passed `../reference/compare_final.py`; maximum
relative error was `2.27e-13` (BCARY, cell 4179). The task's hidden seed and
five-run median were not measured here, so this timing is an indicative local
comparison, not the ranked result.

The fixed 100-cell, 100-step, 30-second seed 1 input also passes for one and
four threads. A 1,000-cell seed 2026 run passed all 211,000 values with maximum
relative error `7.74e-14`. Earlier exploratory 100-cell seed 1 comparison had
maximum relative error `9.7e-13`.

## Build verification

Clean `make test-cpu` builds completed with Apple clang 21 and GCC 15.2. Both
compiler runs reported the same 47/621/12 case split and passed the serial and
four-thread TS1 checks. The target exits with status 1 because the 621
unported cases are reported as FAIL. The `solve_results` tutorial also remains
FAIL: its C run accepts 17 steps for the 200-second interval while MICM accepts
18, although the TS1 benchmark statistics match. Neither build emitted compiler
warnings. The benchmark and library contain no Python linkage. The gate remains
incomplete because the majority of MICM test cases have not been ported.
