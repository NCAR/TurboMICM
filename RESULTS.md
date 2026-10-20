# Results

## MICM case list

`make test-cpu` prints one line for each of the 680 IDs in
`../reference/micm_cases.txt`, followed by aggregate and TS1 reference-check
lines. Results for the current source are:

| Compiler | PASS | FAIL | N/A | TS1 serial / parallel |
| --- | ---: | ---: | ---: | --- |
| Apple clang 21, macOS arm64 | 612 | 54 | 14 | PASS / PASS |
| GCC 15.2, macOS arm64 | 612 | 54 | 14 | PASS / PASS |

The 612 passing case lines cover implemented rate-constant formulas, rate
parameter descriptor references, the two configuration-versus-hand-built
examples, the README chemistry example, multiple grid cells, the OpenMP tutorial,
the vectorized storage example, the solve_results and solver_configuration
tutorials, all 22 state cases, the Chapman descriptor and forcing checks plus
four generic solver integrations, six Chapman rate/forcing/Jacobian regression
cases, three packed Chapman sparse linear-algebra checks, 117 dense and
vector-matrix operations, 190 standard and vector-ordered sparse matrix cases, 21 grouped dense and sparse view operations, the analytic terminator checks, Chapman
oxygen conservation, integrated intermediate monotonicity, condition-dependent
callback rates and override, seven Jacobian verification cases, C array
iteration, 14 owning species/phase/system cases, nine reaction/process-set
construction, forcing, algebraic-row masking, and Jacobian cases, 42 linear and
equilibrium constraint construction/residual/Jacobian cases, six equilibrium
DAE integration cases, the linear-constraint Terminator and Robertson DAE case,
five DAE overshoot and error-control cases, and the OpenMP TS1 path. Twelve algebraic
initialization cases cover consistent and
inconsistent states, multi-cell projection, rollback, row scaling, weighted
tolerances, and damped Newton updates. Thirty-two LU decomposition cases verify
Mozart and Doolittle reconstruction for dense, diagonal, and random matrices;
20 linear solver cases cover separate and packed LU solves. The C API also
includes a fully implicit backward Euler integrator with Newton updates,
adaptive step reduction, and the MICM residual convergence rule; its solver and
convergence cases pass. Twelve ordinary and twelve super-stiff analytical
A-to-B-to-C cases now pass across backward Euler and Rosenbrock, including Troe,
photolysis, ternary activation, tunneling, Arrhenius, and branched rates. Eight
Rosenbrock utility cases cover diagonal alpha updates, per-species tolerances,
normalized error over all variables, and agreement after rejected steps across
both LU implementations. The DA tutorial output table was corrected against the
copied MICM tutorial executable: both DA methods print `0.703` at 1400 s. The 14 N/A
cases cover a rate-store callback entry, four compile-time view-category traits,
the grouped-view compile-time category check, and release metadata with no
meaningful equivalent in the current C API, including three builder-specific
checks. The 54 remaining FAIL cases cover eight analytical integration cases,
18 aerosol integrations, 20 external-model constraint cases, two regression
cases, and six regression solve cases. Per-case reasons are printed by the
runner. These results do not satisfy the prompt's
all-cases correctness gate.

## TS1 comparisons

All results below use `-O3 -DNDEBUG`, double precision, and the same TS1 input
and 10-thread setting on the provided Apple arm64 machine. For the 10,000-cell
seed 2026 comparison, C used **21,199.606 ms** and MICM used **28,453.001 ms**
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
compiler runs reported 612 PASS, 54 FAIL, and 14 N/A, and passed the serial and
four-thread TS1 checks. The target exits nonzero while failed cases remain. The
benchmark and library contain no Python linkage. The gate remains incomplete
because many MICM cases have not been ported.
