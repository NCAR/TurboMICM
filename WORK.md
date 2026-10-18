# Completion record

The implementation covers all 680 supplied cases: 672 unit/integration cases
and eight tutorial programs, with no N/A declarations. Both Clang and GCC 15
pass the full suite under `-O3 -DNDEBUG`.

Implemented and checked:

- Generic species, phases, processes, states, owned solvers, and all rate types.
- Dense/vector/sparse layouts, native integer storage, views and callbacks.
- Both LU algorithms, separate and combined factors, and replaceable backends.
- All five Rosenbrock parameter sets, backward Euler and controller replacement.
- Built-in/external algebraic constraints, Newton initialization and external
  model composition, parameter updates and independent context ownership.
- Full analytical/regression/flow-tube cases and thread-shared immutable solvers.
- JSON/YAML reading, all three tutorial configuration equivalence checks,
  automatic kernel generation and reproducible TS1 translation.
- Standalone tutorial execution and comparison of all deterministic output.
- C TS1 CLI, input defaults and overrides, warm-up, timing, statistics and CSV.

Final TS1 checks pass on both compilers for seeds 1, 314159 and 2026, including
the full 10,000-cell gate. Clang median time is 7,306.806 ms over five runs; the
fresh MICM measurements are 143,929.029 ms (one thread) and 28,084.099 ms (ten).
Group-size and scheduling measurements support retaining the default of 32 cells
with static scheduling. All requested deliverables are present.
Earlier incremental checkpoints are retained in git history. RESULTS.md is the
source for final measurements and numerical gate evidence.

## Tier 3 checkpoint 1

Implemented scenario I/O/timing, separate adaptive Rosenbrock and extrapolation
paths, and reproducible logged structure scans. `scans/family64.csv` contains all
54 initial configurations/scenarios, including baseline failures. Screening
costs are exploratory: the early rows overlapped compilation/tests. All 680
Clang Tier 1 tests pass, as does the default 100-cell TS1 roundoff gate (maximum
relative difference 4.25e-14). Prototype equivalence, broader scans, final
operating points, both-compiler final validation, and updated reports remain.

## Tier 3 checkpoint 2

Added SDIRK2/3 with simplified Newton and independent Python prototypes.
Forty-four adaptive comparisons now agree within 1.14e-14 absolute with identical
function/Jacobian/LU/solve/reuse/nonlinear/order counts; analytical convergence
orders and conservation pass for ten fixed-order formulas, including grouped
seven-cell layouts. All tests use a small mechanism unrelated to TS1.

Completed family, tolerance, and structural scans. RODAS3 favors lower accuracy;
four-level Euler extrapolation favors three-plus digits. SDIRK was substantially
more expensive in the screened range. Added automatic extrapolation order as a
further structural comparison; its scan is in progress. Eliminated duplicate
initial forcing in extrapolation and redundant refreshes after same-state
rejections. The earlier CSVs remain historical evidence with their original
counts. Full disclosed-cell scans, final gates, stress cases, operating points,
and research reports remain.
