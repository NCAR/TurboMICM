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

## Tier 3 checkpoint 3

Completed full 1,000-cell work–precision screening, grouping/error-weight scans,
MICM C++ baselines, and converged references for three new 48-hour scenarios.
The reference pairs agree to 10.23–10.49 digits overall. GCC passes all 680
original cases, all nine original TS1 roundoff checks, all prototype/order/
projection/stiffness/I/O tests, and the positivity/finite gate for five candidate operating points on every
disclosed and stress scenario. The rtol=0.006 extrapolation candidate missed
level three on the tropical 48-hour case (2.838 digits); the recommendation
passed. This failed accuracy result is retained and motivates a tighter
fast three-digit candidate. Detailed CSVs and sanitized validation transcripts
are retained. The candidate recommendation is extrapolation order four,
rtol=0.002, atol=1e-18, fresh Jacobians, groups of 32.

One metadata defect was found in an early tolerance scan: its executable hash
was sampled before a concurrent build finished. The historical start hash is
retained with an explicit warning, not asserted as each row's actual binary.
The scan harness now hashes the executable per run; a stable-build repetition
of that sweep remains to be performed. Final Clang gates, repeated isolated
timings, 10,000-cell scaling, plots/reports, and completion audit remain.

## Checkpoint 4: final driver and numerical selection

The final five numerical configurations are frozen: RODAS3 at rtol 0.01/0.003
and order-four extrapolation at 0.004/0.002/0.001, all atol 1e-18. Generic constant
rate caching and four-step dynamic scheduling preserve bit-identical states and
solver counts in controlled scans. Full compiler validation, serialized five-run
1,000/10,000-cell timings, plots and final reports follow this checkpoint.

## Checkpoint 5: implementation frozen and final gates passed

Integration sources are frozen at commit 07135fb. Both final compiler builds
passed all 680 original cases, all nine TS1 roundoff checks, 44 prototype
comparisons, convergence/stiffness tests, rate dependency tests and scenario
I/O/cache/blocking checks. Final-source GCC qualification passed every declared
point on all six disclosed and three stress scenarios; Clang passed all stress
scenarios. The full final-driver tolerance sweep and serialized five-repeat
operating measurements are being retained under their final scan names.
Reports explicitly distinguish full 1,000-cell accuracy from reference-prefix
accuracy in the 10,000-cell scaling workload.

The 84-run final-driver tolerance sweep and all 210 repeated 1,000-cell runs
are complete. Every declared point meets its level in every repetition. The
recommended point's six median times sum to 4.249 s, versus 10.101 s for matched
C MICM3 and 31.916 s for matched C++ MICM3: 2.38x and 7.51x, respectively. Its
minimum disclosed accuracy remains 3.448 overall / 2.550 worst-species digits.
The 150-run 10,000-cell measurement series is in progress with unchanged binaries.

## Completion

All 210 final 1,000-cell and 150 final 10,000-cell repeated measurements are
complete. Every declared point passes its intended level on the reference
domain and its finite/nonnegative gate on every output cell. Repetitions have
identical output hashes. Ten-thousand-cell accuracy remains explicitly limited
to the independently referenced first 1,000 cells. Generated forcing tables and
active scenario settings match the supplied cases, and input/reference hashes
are retained.

The final reports include all five concrete commands, full per-scenario
accuracy and timing tables, work counters, stress results, per-species errors,
and three visually checked PNG/SVG figures. Recommended extrapolation reaches
3.448 overall / 2.550 worst-species digits in the disclosed set; its largest
individual scored error is 0.998%. Its six 10,000-cell median times sum to
40.364 seconds. The numerical source has not changed since the complete
Clang/GCC validation. The deliverable/provenance/portability audit is retained
in results/completion_audit.txt. No required work remains.
