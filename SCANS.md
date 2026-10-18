# Research log

The experiments compare complete TS1 integration, with unchanged reactions and
species. All accuracies come from the supplied `score.py`; all reported speed
measurements come from C/C++ executables. Python orchestrates runs, checks
independent prototypes and analyzes results. It is never timed as an integrator.
The research is a finite comparison of implemented choices, not a claim of a
globally optimal chemistry solver.

## Controlled design

Each configuration is retained as JSON in `scans/`. Its CSV rows identify the
method/options, scenario, cell/thread counts, accuracy and worst species,
positivity/finite gate, work counters and elapsed integration milliseconds.
Later rows also contain compiler settings, executable hashes and output hashes.
Failures remain in the data. A failed external MICM run can exit before printing
its timing/counters; those fields are left unavailable rather than inventing an
integration time from subprocess wall time.

Within a structural comparison, change one choice from a stated base:
formula family/order, controller, Jacobian age, factor-reuse window,
initial-step reuse, intermediate clipping, extrapolation order selection,
error norm, absolute floor or group size. Tolerance curves then compare work
at measured accuracy. Equal local tolerance does not imply equal global
accuracy, so the family screen alone is not used to declare a speed winner.

| Data | Purpose |
| --- | --- |
| `family64.csv`, `family.json` | Original default plus five adaptive ROS formulas and Euler extrapolation of orders 2–4 at common weights |
| `structure64.csv`, `structure.json` | Fresh versus reused Jacobians, P/PI control, warm starts, stage clipping, step-retention windows |
| `tolerance64.csv`, `tolerance.json` | Historical tolerance screen, including SDIRK2/3; provenance caveat below |
| `tolerance_final64.csv` | Stable-build repetition of the entire tolerance screen after implementation verification and cache improvements |
| `adaptive64.csv`, `adaptive_order.json` | Fixed versus automatic extrapolation order, and factor reuse combined with Jacobian reuse |
| `pi64.csv`, `pi_precision.json` | P/PI work–precision follow-up, avoiding a decision based only on equal tolerances |
| `precision1000.csv`, `full_precision.json` | RODAS3/extrapolation work–precision curves on every disclosed cell |
| `layout256.csv`, `layout.json` | Group sizes 4–128, 256 cells and one thread to isolate layout from core occupancy |
| `error64.csv`, `error_budget.json` | Absolute floors, error norms and a second Jacobian-age check near the recommended accuracy |
| `baselines1000.csv`, `baseline_precision1000.csv` | Original C and supplied C++ MICM three-stage baselines, default and explicit tolerances |
| `operating_gcc1000.csv`, `operating_gcc_final1000.csv` | GCC candidate qualification, before/after an invalid-input guard; includes the subsequently discarded rtol=0.006 point |
| `accurate_gcc1000.csv` | Replacement fast three-digit point, rtol=0.004, on all disclosed cells |
| `stress_references.csv`, `stress_tolerance_gcc.csv` | Converged 48-hour references and follow-up to a failed fast candidate |
| `stress_gcc.csv`, `stress_clang.csv` | Cross-compiler stress qualification; GCC's original candidate file retains the rtol=0.006 failure |
| `rate_cache64.csv`, `rate_cache.json` | Constant-rate caching on/off for RODAS3 and extrapolation; identical states and solver counts |
| `blocking1000.csv`, `blocking.json` | External-step blocks 1/4/16/64 and static/dynamic/guided scheduling; identical states/counts |
| `final_gcc1000.csv`, `final_stress_gcc.csv`, `final_stress_clang.csv` | Final-source cross-compiler qualification of all five points |
| `precision_final1000.csv` | Full tolerance curves with final rate cache and scheduling settings |
| `operating_clang1000.csv`, `operating_clang10000.csv` | Final serialized five-repeat performance/accuracy runs; see RESULTS.md |

The initial 64-cell inputs contain the first 64 supplied cells and exactly their
reference rows. The 256-cell scan uses the analogous prefix. These are screening
sets, not substitutes for the full disclosed validation. The final 1,000-cell
runs use the original scenario files and references directly. The 10,000-cell
inputs reproduce the published recipes, with the supplied first 1,000 cells
verified byte-for-byte. Accuracy is scored on that reference prefix; the
finite/nonnegative gate examines all 10,000 cells. See `tools/scaling.py` and
`results/scaling_inputs.json` for the exact scope.

Solver counters describe group operations; the separately named `rate_updates`
counts individual cell/reaction rate calculations. For example, a solve on 32 cells and one
on four cells both increment `solves` once. Thus operation totals across group
sizes are not interpreted as identical amounts of arithmetic. Timing is the
primary layout comparison; fixed-layout work counts explain numerical changes.

## Findings

**Family and order.** RODAS3 (four stages, three forcing evaluations) gives a
useful low-cost region around two to three digits. The retained ROS3 and ROS4
formulas need more work at the tested accuracies. Four-level linearly implicit
Euler extrapolation uses four LU factorizations and ten solves per trial, yet
its larger accepted steps make it competitive at three-plus digits. Its
seven distinct forcing evaluations share the initial value across rows.
The lower extrapolation orders are cheaper per trial but need more trials at
comparable accuracy. These conclusions concern the tested implementations and
problems, not a universal ordering of solver families.

SDIRK2/3 provide a genuinely nonlinear implicit comparison. Simplified Newton
converges and both methods pass their independent prototype/order tests, but
forcing and nonlinear-iteration work made them substantially more expensive in
the screened tolerance range. Their results are retained rather than replacing
Newton by one iteration and continuing to call the method an implicit solve.
No BDF, multirate, Patankar or mixed-precision performance claim is made.

**Restart cost.** Keeping only the next proposed step across external intervals
is a major improvement. In the initial structural screen, restarting order-four
extrapolation at `dt*1e-4` every time increased summed C time from about 0.99 s
to 3.64 s across the six 64-cell cases. It also changed accuracy, which is why
final selection uses work–precision curves. All chemistry-dependent caches
still reset at each external boundary: changing photolysis never reuses a stale
forcing, Jacobian or factor from the preceding interval.

**Jacobians and factors.** Reusing ROS2 Jacobians preserves its formal order and
saves some work, but ROS2 remains less competitive here. Stale higher-order ROW
Jacobians can noticeably worsen the weakest species: in the initial RODAS3
screen, increasing age from one to four reduced the minimum species digits
from 2.013 to 1.714 for almost no runtime benefit. Extrapolation age two saved
only a few percent near the final operating region. The declared points retain
fresh Jacobians for the clearer accuracy/control behavior. A factor can only
be reused when both its Jacobian and diagonal shift are unchanged. Tests with
an age of four and a step-retention window confirmed actual factor reuse, but
did not produce the best work–precision region.

**Control and automatic order.** PI control increases accuracy and work at the
same local tolerance. Its follow-up tolerance curves did not displace the
selected proportional-controller points. Automatic extrapolation order often
stops at the second-order diagonal on smooth intervals. It saved time at loose
weights but lost enough global accuracy that fixed order four was preferable
at the target levels. This is evidence about the implemented early-acceptance
policy, not a rejection of all variable-order extrapolation algorithms.

**Positivity and error weights.** Intermediate ROS clipping did not improve the
selected frontier. Output projection is retained, measured and charged to the
error norm; ACCURACY.md quantifies its conservation effect. The maximum cell
RMS protects difficult cells from dilution while keeping more useful step sizes
than a full component maximum. A uniform `1e-18 mol m-3` absolute floor provides
margin below the official concentration threshold without species selection.
Floor sensitivity is reported in ACCURACY.md. No method is selected solely
because it passes the finite/nonnegative gate.

**Layout.** At fixed algorithm/weights on one thread, group sizes 4, 8, 16, 32,
64 and 128 took approximately 12.39, 6.83, 5.88, 5.75, 5.93 and 6.30 seconds,
respectively, summed over the six 256-cell scenarios. All were valid and
accurate. Groups of 32 balance SIMD work and workspace/cache size in this
measurement. Final multi-core and 10,000-cell runs retain that size. There are
no architecture flags, intrinsics, vector extensions or compiler tuning sweeps.

**Constant-rate caching and scheduling.** With all conditions and non-photolysis
inputs fixed, a generic dependency list avoids recalculating constant rates.
Across the six 64-cell cases this reduced summed time from 457 to 385 ms for
RODAS3 and from 764 to 689 ms for extrapolation. The 1,000-cell scheduling
scan reduced extrapolation time from 4.855 s (one step, static) to 4.067 s
(four steps, dynamic). Larger blocks gave no further consistent benefit.
The final defaults are rate caching enabled, blocks of four and dynamic
scheduling. These are single-run screening comparisons; repeated final timings
are reported separately. Output hashes and every numerical operation counter
match exactly across cache variants and across all block/schedule variants.
The rate-update count changes as expected. All reactions still contribute at
every solver stage, and every external step still changes photolysis.

**Longer runs changed the selection.** An extrapolation point at `rtol=0.006`
passed level three on every disclosed 1,000-cell case but reached only 2.838
digits in the tropical 48-hour case. Reducing it to `rtol=0.004` restored the
criterion on the stress cases and all disclosed cells. The recommended
`rtol=0.002` point passed throughout. The stress cases informed selection and
are not presented as untouched holdouts. Both their successes and the failed
candidate remain reproducible. There is no access to the experimenter's
withheld scenarios.

## Implementation verification and provenance

`make test-tier3` compares eleven variants with an independent dense Python
prototype on an unrelated nonlinear conservative mechanism, including stiff
rates, Jacobian reuse and both controllers. It checks matching operation counts
and roundoff-level solution differences. Ten fixed-order formulas additionally
show their expected convergence orders against an analytic matrix-exponential
solution; seven-cell layouts match scalar runs. Separate C tests check invalid
states/inputs, backend failures, zero duration and stiff decay with explicitly
accounted projection. Scenario tests compare photolysis interpolation and
snapshot indexing with the supplied C++ runner and compare generated versus
generic kernels for every new CLI method. Both compilers retain all 680 original
cases and all nine original TS1 roundoff-gate variants. Transcripts and numeric
comparisons are in `results/validation-clang/` and `results/validation-gcc/`.

Two development effects are disclosed. The first family scan partly overlapped
compilation/tests, so its elapsed times are exploratory. The early tolerance
scan sampled an executable hash only once, then a build replaced the executable
during startup. Its recorded start hash is retained under that explicit name;
the per-run binary hash is left unknown rather than misattributed. The complete
sweep was rerun as `tolerance_final64.csv` after fixing per-run provenance. No
final performance claim uses the contaminated startup timings.

A later cache improvement shares the initial forcing across extrapolation rows
and retries and avoids redundant same-state Jacobian refreshes. It reduces work
without changing the formulas; prototype agreement was rerun afterward. The
historical CSVs retain their original operation counts. Current implementation
claims use the stable-build repetition and final operating-point runs.

Final runs are serialized and summarized by medians of five. They do not overlap
another build or benchmark from this task; general system load is uncontrolled.
Repeated output is scored again unless its SHA-256 exactly matches an already
scored output for the same scenario and level. Such reuse is recorded explicitly.
Plots and summary tables are generated from the retained CSVs by `tools/report.py`.
