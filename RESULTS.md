# Verification in progress

These are intermediate results, not completion of the requested correctness gate.
The full 680-case suite is not yet ported.

With Clang and `OPTFLAGS=-O3 -DNDEBUG`:

- 32 rate-constant and lambda-override cases passed their ported assertions.
- The C TS1 run on the supplied 100-cell seed-1 input, 100 timed 30-second steps
  after one warm-up, passed all 21,100 values against the supplied reference.
  The largest relative difference was 9.7e-13; the gate is rtol=1e-10, atol=1e-30.
- The same check passed with one thread/groups of 32, four threads/groups of 7,
  and the generic table-driven chemistry evaluator.
- One preliminary generated-kernel run with one thread/groups of 32 took
  405.04 ms. This is a development observation, not a controlled benchmark or a
  ranking measurement. Timed statistics were 3,600 internal/accepted steps,
  zero rejected, 7,200 function calls, 3,600 Jacobians/decompositions, 10,800 solves.

The supplied C++ reference executable was built using `make micm-ref` with the
same optimization flags. Fresh comparative timing, the 10,000-cell run, complete
case coverage, and the portability audit remain outstanding.

With GCC 15 and the same flags, all 32 currently mapped cases passed. TS1 also
passed the supplied 100-cell gate with four threads and groups of seven; the
largest relative difference was 9.69e-13.

The subsequent constraint checkpoint passes 86 mapped cases with both compilers:
32 rate/lambda cases and 54 linear, equilibrium, set, and initialization cases.
These include rollback on failure, damped Newton cold starts, row scaling,
multi-cell indexing, analytic residuals/Jacobians, and finite-difference checks.
The TS1 benchmark has not yet been rerun after that checkpoint.

The matrix checkpoint passes all 420 implemented cases with Clang and GCC 15
under `-O3 -DNDEBUG`: 117 dense/vector, 190 sparse, 27 grouped-view/category,
and the previous 86 chemistry/constraint cases. These cover the original lane
widths, CSR/CSC indices, exact printed matrix output, native integer arrays,
mutable/const callback views, dimension errors, temporary storage, and reductions.
There are 260 cases still unported, and no N/A declarations.

The linear solver checkpoint passes 487 implemented cases with Clang and GCC 15.
This adds all 52 LU/linear-solver cases, ten solver-core cases and five version
cases. Tests preserve the reference random inputs and tolerances, including
1e-12 factor reconstruction checks. A separate L/U backend produces the same
accepted/rejected step counts and concentrations as the combined backend.
The post-change Clang TS1 run (four threads, groups of seven) passes all 21,100
reference values with largest relative difference 9.7e-13. There are 193 cases
still unported; full-suite and performance validation remain outstanding.

The species/system checkpoint passes all 501 mapped cases with Clang and GCC 15.
Fourteen added cases cover the original species, phase-species, phase, and system
assertions. There are 179 cases still unported and no N/A declarations.

The state/process checkpoint passes 550 cases with both compilers. The 49 new
cases comprise state (22), process (4), reaction rate store (11), process set (5),
and finite-difference verification (7). Comparisons use the original tolerances,
including four-ULP matrix entries and the 1e-7 absolute/relative Jacobian gate.
All three random process-set sizes execute and additionally match a scalar
cell-major forcing evaluator. The full suite still has 130 unported cases.

The integration checkpoint passes all 606 mapped cases on Clang and GCC 15 with
`-O3 -DNDEBUG`. All 32 analytical cases execute their full original solver and
matrix-layout combinations, including the super-stiff chains, Robertson, E5,
HIRES and Oregonator reference trajectories. Twelve additional DAE integration
cases check equilibrium residuals, conservation and overshoot bounds. The eight
builder cases, eight-thread shared-solver test, integrated reaction rates, and
both Terminator variants also pass. There are 74 unported cases and no N/A cases.

The Chapman/regression checkpoint passes all 634 implemented cases on both
compilers. Generated and hand-built Chapman mechanisms both pass every regression
against the independent legacy scalar oracle, with original inputs, layouts and
tolerances. Both flow-tube solvers match all 3,601 supplied rows at atol=1e-10;
the original test checks only the first ten. Forty-six cases remain unported.

The external-model checkpoint passes all 672 implemented cases with both compilers.
Eighteen aerosol cases run every original solver/layout variant; twenty external
constraint cases check shared algebraic rows, finite differences, temperature
updates, conservation and independent ownership when rebuilding a description.
Only the eight tutorial programs remain unmapped; there are no N/A declarations.
