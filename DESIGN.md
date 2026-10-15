# Architecture

The public API is `include/chem.h`. A mechanism contains immutable species and
reaction descriptors plus a symbolic sparse factorization plan. Each state owns
its conditions, concentrations, rate constants, tolerances, work arrays, and
linear-solver selection. States sharing one mechanism have no shared integration
workspace. Borrowed descriptor arrays and callback contexts have caller-managed
lifetimes. Named callback rate overrides are applied before parallel state updates.

`tools/configure.py` reads YAML or JSON, follows includes with cycle detection,
normalizes reaction lists and stoichiometry, applies unit conversions, and
optionally reorders species using stable diagonal Markowitz ordering. It emits
C descriptors compiled by the Makefile with CC and OPTFLAGS. The numerical path
is the same for generated descriptors and descriptors assembled in C.

`tools/translate_ts1.py` extracts species, reactions, rate constants, and initial
values from the supplied TS1 declarations. It does not compile or execute that
source. The resulting JSON explicitly selects molar rate units. Counts are
validated during translation; the tests verify that retranslation matches the
committed configuration.

## Numerical components

- `rates.c`: the reference rate formulas and default parameter values.
- `mechanism.c`: descriptor validation, sparsity and fill construction, and
  elimination/derivative operation plans.
- `process.c`: mass-action forcing and analytic negative Jacobian. Repeated
  reactants retain their multiplicity and multiplication order. Algebraic rows
  omit built-in chemical forcing; external model Jacobian contributions remain.
- `linear.c`: in-place Mozart elimination and Doolittle factorization with the
  reference operation ordering, and forward/back substitution without pivoting.
- `parameters.c`: all five Rosenbrock parameter sets.
- `rosenbrock.c`: stage evaluations, the common weighted error norm, acceptance,
  rejection, and step-size policy. Public solves clamp differential species
  afterward, including after partial integration, matching the reference wrapper.
- `backward_euler.c`: Newton iteration, positivity, convergence test, and retry
  reductions from the reference.
- `constraints.c`: linear, equilibrium, and callback residuals/Jacobians, and
  transactional damped Newton initialization. Algebraic errors use the step
  change; stage coupling and alpha shifts include only differential rows.

Storage is species-major inside each cell group: every operation on a species
visits contiguous cells. Sparse entries use the same layout. Elementwise loops
use `restrict` and OpenMP SIMD; error reduction retains a strict order. No
intrinsics, compiler vector types, target tuning, or relaxed floating point is
used.

## Array and sparse utilities

`chem_array.h` and `chem_sparse.h` provide the C counterparts of the reference's
CPU matrix policies. Double, integer, and boolean elements are selected at run
time. Group length zero uses row-major storage; positive lengths interleave rows
within each group. The final group is padded. Sparse matrices store only unique
structural coordinates in CSR or CSC order, with a shared pattern across blocks.
Their values are an array of blocks by nonzeros and use the same grouped layout.

A `ChemView` records its base pointer, logical length, padded capacity, and
strides within and between blocks. Views borrow storage. Group-scoped views
precompute their base and provide contiguous access for interleaved layouts.
The common view operations implement fill, overlapping-safe copy, callback
iteration, and reduction. Strict reduction excludes padding; ordinary reduction
includes it. Read-only flags prevent declaring a const view writable.

`ChemArrayFunction` captures column counts, element types, group length, and
argument mutability. Creation permits different row counts; invocation requires
equal logical row counts among all arguments and compatible captured shapes.
Rows can change between invocations. Callback element pointers preserve aliasing.
Callbacks own their context and must obey the declared mutability mask. The
function metadata is immutable and each invocation owns its argument workspace.

Dense copies and AXPY can operate across layouts directly. Callback functions
require a common effective group length (row-major and group length one are
compatible), preserving the reference's ordering validation. Sparse descriptors
borrow no coordinate arrays after construction: the library owns their maps,
compressed indices, and numeric storage. Clones own independent allocations.

## Tuning parameters

| Parameter | Default | Selection |
|---|---|---|
| Cells per state/group | 32 | benchmark `--group-size N`, or `CHEM_GROUP_SIZE` at build |
| Threads | OpenMP maximum, capped by available groups | `--threads N`, OpenMP environment |
| SIMD length | compiler choice | portable `omp simd`, no prescribed hardware width |
| Group scheduling | static contiguous group assignment | OpenMP schedule in `bench.c` |
| Species reordering | enabled for generated data | generator `--no-reorder` |

Every positive group size and thread count is intended to preserve correctness;
partial groups allocate only their actual cell count. For small inputs the
default group size shrinks enough to keep all requested threads occupied; an
explicit group size is honored as given. Statistics count solver
operations on a whole group, so changing the grouping changes their sums. Solver
tolerances and numerical method constants are accuracy controls, not performance
tuning knobs, and retain reference defaults.

## Extension points

To add a Rosenbrock method, return a `ChemRosenbrock` parameter structure with
its stage tables, function-reuse flags, and embedded estimator order. The stage
loop supports up to six stages. Existing sets are examples in `parameters.c`.

To replace the step-size controller, set `controller` and `controller_context`
on the parameter structure. The callback returns the next-step size factor;
acceptance and rejection history remain in the integrator. A null callback uses
the reference safety factor and bounded power law.

To swap the linear solver, pass a `ChemLinearSolver` to
`chem_set_linear_solver`. The factor callback receives the sparse matrix storage
and cell count; solve overwrites its right-hand side. The public index API maps
rows/columns into that storage. Mozart is the default; Doolittle exercises the
same interface.

To add a rate type, extend the enum and default factory in `rates.c`, add its
calculation branch, and add its named fields and unit conversion in
`tools/configure.py`. Add direct numerical assertions in `tests/rate_cases.c`
and a generated-versus-code mechanism check. A C function/context rate supports
custom calculations without adding an enum value.

External models provide declared Jacobian entries and callbacks for parameter
updates, forcing, and negative Jacobian accumulation. Callback constraints
provide a residual and optional dense gradient with declared dependencies.
