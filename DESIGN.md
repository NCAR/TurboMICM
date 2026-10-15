# Design

`include/chem.h` defines the mechanism, reaction-rate and state API.
`include/chem_solver.h` defines sparse factorization, Rosenbrock parameters and
solver extension points. The library owns copies of input species/reaction
metadata. States and solvers borrow the mechanism and own their numerical
storage. Destroy them before destroying the mechanism. A single state/solver
pair is used by one thread at a time; many pairs can share the same mechanism.

## Configuration and kinetics

`tools/configure.py` reads included CAMP records, validates names and finite
numbers, expands branched reactions, converts rate units and generates a C
builder using the public API. Generated C is compiled by the Makefile with the
selected compiler and optimization flags. The runtime needs no parser or Python.
TS1's data is retained in `mechanisms/ts1.json`; `tools/translate_ts1.py` reproduces
it from the supplied read-only input. The three tutorial configurations are
checked against independent C builders at one and seven cells, including ten
integration intervals. Additional parser error tests remain to be implemented.

`src/rates.c` supplies rate defaults, derived parameters and evaluations.
`src/mechanism.c` performs a stable sort by MICM's rate-type order, retaining
reaction order within each type. Parameterized species have no solution-vector
entry; their values multiply rate constants once per reactant occurrence.

Each variable is contiguous across cells: `values[variable * cells + cell]`.
Forcing evaluates each reaction's rate product into a cell workspace, then
applies reactant losses and product yields in source order. OpenMP SIMD acts on
contiguous cell loops. Dense Jacobian evaluation is available as a reference
path. Solver construction prepares a sparse derivative plan ordered by
independent variable and reaction, avoiding name lookups during integration.

State copies own independent buffers and borrow their mechanism. Copy assignment
allocates before replacing the destination; moves transfer the owning pointer.
`chem_state_set_concentration` and `chem_state_set_custom` require exactly one
value per cell. Their `chem_state_fill_*` counterparts broadcast explicitly.
Bulk named setters use `ChemNamedValues`; the custom-matrix setter transposes
checked row-major input into the library layout. Named concentration accessors
return contiguous columns for mutable or const access.

## Matrix utilities

`include/chem_matrix.h` exposes owning dense matrices and borrowed array/column
views. Group zero selects row-major storage; a positive group size stores each
column contiguously within the group and pads the final group. Real and native
integer storage share shape and ownership operations. Scalar fill/min/max touch
allocated padding, while element arithmetic visits logical entries.

Reusable array functions capture a callback and caller-owned context, with a
schema derived from prototypes. Invocation checks types, kinds, columns,
mutability and agreement of row counts before calling user code. Row counts may
change between calls. Views support copies, fills, row callbacks, and reductions
with an explicit strict flag controlling padding. Const views reject writes.
These utilities provide the C equivalent of MICM's matrix-policy operations.

## Linear algebra and integration

`src/sparse.c` adds symbolic fill once and stores a plan for Mozart elimination.
It implements in-place Mozart and Doolittle decomposition and forward/backward
substitution, retaining the original arithmetic order. `chem_markowitz` selects
the minimum degree-product pivot and accounts for fill; ties use original index
order. `chem_mechanism_reorder` applies this permutation before state creation.
`chem_lu_factor_copy` preserves the matrix and puts factors in a separate
caller-owned buffer with the same combined layout. `chem_lu_lower` and
`chem_lu_upper` expose the triangular factors, including the implicit unit lower
diagonal. Separate-factor and in-place paths share the arithmetic kernels.

`src/rosenbrock_parameters.c` contains the five original coefficient sets.
`src/rosenbrock.c` implements the common stage loop, weighted RMS error and
step-size controller. Initial step size, rejection behavior, tolerances and
statistics follow MICM's ODE path. Backward Euler in `src/backward_euler.c` shares the sparse linear algebra and
implements the original Newton convergence criterion and step reductions.
Built-in linear and equilibrium constraints are installed with
`chem_mechanism_set_constraints` before creating states or solvers. Definitions
are copied and species names resolved to indices. Algebraic rows replace kinetic
rows in forcing and Jacobians, have zero mass-matrix diagonal, and use their
step change in the error estimate. Equilibrium constants are updated alongside
rates from the original Van't Hoff formula.

`chem_initialize_constraints` performs the original weighted Newton projection
and damped line search before integration. It modifies only algebraic variables
and restores the full input state on failure. Projection uses the same sparse
factor/solve callbacks as integration. `chem_add_mass_diagonal` and the finite-
difference utilities in `chem_jacobian.h` expose independently testable numerical
operations. External and aerosol models remain to be implemented; this is not
yet the complete CPU port.

To add a Rosenbrock method, add a method enum and a coefficient-set branch in
`chem_rosenbrock_parameters`, then test the coefficients and stiff integration
behavior. The stage implementation reads coefficients from the passed structure
and supports up to six stages without method-specific branches.

To replace step-size control, supply a `ChemStepController` through
`chem_solver_set_controller`. It returns the multiplicative step adjustment;
acceptance/rejection and bounds remain in the solver. The default implements
MICM's safety factor, local order and minimum/maximum adjustment factors.

To replace linear algebra, provide factor and solve functions and an explicit
context through `chem_solver_set_linear_solver`. They receive the symbolic
pattern, cell count and matrix storage. The default factor is Mozart. A custom
adapter can select Doolittle or another experimental implementation without
changing the stage formulas.

To add a rate type, add its parameter member and tag in `chem.h`, defaults,
preparation and evaluation in `rates.c`, a configuration mapping, and independent
reference tests. A callback rate takes conditions and a context pointer; the
caller manages that context's lifetime.

## Performance parameters

| Parameter | Default | Control |
| --- | --- | --- |
| Cell group size | 32 | `CHEM_GROUP_SIZE` environment variable; positive integer |
| Threads | OpenMP maximum, normally all cores | benchmark `--threads N` or OpenMP environment |
| Group scheduling | OpenMP static | fixed policy; each group has its own adaptive state |
| SIMD width | compiler-selected | portable `omp simd`, no target-specific width |
| LU algorithm | Mozart | linear-solver callback; Doolittle is also provided |

A final short group contains only its real cells. Error normalization includes
every real cell and excludes no species. Each group performs the same warm-up
and timed steps. Statistics sum over groups, so group size changes their totals.
The benchmark uses the same generic mechanism and solver functions as the tests;
there is no hand-written TS1 numerical kernel. Generated C currently supplies
mechanism data, not specialized arithmetic.
