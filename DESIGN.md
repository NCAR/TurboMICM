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
integration intervals. Forty-nine parser validation checks cover malformed
records, names, units, stoichiometry, include cycles and finite numbers. Three
additional generated rate types are checked against independent formulas, and
all 379 generated TS1 initial fields match the reference exactly by name.

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

`include/chem_system.h` adds owning species properties, affine parameterizations,
phase entries with optional diffusion/density, systems, and named reaction
metadata. Constructors and assignment helpers make deep copies. Name queries
exclude parameterized species and optionally apply caller-defined ordering.
`chem_system_mechanism` resolves names into the same indexed reactions used by
the low-level API; the mechanism owns a metadata snapshot that keeps parameter
contexts alive. Callers can destroy all source metadata after construction.

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

`include/chem_sparse_matrix.h` adds declared-structure CSR/CSC matrices with a
builder, deduplicated entries, standard or grouped block layouts, checked lookup,
diagonal operations and printing. This structure does not add the diagonal or LU
fill. Its values reuse `ChemMatrix` storage and its block views reuse `ChemView`.
Sparse functions allow different matrix structures, enforce matching group
widths, and check equal block/row/vector lengths on each invocation.

`chem_view_group` narrows a raw view to one group with a precomputed base pointer
and distinct grouped category. Grouped fill/copy and default grouped row loops
include the allocated tail padding; `chem_views_for_each_strict` excludes it.
The logical length remains available for bounds-checked element access. Raw
views and caller-owned vectors retain logical-only iteration.

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

`chem_sparse_multiply` evaluates an unfactored sparse matrix against column-major
cell arrays, using the same sparse index plan. Its output must not overlap either
input. This also lets the old Chapman matrix/vector fixture exercise the generic
library without introducing a mechanism-specific kernel.

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
operations.

`chem_external.h` defines external process and constraint callbacks. The model
constructor copies added species, resolves named variables, parameters and
Jacobian entries, and owns one context per model. Plain contexts are copied by
size; contexts with pointers supply clone/destructor hooks. Bindings are passed
into callbacks and remapped by mechanism reordering. Callbacks use column-major
state arrays and `chem_model_jacobian_add`, which accepts positive derivatives
and applies the dense or sparse destination's sign convention.

The composed forcing/Jacobian paths first evaluate kinetic and external
processes, clear algebraic rows, then apply built-in and external constraints.
External process structural entries remain present even on algebraic rows, so
model derivative callbacks can always address their declared entries. Constraint
projection and the mass matrix use the combined algebraic mask. Separate process
and constraint APIs support finite-difference verification. Parameter updates
run for every model, including a constraint model preceded by a process-only
model. Independent states can share the immutable mechanism and model contexts;
callbacks must keep mutable per-cell values in the state.

The integration tests implement two aerosol models through this public callback
interface, with added phase species, number concentrations, temperature-updated
parameters and user-set parameters. Both solvers use the same mechanism and
sparse kernels. All 680 reference cases are exercised by the test runner.

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

## Solver plans

`chem_builder.h` offers an owning high-level construction path over the same
mechanism and numerical solvers. A builder borrows input metadata during build;
its resulting plan owns copied system, reactions, constraints and model contexts.
It resolves names, optionally rejects unused variables, reorders the mechanism
and finalizes model indices against the sparse pattern. Dense/sparse storage
descriptors expose the canonical cell layout and the configured cell limit.
A full model finalizer receives both descriptors; a sparse-only callback remains
available for adapters that do not use dense storage information.

Plans create independent states and workspaces, enforce an optional cell limit
and retain Rosenbrock or backward Euler options. Passing an override to
`chem_solver_plan_solve` updates the retained options; passing NULL reuses them.
Free states/workspaces before their plan. Shared plans may be read concurrently,
but replacing retained options requires serialization. The lower-level APIs
remain available when the caller manages mechanism ownership and parameters.

## Performance parameters

| Parameter | Default | Control |
| --- | --- | --- |
| Cell group size | 32 | `CHEM_GROUP_SIZE` environment variable; positive integer |
| Threads | OpenMP maximum, normally all cores | benchmark `--threads N` or OpenMP environment |
| Group scheduling | OpenMP static | `OMP_SCHEDULE`; runtime static, dynamic or guided policy |
| SIMD width | compiler-selected | portable `omp simd`, no target-specific width |
| LU algorithm | Mozart | linear-solver callback; Doolittle is also provided |

A final short group contains only its real cells. Error normalization includes
every real cell and excludes no species. Each group performs the same warm-up
and timed steps. Statistics sum over groups, so group size changes their totals.
The benchmark uses the same generic mechanism and solver functions as the tests;
there is no hand-written TS1 numerical kernel. Generated C currently supplies
mechanism data, not specialized arithmetic.

## Build-time mechanism specialization

`tools/specialize.py` shares `tools/configure.py`'s validation and normalized
mechanism representation. Its passes are:

1. Stable rate-type sorting and custom-parameter indexing, excluding
   parameterized species from the solution vector.
2. Markowitz ordering with the same ties and unsigned-cost ordering as the
   generic implementation, followed by complete symbolic Mozart fill.
3. Ordered forcing and derivative operations, with literal matrix/species
   indices and yields. Repeated reactants remain repeated operations.
4. Emission of fixed-parameter rate formulas, pivot operations, forward/backward
   substitution and mass-diagonal updates as C11 cell loops.
5. Instantiation of the common Rosenbrock template with direct chemistry calls,
   a known three-stage loop and a fused final-state/error-estimate loop.

`Plan` owns the symbolic representation; `Emitter` handles C loops and function
boundaries. Separate rate, chemistry, factor, solve and integrator emitters keep
numerical expressions separate from formatting. Adding an existing rate-type
reaction requires only a configuration edit. Adding a new rate type requires
reader, generic evaluator and specialized rate-expression support, plus tests.
Generated source dependencies include every transitive CAMP input.

Species ordering and sparse offsets are known during generation. The runtime
builder still creates ordinary generic-library metadata and reproduces the
ordering once, allowing shared allocation, named setters and CSV I/O. This setup
is outside benchmark timing. Generated kernels never traverse reaction objects,
derivative plans or sparse index maps in their full-width arithmetic. The two
prepared branched constants and surface molecular-speed constant are read from
immutable metadata; C computes them once so Python/libm differences cannot enter
these constants.

Temperature, pressure, density, all custom parameters, concentrations, interval
lengths and tolerances remain runtime inputs. A per-state exact-input cache holds
copies of the three condition fields and all custom parameters. A cache hit
requires byte-identical inputs, including signed zero. A miss evaluates rates
normally and replaces the snapshot. Allocation failure simply disables caching.
The generic updater invalidates the cache; state copies start without a cache,
and destruction frees it. No Jacobian or factorization is cached. Concentrations
do not belong in this cache because CAMP rate constants depend only on conditions
and custom parameters; concentration products are recomputed at every stage.

The generated mechanism must remain immutable and the solver must be attached to
that same mechanism and ordering. Generic external processes, constraints and
runtime mechanism edits remain available through the generic API. The supported
configuration reader does not construct these callback-based features. The
specialized three-stage path removes their empty checks. Other Rosenbrock stage
counts delegate to the generic integrator, retaining the attached exact algebra.

The common template preserves MICM's stage equations, ordered RMS sum, acceptance,
rejection, step-size bounds, statistics and final concentration clamp. Fusing the
three final-state additions and the three error additions changes traversal, not
the arithmetic order for an element. SIMD applies only across independent cells;
there is no floating-point SIMD reduction in the adaptive error norm. The LU
singularity check reduces an integer flag. No precision or compiler math mode is
changed.

### Parameters

| Setting | Default | Meaning and valid alternatives |
| --- | --- | --- |
| `GROUP_SIZE` | 32 | Positive default group size for specialized drivers; also default generated stride |
| `GEN_WIDTH` / `--width` | `GROUP_SIZE` / 32 | Positive fixed stride, or 0 for variable-stride generated kernels |
| `GEN_CHUNK` / `--chunk` | 32 | Positive operations or pivots per generated helper function |
| `GEN_FUSION` / `--fusion` | 1 | Positive operations fused into one SIMD cell loop, bounded by a helper's size |
| `GEN_OPTIONS=--split-lu` | off | Separate cell loop for each LU update instead of fusing a pivot |
| `GEN_OPTIONS=--no-cache` | off | Recompute rates even when all inputs are unchanged |
| `GEN_OPTIONS=--generic-stages` | off | Dynamic stage loops and separate final-state/error loops |
| `GEN_OPTIONS=--joint-initial` | off | Combine the initial forcing and Jacobian pass, preserving each entry's accumulation order |
| `CHEM_GROUP_SIZE` | build default | Positive runtime group-size override |
| `--group-size` | environment/build default | Scenario-specific runtime group override |
| `--threads` | OpenMP maximum | Positive thread count, bounded by `INT_MAX` |
| `OMP_SCHEDULE` | static | OpenMP runtime schedule and optional chunk |
| `OPTFLAGS` | `-O3 -DNDEBUG` | User/compiler experiment flags; no extra optimization flags are added |
| `OMP_PREFIX` | `brew --prefix libomp` on macOS Clang | OpenMP include/library installation |

Combine generator switches in `GEN_OPTIONS`, for example
`make GEN_WIDTH=0 GEN_CHUNK=17 GEN_FUSION=3 test-specialized`.
All positive sizes must fit the available address space/memory; allocation failure
is reported normally. No size changes the chemistry. Group sizes affect adaptive
error averaging, so scenario comparisons must use the same grouping as MICM.
The fixed-stride kernels dispatch mismatched/tail groups to generic kernels;
variable-stride generation specializes every group size.

No SIMD length or thread count is hard-coded into the kernels. The C compiler
chooses SIMD instructions from ordinary C and `omp simd`. Helper chunk sizes,
loop fusion, cell groups and scheduling are all exposed above.
