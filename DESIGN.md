# Architecture and extension guide

## Model construction and execution

`include/chem.h` exposes a descriptor interface, an opaque `ck_model`, and a
public `ck_state` containing species-major arrays. All chemistry units are
mol/m3, K, Pa and seconds. Species and reactions are arbitrary caller-provided
arrays. A reaction retains each reactant occurrence and each product/yield
entry, including duplicates. Parameterized species contribute to the rate but
are excluded from integrated state variables.

`src/model.c` binds descriptors once: it validates species references, assigns
variable/custom-parameter indices, stably sorts rate types in reference order,
computes the diagonal Markowitz permutation, builds the sparse fill pattern,
and prepares derivative and LU operation schedules. Forcing and Jacobian loops
in `src/process.c` use those schedules for every mechanism. There is no TS1
branch or handwritten mechanism-specific kernel. Construction uses a dense
symbolic pattern; timed integration operates on compact sparse values.

The model copies its descriptor header, species descriptors, rates, bound
reactions/constraints and numeric schedules. Names, property tables, rate/species callback contexts, original
reaction descriptor arrays and initial-default arrays are immutable borrowed
storage and must outlive the model and all states. Generated descriptors have
static lifetime. `ck_model_free` releases owned storage. A model is read-only
while solving; change lambda functions only before concurrent use.

`src/state.c` owns conditions, concentrations, calculated rates, custom
parameters, tolerances and solver workspaces. Allocate with `ck_state_create`,
release with `ck_state_free`; clones have independent numeric storage. A model
may be shared across threads, with one state per concurrently solved group.
Each array is `[variable * cells + cell]`. Custom parameters and rate constants
use the same convention. `ck_update_rates` must be called after environmental
or custom-rate inputs change. Solves hold those calculated rates constant.

`ck_forcing` and `ck_jacobian` clear their destinations and include chemical and
external processes. Jacobians store **negative** derivatives. The accumulation
variants `ck_add_forcing` and `ck_subtract_jacobian` handle chemical reactions
without clearing pre-existing data. Algebraic rows exclude kinetic terms;
constraint callbacks then replace their forcing and supply constraint
Jacobian terms. The symbolic pattern retains external process entries even
when a constraint controls that row, so finalizers can bind every declaration.

`src/constraints.c` supports linear conservation and nonlinear equilibrium
constraints, including Van't Hoff temperature dependence, non-unit exponents,
clipped negative concentrations and correct zero-concentration derivatives.
Initialization solves only algebraic corrections, uses the reference's
weighted correction criterion and optional backtracking, and restores the
caller's original state on failure. Rosenbrock uses a zero mass diagonal for
algebraic rows. Its error norm includes every variable and substitutes the
algebraic step change for the embedded error estimate. Public post-solve
clamping affects differential variables only.

## Solvers and numerical order

`src/parameters.c` retains the two-, three-, four-, four-stage DAE and six-stage
DAE Rosenbrock tables. `src/solver.c` performs the original stage formulas,
adaptive error scaling, rejection handling and step-size limits. Defaults are
rtol `1e-6`, atol `1e-3`, initial internal step `1e-6 * DT`, minimum step
`1e-15 * DT`, and the reference method-specific coefficients. Parameters can
be overridden per call. `ck_solve_raw` exposes the underlying integration;
`ck_solve` includes external initialization and public post-solve handling.
`src/backward_euler.c` retains Newton convergence, positivity, failure step
reductions and accepted-step growth. Both solvers return explicit status and
statistics instead of exceptions.

Both LU algorithms are unpivoted, as in MICM. Mozart uses below-diagonal
reciprocal multiplication and scheduled updates. Doolittle uses the original
upper-row/lower-column order. Packed and separate-L/U interfaces share these
algorithms. Reaction, derivative, stage and reduction ordering is preserved;
there is no fast-math reassociation. SIMD is confined to independent cell lanes
using `restrict` and `#pragma omp simd`. All floating-point storage is `double`.

To add a Rosenbrock parameter set, add a table initializer in
`src/parameters.c`, select a new method ID in `ck_rosenbrock_parameters_init`,
and add numerical table/trajectory checks. For methods of at most six stages,
the integrator needs no change. More stages require extending the fixed stage
and triangular coefficient capacities in `chem.h`, `internal.h`, `state.c`
and solver validation. The generic stage loop already uses `parameters.stages`.

To replace the step-size factor, set `parameters.step_factor` and
`controller_context`, or change `ck_next_step_factor`. The callback receives
the parameters and normalized error and returns a multiplicative factor.
Acceptance, rejection recovery and endpoint bounds remain explicit in
`ck_solve_raw`; change that loop when experimenting with those rules too.

To replace linear algebra, provide `ck_linear_operations.factor`, `.solve` and
`.context`, and set the Rosenbrock or backward Euler parameter's `linear`
pointer. Callbacks receive the model, cell count, packed matrix and RHS.
`ck_jacobian_size`/`ck_jacobian_index` describe storage; callers can convert to
another representation inside their implementation. Null operations retain
the built-in algorithms. Callback dispatch is covered by the solver tests.

To add a rate type, extend `ck_rate_type`/`ck_rate` in `chem.h`, its defaults,
preparation and evaluation in `src/rates.c`, and the stable priority mapping
in `src/model.c`. Add schema/default/unit mappings in `tools/configure.py` if
it is configurable, followed by an independent formula test and a
configuration-vs-C comparison. `CK_LAMBDA` and `ck_model_set_lambda` already
allow C callbacks without a new rate enum. `src/state.c` multiplies
parameterized reactants after rate evaluation.

External processes use `ck_external_model`: declare added species, custom
parameters, algebraic species, named Jacobian entries and callbacks. Finalize
once after names and sparse indices are bound. Update per-state parameters,
accumulate forcing/negative Jacobians, and optionally supply constraint
initialization/residual/Jacobian callbacks. Supply a clone/destructor pair for
owned per-model context, otherwise context is borrowed. Context shared across
states must be immutable while solving; mutable cell data belongs in the
state. Test fixtures demonstrate process-only, constraint-only, composed,
temperature-dependent and aerosol models.

## Configuration to C

`tools/configure.py` implements CAMP semantics; PyYAML and Python's JSON module
provide serialization only. It resolves relative `camp-files` includes, rejects
cycles, collects `CHEM_SPEC` and `MECHANISM` entries, validates references and
numeric fields, expands branched reactions, converts rate units, and emits
static C arrays. The generated mechanism is passed through the same
`ck_model_create` binding path used by C-built mechanisms. Python never runs
inside the library or benchmark.

The format follows the three committed tutorial configurations. Default rate
concentration units are molecule/cm3, with reaction-order conversions based on
`6.02214076e17` molecules/cm3 per mol/m3. TS1 declares `mol m-3` because the C++
source already converted its constants. Supported configurable types include
Arrhenius, Troe, ternary activation/JPL19, branched, tunneling, Taylor,
reversible, surface and user-defined rates, plus photolysis, loss and emission
aliases. C function callbacks remain C descriptors.

Documented extensions are `concentration units`, species
`parameterization: air-density`, and `initial state` (temperature, pressure,
concentrations and custom parameters). Initial concentrations are SI mol/m3;
custom values use the corresponding parameter's units. Ordered reactant lists
and product `{species, yield}` lists preserve repeated terms and update order.
`tools/import_ts1.py` translates only source text, checks the 210 species,
547 reactions and 209 initial concentrations, and writes committed JSON.
Normal builds have no MICM source dependency. `tools/import_parameters.py` and
`tools/import_testdata.py` retain one-time table/fixture provenance.

## Storage, parallelism and tuning

The benchmark partitions the input into contiguous independent cell groups.
Each group has its own adaptive step history, error norm and workspace; its
cells advance together. Default groups shrink to at most ceil(cells/threads) to use available threads
for small inputs; explicit group sizes are used as requested. The final partial
group has its actual size, so unused
lanes never affect error control. OpenMP distributes groups, and contiguous
cell lanes expose portable vectorization. There is no cross-cell reduction
except each group's ordered error norm, and no shared mutable solver state.
The warm-up creates the thread pool before timing. Timing excludes CSV parsing,
model/state construction and final output. Statistics are summed over timed
group solves, so their totals depend on group count.

| Parameter | Default | How to change |
|---|---|---|
| Cell group size / contiguous SIMD span | 32 | `make GROUP_SIZE=N cpu`, `CK_GROUP_SIZE`, or `--group-size N` |
| Threads | OpenMP maximum, capped by groups | `OMP_NUM_THREADS` or `--threads N` |
| Group scheduling | static | `OMP_SCHEDULE`, using standard OpenMP runtime syntax |
| Utility matrix group width | caller-specified; 0 means row-major | `ck_matrix_create`/`ck_sparse_create` argument |
| LU algorithm | Mozart | solver parameter `lu`, or replacement callbacks |

Valid group sizes and widths are positive (utility width 0 selects row-major).
Changing them changes rounding and adaptive history, not chemistry or solver
formulas. No target-specific vector length, CPU affinity, instruction,
optimization flag or OS API is used.

`chem_matrix.h` provides element-size-generic row-major/AoSoA dense matrices,
CSR/CSC sparse matrices and padded groups. `ck_view` represents a column,
sparse element over blocks, or vector with explicit mutability. Group views
precompute contiguous bases. Iteration, fill/copy, reduction and prepared array
callbacks validate shape/order/writability at runtime. `chem_verify.h` provides
central/one-sided finite differences, negative-Jacobian comparisons and missing
sparsity detection.

## Specialized code generation

`tools/specialize.py` imports the validated CAMP reader. `bind` applies stable
rate-type sorting, removes parameterized species from integrated variables,
computes the reference Markowitz permutation, performs symbolic fill and creates
fixed-index derivative, factorization and triangular-solve schedules. Repeated
reactants/products remain separate and preserve their update order. No TS1 name,
reaction count, initial concentration or photolysis assumption enters binding.

`rate_body` emits fixed-parameter C formulas for every supported configurable
rate, keeping temperature, pressure, independently supplied air density and all
custom parameters as inputs. Parameterized reactants multiply the calculated
rate for every occurrence. Prepared branch constants use the C model's setup
calculations to avoid introducing Python/libm rounding differences. All other
rate constants are calculated by generated C; Python never runs at integration
time. Initial concentrations, initial custom values and tolerances are defaults,
never constants in the numerical kernels.

`Emitter.function` batches scalar arithmetic units inside independent-cell
OpenMP SIMD loops, splits helpers at a configurable boundary, and can emit both
constant-stride and dynamic-stride versions. Full groups use the preferred width;
any other width uses the dynamic version with its actual cell count. A separate
C translation unit dispatches helper arrays to bound compiler inlining cost.
Even constant-stride helpers retain runtime loop bounds so the compiler does
not expand the entire cell group into instructions. Optional tiling changes
only the independent-cell traversal, while retaining the full group stride.
The schedules and arithmetic are generated; helper dispatch has no sparsity
lookup or chemistry decisions.

`src/specialized_solver.h` is included into each generated translation unit.
It explicitly performs MICM's three stages, copying the shared forcing before
adding stage terms, and fuses independent array passes while retaining each
element's arithmetic order. It recalculates the Jacobian after rejection,
refactors for every attempted step and applies the public nonnegative clamp.
The initial forcing and Jacobian may be computed together, grouping derivative
jobs by reaction. For any matrix entry, its column and the order of contributing
reactions remain fixed; independent-column interleaving does not change its
arithmetic. Within each LU pivot, updates target distinct entries and read only
the completed pivot row/column, so row and column traversal orders perform the
same arithmetic per entry. The ordered error reduction is shared with the
generic solver. The specialized API selects the default three-stage coefficients and controller; the full generic
API continues to support other methods, custom controllers and constraints.

| Build parameter | Default | Meaning |
|---|---:|---|
| `GROUP_SIZE` | 32 | Driver's default consecutive cell-group width |
| `SPECIALIZE_WIDTH` | `GROUP_SIZE` | Preferred constant stride; 0 emits only dynamic paths |
| `SPECIALIZE_BATCH` | 8 | Chemistry reactions/derivative jobs per SIMD cell loop |
| `SPECIALIZE_ALGEBRA_BATCH` | 8 | Factorization/triangular operations per cell loop |
| `SPECIALIZE_SPLIT` | 32 | Arithmetic units per helper before starting another helper |
| `SPECIALIZE_TILE` | 0 | Cells per kernel tile; 0 uses the full group |
| `SPECIALIZE_FUSE` | 0 | Combine the initial forcing and Jacobian pass (0 or 1) |
| `SPECIALIZE_JOINT_BATCH` | 1 | Combined reaction/derivative units per cell loop |
| `SPECIALIZE_LU_ORDER` | column | Independent updates within each pivot, row or column |

All batching/splitting values must be positive; preferred width and tile may be zero.
Every setting preserves arithmetic dependencies and supports all positive
runtime group widths. Group width is also selectable with `--group-size` or
`CK_GROUP_SIZE`; threads and scheduling use the existing OpenMP settings.
Generator options `--batch`, `--algebra-batch`, `--split`, `--width`, `--tile`,
`--fuse`, `--joint-batch` and `--lu-order` expose the same controls for standalone
generation. Make tracks these values separately from compiler flags and regenerates/rebuilds when they change. Configuration
includes are Make dependencies for the bundled mechanisms; a new mechanism's
rule must likewise list its included files. `make clean` removes every emitted
source. No compiler optimization flags beyond `OPTFLAGS` are added.

`src/bench_common.h` shares CSV validation, state grouping and final output
between both benchmarks and the scenario driver. The scenario driver retains
original custom values independently of the current photolysis-scaled values.
Interpolation matches MICM at knots and outside table endpoints. At snapshot
boundaries it completes the parallel solve block, pauses timing, copies states
and resumes; file output follows integration.
