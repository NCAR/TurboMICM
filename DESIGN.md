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
