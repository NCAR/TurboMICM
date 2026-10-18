# Architecture

The public interface is `include/chem.h`. A `chem_mechanism` is an immutable
array description of species, reactions, rates, and initial conditions. It can
be initialized in C or emitted by the configuration compiler. A `chem_model`
compiles this description into species maps, custom-parameter slots, Jacobian
jobs, and a symbolic sparse LU plan. The description must outlive the model.
Each `chem_state` owns its concentrations, rate parameters, conditions,
tolerances, and work arrays. Models can be shared between threads; states are
independent. Lambda callbacks and their contexts remain caller-owned.

The numerical arrays store a contiguous column of cells for each species,
reaction, or sparse element. The solver acts on a group of cells with a single
step-size controller, as in the reference. Species use diagonal Markowitz
ordering. The LU plan includes symbolic fill and an explicit operation sequence
for either Mozart elimination or Doolittle elimination. The solve uses forward
then backward substitution. Scalar operations keep their reference order;
OpenMP SIMD directives operate across independent cells.

## Configuration compilation

`tools/configure.py` reads JSON or YAML, resolves `camp-files` relative to their
containing file, checks species references, expands branches and multiplicities,
and converts rate units. `tools/generate.py` emits the immutable C description
and chemistry/Jacobian kernels. These kernels are generated from reaction terms
and the generic symbolic sparsity algorithm. There is no hand-written
mechanism-specific kernel. `tools/translate_ts1.py` transcribes the supplied TS1
C++ declarations without compiling or executing them; `configs/ts1.json` is its
retained output. All emitted C is compiled by Make with `CC` and `OPTFLAGS`.

The table-driven C evaluator implements the same operations and is selected by
the benchmark's `--generic` option. Generated kernels currently support kinetic
mechanisms; models containing constraints use the generic evaluator.

## Arrays and views

`chem_array.h` supplies owned dense matrices and vectors in native double, int,
or size_t storage. A vector length of zero selects row-major storage; any
positive length selects padded groups of contiguous column lanes. The solver's
fixed-cell SoA arrays remain its compact numerical workspace. `chem_sparse.h`
adds CSR/CSC symbolic structure and stores the nonzero columns in a `chem_array`.
Sparse builders deduplicate elements and report out-of-range coordinates.

Array callbacks record argument types, layouts, and column/nonzero counts at
creation. Row/block counts are checked on invocation and may differ from those
used at creation. All arguments must have equal logical row/block counts and
compatible group widths. The callback receives group descriptors; column and
sparse-element spans resolve the base pointer once. Raw descriptors can also be
bound to a group. Const spans reject writes. Mapping, filling, copying, and
reductions work across dense columns, sparse elements, vectors, and caller-owned
temporary spans. Strict reduction visits only active lanes. Non-strict operations
include initialized trailing padding, matching the reference policy semantics.

The utility arrays expose native storage for exact integer/index operations.
The numerical map/get/set conveniences convert values to double, so callers
handling indices above the exact-double integer range should use native storage.

## Numerical extension points

`chem_ros_parameters` contains the stage count, flattened lower-triangular A/C
coefficients, solution and error coefficients, stage flags, and controller
constants. The five reference parameter sets are in `src/parameters.c`. To add
a parameter set, return another initialized structure; the stage loop does not
need a new method-specific branch. Six stages are currently the maximum, matching
the reference workspace contract.

`chem_ros_solve` accepts a step-controller function and context. The default
`chem_ros_controller` computes the reference bounded safety-factor proposal.
Acceptance, rejection history, minimum/maximum steps, and failure detection
remain explicit in `src/rosenbrock.c`.

`chem_be_parameters` contains backward Euler's iteration limit and five reduction
factors. Its Newton iterations, positivity clamp, convergence checks, and
unconverged/nonfinite result handling are in `src/backward_euler.c`.

`chem_lu_algorithm` selects Mozart or Doolittle operation generation.
`chem_factor.h` adapts the same operation plans to grouped CSR/CSC arrays, with
combined or separate L/U storage. Separate factors have a unit lower diagonal.
`chem_state_set_linear_backend` replaces factorization and solution per state;
Rosenbrock, backward Euler, and constraint initialization all use these hooks.
Callbacks and their context remain caller-owned; copying a state copies this
binding, so independent backend workspace requires rebinding the copied state.

A rate type consists of a tag, a parameter structure in `chem_rate`, defaults,
optional derived-field preparation, and a calculation branch in `src/rates.c`.
The configuration emitter adds its field names and unit conversions separately.
Custom functions can use the lambda-rate callback without adding a tag.

## Performance parameters

| Parameter | Default | Selection |
| --- | --- | --- |
| Cell-group size | 32 | `--group-size G`; build macro `CHEM_GROUP_SIZE` |
| Thread count | OpenMP processor count | `--threads N` |
| Group scheduling | static contiguous blocks | `--schedule static`, `dynamic`, or `guided` |
| Schedule chunk | OpenMP policy default | `--chunk K`, positive number of groups |
| SIMD width | compiler-selected | portable `omp simd`, no hardware-specific width |
| Chemistry evaluator | automatically generated | `--generic` for the generic table evaluator |
| LU method | Mozart in benchmark | library `chem_lu_algorithm` |

There are no target-specific flags, intrinsics, vector extensions, or alternative
floating-point modes. Group size affects work distribution and the common error
norm, not the method or requested tolerances. Every positive group size and
thread count remains independent of the numerical formulas. Correctness checks
include one-cell groups, irregular tails, groups larger than the domain, and all
three scheduling policies.

## Owned species metadata

`chem_system.h` supplies typed properties, owned species metadata, phases, and
systems. Create/clone destinations are zero initialized; replacement constructs
a full copy before releasing the previous value. Property strings and phase
species are deep copied. The compact `chem_species` descriptor borrows its name
from metadata and copies affine parameter coefficients by value. Thus descriptors
can be assembled into mechanisms without retaining a pointer to a metadata
object for parameter evaluation; the name and any callback context must still
outlive the model. Phase names optionally prefix nonparameterized species names;
system names omit the gas prefix and accept a validated permutation.

## State assignment and process ownership

`chem_state_copy` duplicates values and workspace; `chem_state_move` transfers
ownership and clears the source. Both accept zero-initialized destinations or
replace existing states. `chem_state_assign_*` requires exactly one value per
cell, matching reference state assignment. Batch assignments validate all inputs
before writing. The lower-level `chem_state_set_*` functions explicitly support
broadcasting one value to all cells. Column access returns borrowed mutable or
const pointers; bulk custom input is row-major and is copied into SoA storage.

`chem_process.h` owns phase information, reactant/product lists, and rate labels.
`chem_process_model_create` resolves process species by name against a copied
system, stably sorts reaction types, and compiles an owned mechanism. Its system
and all input processes can be destroyed after creation. Low-level
`chem_model_create` preserves the supplied reaction order, so callers that require
reference builder ordering use the owned builder or `chem_reactions_sort`.
Kinetic sparsity and used-species queries exclude solver fill and extra diagonals.

## Owned solver and concurrent states

`chem_solver_create` combines an owned process model with Rosenbrock or backward
Euler options and one-time setup callbacks. The normal `chem_solver_solve` entry
point is const: threads may share a solver while owning separate states and
backend workspaces. Explicit parameter-override entry points update solver
options and require exclusive access. Setup descriptors are borrowed only during
creation. Preferred vector length describes one group; states can span multiple
groups. Post-solve clamping applies to differential variables only, matching the
reference wrapper; algebraic variables remain governed by their constraints.

## External models

`chem_external.h` describes extra species, named parameters, kinetic sparse
entries, algebraic rows and constraint dependencies. Composition merges these
before Markowitz ordering and LU planning. Each context is either an explicit
borrow, a byte-copied value, or owned through clone/destroy callbacks. Finalize
callbacks bind names to reordered indices independently for each model; update
callbacks run when state parameters are refreshed. Kinetic and constraint
callbacks add to SoA destinations. Kinetic evaluation preserves algebraic rows,
including their incoming values. Workspace belongs to each state, so a model can
be shared between threads without callback workspace races.

The benchmark also exposes `--schedule static|dynamic|guided` and `--chunk K`,
implemented with OpenMP runtime scheduling. These affect which thread handles a
complete cell group; they do not share a group's solver workspace between
threads or reorder arithmetic inside it. A zero internal chunk selects OpenMP's
default for that scheduling policy; an explicit chunk must be positive.

To replace the integrator, add a function taking a `chem_state`, interval and
parameter object and returning `chem_result`, then select it in
`chem_solver_solve` (or call it directly). Keep its temporary data in the state
or an explicitly owned backend context. The state-level forcing/Jacobian,
constraint initialization, factor and solve entry points allow a new integrator
to reuse chemistry and algebra without changing mechanism descriptions. The
benchmark's integrator call is explicit in `src/bench.c`; experiments can replace
that call while retaining CSV initialization, cell grouping and timing.

# Tier 3 integrators and scenario interface

`include/chem_adaptive.h` adds an autonomous ODE interface. Create one
`chem_adaptive_workspace` per state, set the state's `rtol` and species `atol`,
then call `chem_adaptive_solve(state, dt, parameters, workspace)`. Destroy the
workspace before destroying the state/model. The workspace borrows the model
and checks the cell count. It owns a saved sparse negative Jacobian and twelve
species-vector slots. It does not allocate during integration. A workspace is
not shared between threads. Existing Rosenbrock, backward Euler, constraints,
and owned-solver APIs retain their original behavior.

The new methods evaluate every reaction and evolve every variable species.
There is no reaction reduction, species selection, trained approximation,
mechanism-specific threshold, or TS1-specific integrator code. The new API is
for autonomous differential equations: it explicitly rejects algebraic
constraints. Temperature and reaction parameters must remain fixed within a
call; the scenario runner enforces this and starts a new call at each forcing
update. External kinetic callbacks must represent the same autonomous ODE.

## Runtime selection

```
bin/ts1_scenario [--threads N] [--method NAME] [--option KEY=VALUE ...] SCENARIO
```

| Method | Formula and implementation |
| --- | --- |
| `micm` (default) | Original MICM-compatible three-stage Rosenbrock, original defaults |
| `micm2`, `micm4`, `micm4dae`, `micm6` | Other retained MICM parameter sets, original driver |
| `ros2`, `ros3`, `ros4`, `rodas3`, `rodas4` | Retained coefficient tables in the new adaptive driver; respectively 2/3/4/4/6 stages and order 2/3/4/3/4 |
| `extrap` | Fixed-Jacobian linearly implicit Euler extrapolation, order 2–4; optionally chooses order during a trial |
| `sdirk` | Stiffly accurate singly diagonally implicit Runge–Kutta, order 2 or 3, with simplified Newton iteration |

The default method and its tolerances are intentionally unchanged even though
its default absolute tolerance is unsuitable for several disclosed long runs.
Use an explicit configuration from OPERATING_POINTS.md for atmospheric runs.
The methods are library functions, not Python implementations behind a CLI.

## Formula details

The chemistry interface supplies `-J`. A factorization of `alpha I - J` is
formed in the state's original sparse storage; the unmodified `-J` stays in
separate workspace. All arithmetic is double precision, including error
norms, factors, and triangular solves. There is no mixed-precision or floating
point flag change.

**Rosenbrock.** The coefficient tables and stage equations are exactly those
already retained in `src/parameters.c`. The new path changes control and cache
policy, not their coefficients. ROS2 permits approximate Jacobians without
losing its formal second order [REFERENCES.md, source 3]. Reusing Jacobians in
higher-order ROW formulas is only an experimental approximation; no retained
order claim is made for that option. The declared operating points use fresh
Jacobians. Stage clipping is a separately selectable experiment and is off by
default. It is not required for output positivity.

**Euler extrapolation.** For each row with `q = 1,...,p`, start at the same
initial state and take `q` substeps of length `d = h/q`:

```
( I/d - J_saved ) delta = f(v)
v <- v + delta
```

Each row needs one factorization and `q` solves. The initial forcing is shared
by all rows and all retries at the same initial state. The remaining substeps
have new forcing evaluations. Neville extrapolation in `1/q` eliminates the
successive substep-size error terms. In the in-place table, descending `j`
updates `T[j-1] = T[j] + (T[j]-T[j-1]) / (q/j-1)`. The returned diagonal has
formal nonstiff order `p`; its difference from the preceding diagonal estimates
local error with exponent `p`. A successful fixed order-four trial needs four
factorizations, ten solves, and seven forcing evaluations (six if its initial
forcing was already cached). It uses every reaction in every evaluation.

`adaptive-order=1` returns the first acceptable diagonal of order at least two,
up to the selected maximum. Otherwise it always constructs the selected order.
This simple variable-order policy is an experiment, not an optimal-work order
selector; its poorer work–precision performance is reported in SCANS.md.

High nonstiff order is not a claim of uniform accuracy for arbitrarily stiff
nonlinear equations. On the linear scalar equation, the stability function is
a weighted sum of `(1-z/q)^(-q)`. Some extrapolated solutions become slightly
negative for very stiff decay. Their projection and associated conservation
change are explicitly measured; see ACCURACY.md and the stiff-decay test.

**SDIRK.** The two-stage formula has diagonal `g=1-1/sqrt(2)` and last row
`[1-g,g]`. Its embedded first-order weights are `[1,0]`. For three stages,
`g=0.435866521508459...` solves `6g^3-18g^2+9g-1=0`. The rows are
`[g,0,0]`, `[(1-g)/2,g,0]`, and
`[-1.5g^2+4g-1/4, 1.5g^2-5g+5/4, g]`. The last row supplies the solution weights.
The embedded weights are `[1-b,b,0]` with `b=(1/2-g)/((1-g)/2)`, giving order two.
The method families and diagonal root are documented in REFERENCES.md, source 4.

A stage with known base `b` solves `z=b+h*g*f(z)`. Simplified Newton solves
`(I/(h*g)-J_saved) delta = f(z)+(b-z)/(h*g)` and updates `z`. The iteration stops
when the maximum cell RMS correction in the integration tolerance scale is at
most `newton-tol`. It fails after `newton-iterations` or after detected loss of
contraction. Failure retries with a fresh Jacobian, then a smaller step. No
partial Newton stage is accepted. The final stage is the returned solution.

## Error, positivity, control, and caches

Each component has scale `atol[i] + rtol*max(abs(y[i]),abs(candidate[i]))`.
Before accepting a trial, divide the maximum of the embedded error magnitude
and the magnitude of a negative candidate by this scale. The default norm is
the largest per-cell RMS, so adding easy cells cannot dilute a difficult cell's
error. Other norms are controlled experiments. Acceptance requires a finite
norm at most one. Nonfinite candidates are rejected and never made finite by
clipping. Accepted negative values are set to zero and their number and sum are
recorded. This projection can change conservation; it is a bounded local
numerical correction, not an exact conservation algorithm.

The proportional controller proposes `h_new/h = safety * error^(-1/order)`,
with bounded decrease/growth and no growth immediately following rejection.
The PI option uses exponents `-0.7/order` for current error and `0.4/order` for
previous accepted error. Rejection with a stale Jacobian first refreshes it at
the same step size; rejection with a Jacobian already at the current state
reduces the step. Nonfinite trials reduce the step by 0.1. Every attempt counts
toward the step limit. Fatal backend errors propagate without being retried.

Only the next proposed step length persists across calls. Jacobians, factors,
initial forcing, and PI history are invalidated at a new call, including every
photolysis discontinuity. Within a call, `jacobian-age` bounds the number of
accepted steps sharing a Jacobian. `reuse-window` can retain an unchanged step
when its proposed factor is close to one; ROS/SDIRK can then reuse a
factorization if the Jacobian is also unchanged. Extrapolation factors each
row separately. No factor is reused after changing its diagonal shift.

## Parameters and valid domains

All speed-related choices are runtime parameters; no CPU-specific tuning flags
are used. Numbers outside these domains are rejected. A restrictive step limit
or impossible requested precision returns an explicit solver failure.

| Option | New-path default | Domain and effect |
| --- | --- | --- |
| `rtol` | `1e-4` | Finite positive relative local tolerance; MICM path default remains `1e-6` |
| `atol` | `1e-18` | Finite positive uniform absolute tolerance, mol m-3; original path retains species defaults (`1e-3` for TS1) |
| `group-size` | `32` | Any positive cell count; final group may be shorter |
| `generic` | `0` | `0/1`, generated/table-driven chemistry kernels |
| `rate-cache` | `1` | `0/1`, retain rate constants whose inputs are unchanged |
| `step-block` | `4` | Positive external steps per group visit, capped at each snapshot |
| `schedule` | `dynamic` | `static`, `dynamic`, or `guided` OpenMP group assignment |
| `chunk` | `0` | Zero selects OpenMP policy default; positive through `INT_MAX` groups |
| `max-steps` | `100000` | Positive attempt limit per group and external step; original path retains `1000` |
| `h-start` | `0` | Nonnegative seconds; zero selects `dt*1e-4` (original path: `dt*1e-6`) |
| `h-max` | `0` | Nonnegative seconds; zero selects the external interval |
| `warm-start` | `1` | `0/1`, retain proposed next step across calls |
| `safety` | `0.9` | Strictly between zero and one |
| `growth` | `5` | Finite value at least one, maximum step growth factor |
| `controller` | `0` | `0` proportional; `1` PI |
| `norm` | `1` | `0` global RMS; `1` maximum cell RMS; `2` component maximum |
| `jacobian-age` | `1` | Positive maximum accepted steps between updates |
| `reuse-window` | `0` | `[0,1)`, retain step if proposed factor differs from one by at most this value |
| `clip-stages` | `0` | `0/1`, optional ROS intermediate-stage projection |
| `order` | `3` | `2..4` for extrapolation; `2..3` for SDIRK |
| `adaptive-order` | `0` | `0/1`, extrapolation early acceptance at a lower order |
| `newton-iterations` | `8` | Positive SDIRK iteration limit per stage |
| `newton-tol` | `0.03` | `(0,1)`, SDIRK correction stopping tolerance |
| `--threads` | OpenMP processor count | Positive integer through `INT_MAX`; default uses all cores |

The compatibility methods use `rtol`, `atol`, `h-start`, `h-max`, `max-steps`,
`group-size`, `generic`, and all four scenario scheduling/rate-cache options;
options for the adaptive numerical driver do not change their original driver. `order` and nonlinear iteration parameters affect only the
methods described above. The original Tier 1 benchmark retains its original
static scheduling default and separate command-line scheduling/chunk controls.

## Scenario contract and accounting

`src/scenario.c` parses the documented key/value format, resolves paths relative
to the scenario file, interpolates the photolysis table, and counts CSV rows.
It preserves each cell's input custom parameters and updates all labels starting
with `j` from those original values. It updates affected rates at every external
step; no warm-up is performed. It validates returned states at every external
step, including those without snapshots.
The timed intervals include interpolation, parameter/rate updates, integration,
validation and OpenMP overhead. Input, setup and CSV serialization are untimed.
Both output formats use 17 significant digits and original cell order. Snapshots
include the final step even when not divisible by the snapshot interval.

The first complete rate calculation remains inside the timer. With `rate-cache=1`,
subsequent calls use a generic reaction dependency list: USER rates depend on
one custom slot and SURFACE rates on both slots. Conditions and unmarked custom
inputs are frozen by this scenario contract. Unknown lambda/species callbacks,
external composition or constraints conservatively select every rate. The API
requires full initialization before selected updates and a rebuilt list if the
model or varying inputs change. Every reaction remains in every forcing and
Jacobian evaluation. `rate-cache=0` recalculates all rate constants every step.
Unit tests cover every custom mask on unrelated mechanisms, both surface slots,
third bodies, unknown callbacks and invalid indices.

Step blocking advances each independent cell group for up to `step-block` external
steps before visiting another group. Each external boundary still updates forcing
and resets the integrator caches as specified above. Blocks stop at snapshot
boundaries, and photolysis interpolation remains timed. This reorders independent
trajectories without changing their arithmetic; it does not imply that a coupled
transport model can skip its own synchronization boundaries. Dynamic/guided/static
scheduling and blocking produce bit-identical states and operation counts in the
controlled scans and scenario tests.

Standard counters count group operations, as in Tier 1, not per-cell FLOPs.
Compare counts only with group sizes stated. `internal_steps` counts attempted
macro steps, including rejected ones; substeps are reflected in forcing and
linear-solve counts. Additional counters report Jacobian/factor reuses,
positivity/nonfinite rejections, projected values, summed projection magnitude,
Newton iterations/failures and the selected extrapolation orders. `rate_updates`
counts individual cell/reaction rate evaluations (unlike the group counters);
`rate_cache` records whether dependency caching was enabled. Projection
magnitude sums concentrations over all affected species and cells; it is not an
elemental mass or a reaction-rate modification.
The compatibility path retains MICM's rejection-counter convention, which
excludes rejections before the first accepted step of a call. Factorization and
function/solve counts are more direct measures of work across those drivers.

## Adding and verifying a method

Add a method enum and stage implementation to the adaptive API, using the
state's kinetic kernels and linear-backend hooks. Preserve the separate original
driver. Add a CLI name or option, an independent prototype, small-mechanism
agreement with operation counts, and scenario configurations in `scans/`.
`make test-tier3` runs prototype, analytic-order, conservation/projection,
stiff-decay, error-path, scenario I/O and generic-kernel checks. `make test-cpu`
continues to run the 680 original cases. Record new structural runs with
`tools/scan.py`, then rescore final operating points on all disclosed cells.
