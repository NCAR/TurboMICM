# Design

## Objects and ownership

`include/chem.h` is the chemistry/solver API; `include/chem_array.h` is the
standalone dense/sparse array and view API. There is no process-global chemistry
state. Return values are `CH_OK` or a negative `ch_status`; constructors return
an owned pointer or NULL, and builder insertion functions return a nonnegative
ID or an error. `ch_status_string` explains failures. Free functions accept NULL.

A model owns species metadata, custom parameter defaults, reaction terms,
constraints, symbolic sparsity, and execution plans. Builder species IDs remain
valid until finalization; numerical state indices are assigned by finalization.
Use name lookup functions after finalization rather than retaining a builder
ID as a state offset. `ch_model_clone_builder` deep-copies an unfinished builder.
After a failed finalization, free the model and build a fresh one.

A finalized model can be shared by independent states. Do not change rates,
metadata, or callback contexts while any thread is using that model. Callback
contexts are borrowed, including across builder and state copies; applications
own their lifetime and synchronization. Strings returned by getters are
borrowed. Ordinary metadata and reaction arrays are copied on input.

A state owns conditions, concentrations, custom parameters, rates, tolerances,
and all integration workspaces. The model must outlive its states. State cloning
copies physical data and settings, with fresh work arrays. A state cannot be used
by two concurrent solves. Concentration buffers can swap after an accepted step:
reacquire `ch_state_values()` after each solve rather than retaining its pointer.

## Layout and execution

Chemistry arrays are species/reaction/sparse-entry major, with the cells of one
group contiguous: `values[species * cells + cell]`. A state is one independent
adaptive group, including a possibly short final group. There are no inactive
padding lanes in solver states. Small inner loops use `restrict` and OpenMP SIMD;
there are no architecture-specific kernels, intrinsics, or precision changes.

`src/model.c` finalizes a mechanism once. It separates parameterized reactants,
optionally performs the stable diagonal Markowitz reordering, derives the
negative Jacobian's sparse pattern, computes LU fill, and records the indices
needed by each derivative and factorization operation. Kinetic accumulation
preserves repeated reactants and ordered product contributions. No name lookups,
sparsity searches, or allocations occur in the normal ODE solve loop.

`src/state.c` computes rates, forcing, and the negative analytical Jacobian.
Rates include the product of parameterized-reactant values. They are refreshed
by `ch_update_rates`; call that function when conditions, custom parameters,
rate definitions, or external parameterizations change. Native TS1 computes
rates before warm-up, since these inputs remain fixed during the benchmark.

`src/sparse.c` supports Mozart outer-product and Doolittle row-wise elimination,
with identical symbolic sparsity and unit lower diagonals. The public packed
matrix stores actual U diagonal values. Solver states retain the pivot inverses
already computed during factorization for reuse by each stage's back solve.
Standalone `ch_lu_solve` accepts an ordinary packed LU and divides by its diagonal.
Out-of-place factorization copies the input into the caller's LU buffer first.
There is no pivoting or substitute linear algebra package.

`src/array.c` implements the utility matrix equivalents used by the CPU cases:
dense and sparse row/column ordering, ordinary and grouped layouts, real and
integer data, views, strict padding operations, mapping, reduction, copying,
printing, and reusable callback dispatch. Group width zero selects ordinary
storage; any positive width selects an interleaved layout. Logical lengths and
storage capacities are distinct, especially in the final group.

## Numerical methods

`src/parameters.c` contains all five MICM coefficient sets: two-, three-, and
four-stage Rosenbrock, four-stage DAE (selector 43), and six-stage DAE. The
solver engine receives a `ch_rosenbrock_parameters` structure, not a method ID.
Its stage equations, function reuse, weighted RMS error, rejection rules, and
step-size limits follow MICM. The default controller evaluates
`min(factor_max, max(factor_min, safety / error^(1/order)))`. Summation of the
normalized error is kept ordered; it is not an OpenMP reduction.

For default Rosenbrock parameters, minimum/start/maximum steps are derived from
the requested interval as `1e-15 * dt`, `1e-6 * dt`, and `dt`. State defaults are
relative tolerance `1e-6` and species absolute tolerance `1e-3`. TS1 uses those
defaults and the unmodified three-stage table. The benchmark never overrides
tolerances, the error estimate, or the controller.

`src/backward_euler.c` retains MICM's Newton iterations, positivity treatment,
convergence criteria, failed-step reductions, and step growth after two successes.
Exhausting reductions returns `CH_UNCONVERGED` with the last trial and reported
elapsed time, corresponding to `AcceptingUnconvergedIntegration`. NaN and infinity
are reported separately. Rosenbrock's `CH_MAX_STEPS` is a separate status.

`src/constraints.c` handles linear and equilibrium constraints and external
constraint callbacks. Algebraic rows replace kinetic rows, and the mass-matrix
shift applies only to differential rows. Consistent initialization uses weighted
Newton iterations and backtracking, preserving differential variables and
rolling back on failure. Each solve reinitializes constraints. Positivity
postprocessing applies to differential variables. Backward Euler remains an ODE
method, as in MICM.

Statistics count work per independent group, not per cell. Thus changing group
size changes summed statistics even when every cell follows the same nine
internal steps. `ch_stats_add` sums counters; its `final_time` is summed too.
For a single solve, `final_time` reports progress and permits resuming after the
maximum-step limit, as in the analytical test drivers.

## Configuration to native code

`tools/configure.py` has two stages: `read_config` normalizes a mechanism into a
Python dataclass; `generate` emits ordinary C constructor and initializer calls.
PyYAML handles YAML syntax, and the standard JSON module handles JSON syntax.
The chemistry schema, inclusion, validation, defaulting, branching, and units
are implemented here without MICM or MUSICA code dependencies.

Relative `camp-files` paths resolve against the including file; cycles fail.
`CHEM_SPEC` records define species and optional surface properties; `MECHANISM`
records define reactions. Unknown species, unsupported rate types, duplicate
species, invalid quantities, missing surface properties, and nonfinite numeric
values fail generation. Mechanisms default to molecule/cm3 rate units and are
converted to native mol/m3 units according to reaction order and rate type.
Both branch channels become separate reactions. Photolysis, emission, and
first-order loss use named custom rate parameters.

The generated C contains data and calls to the same model API used by hand-built
tests. It finalizes with reordering enabled and unused species permitted.
Generated files live under `build/generated` and are compiled with `CC` and
`OPTFLAGS`. No generated mechanism-specific numerical solver exists. Changes to
the reader or configuration trigger regeneration.

`tools/translate_ts1.py` parses the supplied TS1 data declarations, verifies 210
species and 547 reactions, and emits the committed `configs/ts1.json`. It reads
initial concentrations, parameters, temperature, and pressure from the source.
TS1 has 209 evolving species because `M` is parameterized by air density. The
format extensions are specified in [DECISIONS.md](DECISIONS.md). All three
tutorial configurations generate C and are compared with hand-built mechanisms
over ten solves in three differing cells.

## Tuning parameters

All parameters below leave each group's equations and convergence tests intact.
Group size may affect adaptive step choices through its RMS error; correctness
does not depend on a fixed size, SIMD width, or worker identity.

| Parameter | Default | Selection and valid values |
| --- | --- | --- |
| Benchmark group size | 32 cells | `make GROUP_SIZE=N` changes the compiled default; `CH_GROUP_SIZE=N` overrides it; `--group-size N` overrides both. Any positive representable size that fits memory; capped at the cell count. |
| Benchmark threads | OpenMP maximum | `OMP_NUM_THREADS` sets the runtime default; `--threads N` takes precedence. Positive integers through `INT_MAX`, capped at the number of groups. |
| Benchmark scheduling | `dynamic,1` | Standard `OMP_SCHEDULE`, e.g. `static`, `dynamic,4`, or `guided,1`. A group is assigned as a whole, including all its timed steps. |
| `ch_solve_many` threads | OpenMP maximum | Function argument; nonpositive means default, positive is capped at state count. Its schedule is static. |
| Utility-array group width | Caller supplied | Zero for ordinary row/block storage; any positive width for grouped storage. Independent of solver group size. |
| LU algorithm | Mozart (0) | `ch_state_set_lu(state,0)` or Doolittle (1). A callback pair may replace both. |
| Species reordering | Caller supplied | `ch_model_finalize(model,reorder,ignore_unused)`; generated constructors enable reordering. |

No hard-coded vector width, processor affinity, architecture flags, or compiler
optimization experiments are used. Larger groups reduce dispatch overhead but
increase working-set size. The measured 16/32/64/128/256 sweep favored 32; this is
a configurable default, not a numerical requirement. The primary runtime uses
all cores and dynamic group assignment to allow different group costs.

## Extension points

**Add a Rosenbrock table.** Create a `ch_rosenbrock_parameters`, initialize the
common controls (usually by copying a factory result), then fill `stages`,
`new_function`, `a`, `c`, `m`, `e`, `alpha`, `gamma`, and `estimator_order`. Lower
triangular stage coefficients use offset `stage*(stage-1)/2 + previous_stage`.
Pass the structure directly to `ch_solve`. Adding a named factory selector only
requires another case in `src/parameters.c`; the engine does not switch on it.
The current workspace supports up to six stages. A larger method requires
expanding the parameter arrays, stage bound, and state stage allocation together.

**Replace the step-size controller.** Assign `parameters.controller` and
`controller_context`. The callback receives the normalized error and parameter
table and returns a positive finite step multiplier. Acceptance, rejection
history, and min/max bounds stay in the engine. Change `src/rosenbrock.c` when
experimenting with acceptance logic itself. Contexts used across threads must be
immutable or synchronized; per-state counters need per-state contexts.

**Replace the linear solver.** Call `ch_state_set_linear_solver` with factor and
solve callbacks and a context. Factorization receives the sparse symbolic plan,
cell count, mutable packed matrix, and a scratch array of `cells` doubles. The
matrix is `alpha * mass - Jacobian`; factorization may overwrite it. A nonzero
factor status stops integration. Solve receives the factorized matrix and
overwrites the RHS. Public CSR getters expose dimension, row offsets, column
indices, diagonal indices, and fill size. Additional backend storage belongs in
the context. `tests/native_cases.c` demonstrates installing Doolittle through
this interface while replacing the controller, checks invocation counts against
solver statistics, and compares the solution with the default backend.

**Add a rate type.** For a runtime function, use `CH_LAMBDA` with
`ch_rate_function` and context; named custom parameters are available to it. For
a built-in formula, append an enum value, define defaults and evaluation in
`src/rates.c`, and extend constructor validation and rate-group ordering in
`src/model.c` / `src/inspect.c`. Add the schema fields and unit conversion to
`tools/configure.py` if configuration should support it. Add numerical rate
checks and a small mechanism test. For a new external process family, use
`ch_external` to declare sparsity, update parameters, and accumulate forcing
and negative Jacobian contributions without changing mass-action kinetics.

The `ch_rate.p` fields are:

| Rate | Parameter slots, in order |
| --- | --- |
| Arrhenius | A, B, C, D, E |
| Troe / ternary activation | k0_A, k0_B, k0_C, kinf_A, kinf_B, kinf_C, Fc, N |
| Ternary JPL19 | k0_A/B/C/D, kinf_A/B/C/D, kint_A/B/C/D, Fc, N |
| Tunneling | A, B, C |
| Branched | X, Y, a0, n, branch (0 alkoxy / 1 nitrate); slots 5–6 prepared internally |
| Surface | diffusion, mean-speed coefficient, reaction probability; first two derived from species properties |
| User defined | scaling factor; `custom_index` selects its parameter |
| Taylor | Arrhenius A/B/C/D/E followed by 1–16 polynomial coefficients |
| Reversible | A, C, reverse rate |

Initialize rates with `ch_rate_default`, then set the desired fields. Model
insertion prepares derived constants. If changing branched fields after
finalization, call `ch_rate_prepare` before updating rates. Surface custom
parameters are adjacent radius and particle-number-density entries.
