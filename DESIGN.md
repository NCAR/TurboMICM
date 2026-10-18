# Design

## Scenario runner and added ODE methods

`bin/ts1_scenario [--threads N] [--method NAME] [--option KEY=VALUE ...] SCENARIO`
reads the supplied scenario format. Its default, `rosenbrock3`, uses the legacy
solver and unchanged parameters/tolerances. The other legacy names are
`rosenbrock2`, `rosenbrock4`, `rosenbrock4dae`, and `rosenbrock6dae`. They accept
only `rtol`, `atol`, `h_start`, `h_max`, `max_steps`, and `group_size`. Omitted
tolerances retain their original values. A scenario has no warmup.

Added methods use `ch_integrator`, implemented in `src/integrator.c`:

| Name | API selector | Formula | Stages / fresh RHS per trial |
| --- | ---: | --- | ---: |
| `ros2` | 2 | L-stable Rosenbrock 2(1) | 2 / 2 |
| `ros3` | 3 | L-stable Rosenbrock 3(2) | 3 / 2 |
| `ros4` | 4 | L-stable Rosenbrock 4(3) | 4 / 3 |
| `rodas3` | 43 | Stiffly accurate Rosenbrock 3(2) | 4 / 3 |
| `rodas4` | 6 | Stiffly accurate Rosenbrock 4(3) | 6 / 6 |
| `rodas31` | 431 | RODAS3.1 published 3(2) coefficients | 4 / 3 |
| `sdirk2` | 102 | L-stable SDIRK 2(1), chord iteration | 2 / iteration-dependent |

The first five reuse Tier 1 mathematical tables with a new adaptive engine.
RODAS3.1 is independently transcribed from the mathematical table cited in
[REFERENCES.md](REFERENCES.md). The API is mechanism independent and integrates
every reaction and evolving species. Added methods reject DAEs; the retained
solvers continue to support them.

Rosenbrock trials factor `A = I/(gamma*h) - J`, with the analytical Jacobian.
Stage right sides are `f(y + sum(a_ij*k_j)) + sum(c_ij*k_j/h)`; the candidate is
`y + sum(m_i*k_i)` and embedded error is `sum(e_i*k_i)`. One LU serves all
stages. The initial RHS is exactly reused on rejected trials because state and
rates have not changed. There is no approximate RHS reuse. Fresh-Jacobian runs
assemble directly into the factor workspace; only lagged runs allocate and
copy a preserved unfactored Jacobian.

SDIRK2 uses `g = 1 - 1/sqrt(2)` and solves
`Y1 = y + g*h*f(Y1)`, then `Y2 = y + (1-g)*h*f(Y1) + g*h*f(Y2)`.
Chord iterations share `I/(g*h) - J(y)` across both stages. Each correction
solves the scaled residual equation; its normalized size must be at most
`newton_tolerance`. Failure rejects the trial. The candidate is `Y2`; the
first-order embedded candidate is `y + (Y1-y)/g`. Store the error as
`Y2 - [y + (Y1-y)/g]`, with controller order 2, so adding a projection to `Y2`
adds the same signed amount to the error. A stiff scalar regression checks
the acceptance decision when that projection moves the error through one.
This inexpensive estimator is not L-stable; fewer stages do not guarantee lower
cost on stiff transients.

### Error, positivity, controllers and reuse

Species scales are `atol[i] + rtol*max(abs(y_i),abs(candidate_i))`.
`norm=rms` takes the RMS across species/cells, `cell` takes the largest cell's
species RMS, and `max` takes the largest component. Cell norms prevent dilution
of a difficult cell by easy cells. Species accumulate in order; SIMD operates
across cells. Infinity and NaN are explicitly rejected.

Candidates below `-positivity*atol[i]` reject the trial. Smaller negative
components are projected to zero and the correction is added to the embedded
error **before** acceptance. This is bounded projection, not a proof of a
positivity-preserving discretization or exact conservation. It guarantees
nonnegative accepted states in the native ODE algorithm; callbacks remain the
caller's responsibility. Rejected trials leave the last accepted state intact.
Projection counts, total concentration increment, and maximum increment are
reported. Conservation is separately tested on analytical cases.

For normalized error `e`, floored at `1e-12` for prediction, estimator order
`q`, and safety `s`, growth factors are:

| Controller | Factor after a successful previous step |
| --- | --- |
| `i` | `s * e^(-1/q)` |
| `pi` | `s * e^(-(1-beta)/q) * e_previous^(beta/q)` |
| `h211b` | `s^(2/(b*k)) * (e*e_previous)^(-1/(b*k)) * (h/h_previous)^(-1/b)` |

The filter follows equation (6) of Dreger et al. with target error `s` and
common growth bounds. `filter_k=0` selects `k=q`. Without valid history, and
after rejection, use the integral rule. Every controller keeps the same error
acceptance condition `e <= 1`. Growth is bounded by `factor_min`/`factor_max`.
The first accepted trial after rejection cannot grow. A fresh-matrix rejection
shrinks h; a stale matrix gets one retry at the same h with a fresh Jacobian.
No minimum-step forced acceptance is permitted. Time stagnation or exhausted
`max_steps` returns failure, a valid partial state, and `stats.final_time`.

`jacobian_steps > 1` lags a full analytical Jacobian for at most that many
accepted steps. This is **not** a formal W-method and can reduce order. With
`lu_band > 1`, hold h when proposed growth lies in `[1/lu_band, lu_band]` and
the Jacobian remains reusable, allowing exact reuse of that matrix's LU.
Changed h always requires refactorization. Every external call invalidates
matrices because rates, conditions or concentrations may have changed.
`warm_start=1` carries only a step prediction across calls. The final shortened
step does not erase a larger successful prediction. `reset` discards history.

### Added-method runtime options

Invalid/nonfinite values fail at setup. Resource limits too restrictive to
finish an interval cause explicit failure. No parameter depends on a species
name, mechanism, input seed, or assumed SIMD width.

| Key | Default | Domain / meaning |
| --- | ---: | --- |
| `rtol` | `1e-4` | Nonnegative shared relative tolerance |
| `atol` | `1e-18` | Nonnegative mol m-3; not zero together with rtol |
| `group_size` | `32` | Positive integer; short last group supported |
| `controller` | `i` | `i`, `pi`, `h211b` |
| `norm` | `cell` | `rms`, `cell`, `max` |
| `jacobian_steps` | `1` | Positive integer |
| `lu_band` | `1` | At least 1; 1 disables holding |
| `warm_start` | `1` | 0 or 1 |
| `h_start` | `0.01` | Positive seconds |
| `h_max` | `0` | Nonnegative seconds; 0 means the external interval |
| `max_steps` | `100000` | Positive integer; all trials count |
| `safety` | `0.9` | Between 0 and 1 |
| `factor_min` | `0.2` | Between 0 and 1 |
| `factor_max` | `6` | At least 1 |
| `pi_beta` | `0.4` | 0 through 1, PI only |
| `filter_b` | `4` | At least 1, H211b only |
| `filter_k` | `0` | Nonnegative; 0 selects estimator order |
| `positivity` | `0.1` | Nonnegative multiple of atol |
| `newton_iterations` | `8` | Positive integer; per SDIRK stage |
| `newton_tolerance` | `0.03` | Between 0 and 1; SDIRK correction norm |

These tolerance defaults belong to the scenario runner; library callers set
state tolerances explicitly. `--threads N`, `OMP_NUM_THREADS`, `OMP_SCHEDULE`
and `CH_GROUP_SIZE` retain Tier 1 meanings. By default all cores are used,
capped by group count. CLI group size overrides environment, then build
`GROUP_SIZE`. All numerical storage and factors are double precision.

Timing sums only photolysis interpolation, parameter/rate updates, and group
solves. Input, setup, diagnostic checks, snapshot copies, and CSV writing are
outside the timer. Both CSV files use 17 significant digits, preserve cell
order, and include the final state; snapshot steps are multiples of the stated
interval and always the last step. Photolysis is interpolated once at each
external step's start and held constant throughout that solve.

The workspace borrows one state; free it first and never solve the same state
concurrently. Reacquire state value pointers after solving. Added statistics
count reuse, clipping, and nonlinear work. `clipped_amount` sums concentration
increments across all cells/species/accepted steps; `maximum_clip` takes their
maximum. Legacy statistics still count work per group, not per cell.

### Adding a method

Add a Rosenbrock table in `table` and a runner name in `method_id`, or add a
family-specific trial routine writing `new_y`/`error` without modifying `y`.
Report RHS, Jacobian, LU, solve and nonlinear work. State workspaces hold six
stage vectors; larger formulas must allocate additional storage. Add an
independent prototype, order/analytical/failure tests, and all-scenario scans
through `tools/scan.py`. Only measured configurations enter the operating-point
list. Legacy formulas and arithmetic remain independently testable.

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
relative tolerance `1e-6` and species absolute tolerance `1e-3`. The Tier 1 TS1 benchmark uses those
defaults and the unmodified three-stage table. That benchmark never overrides
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
increase working-set size. The Tier 1 16/32/64/128/256 sweep favored 32; the final
Tier 3 comparison also selected 32, with a modest advantage over 16 described in
[SCANS.md](SCANS.md). This is a configurable default, not a numerical requirement.
The primary runtime uses all cores and dynamic group assignment to allow
different group costs. The retained benchmark schedules a group's complete
history at once. The scenario runner schedules groups at each external step,
retaining each group's private solver history.

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
