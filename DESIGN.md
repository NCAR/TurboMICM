# Architecture

The mechanism owns species, reaction, rate-constant, and constraint data. `chem_model_reorder` optionally applies the reference diagonal Markowitz ordering. Finalization establishes the stable reaction order, custom-parameter labels, algebraic mask, sparse Jacobian structure, fill pattern, and sparse factorization/substitution plans. Mutable state and solver workspaces belong to each cell group; the finalized mechanism is shared across threads.

Species and sparse matrix elements are stored in separate contiguous arrays over cells: `values[variable * stride + cell]`, and `jacobian[element * stride + cell]`. A state normally has `stride == cells`. Six stage vectors and three other work vectors share one allocation. Solver calls preserve caller-visible state pointers. OpenMP SIMD loops traverse cells; outer OpenMP loops distribute independent groups. Each group has its own error norm and adaptive history. Thus statistics scale with group count, rather than counting every cell separately.

The tuning parameters are:

| Parameter | Default | Set through | Effect |
| --- | --- | --- | --- |
| Cell group size | 32 | `make GROUP_SIZE=N`, or `--group N` | Cells sharing workspaces and an adaptive controller |
| Worker count | `omp_get_max_threads()` | `--threads N`, `OMP_NUM_THREADS` | Concurrent independent groups |
| Group scheduling | static | `OMP_SCHEDULE`, e.g. `dynamic,1` | Assignment of groups to workers |
| Generated kernel chunk | 32 reaction/derivative blocks | `make KERNEL_CHUNK=N` or generator `--kernel-chunk N` | Bounds generated function size and compiler work |
| SIMD width | compiler chosen | Compiler's ordinary OpenMP implementation | No fixed vector width or architecture assumption |

Positive group and chunk sizes, including partial final groups, use the same formulas. Scheduling never moves a cell between adaptive groups. Changes to grouping can change the roundoff and adaptation history, as in the reference OpenMP tutorial. No tuned value changes tolerances or coefficients.

`tools/configure.py` reads YAML or JSON, resolves `camp-files` relative to each source, validates species/reaction references, expands branches and repeated reactants, and normalizes rate units. Configuration rate coefficients default to molecule/cm3 units, as demonstrated by the tutorial's by-hand conversions; `concentration units: mol m-3` selects coefficients already in SI. Surface properties and custom parameter values are SI. The ordered TS1 sequences retain the exact source reaction arithmetic. See `DECISIONS.md` for the extensions.

The reader emits public-API model construction and optional forcing/Jacobian kernels. It calculates the same permutation and symbolic fill as the library so generated constant indices are verified against the finalized model. It performs no scientific parameter fitting or solver specialization. A generated loader `chem_load_NAME(int specialized)` selects either generated or generic evaluation for the same model. Generated rate evaluation, LU, controller, stages, and state storage remain ordinary library functions. The benchmark selects generated forcing/Jacobian kernels.

For another mechanism, run:

```sh
.conda/bin/python tools/configure.py mechanism.yaml generated/example.c --name example
```

Add that generated source to a Makefile target compiled with `$(CC) $(CFLAGS)`. Call `chem_load_example(1)`, create a state, assign conditions/concentrations/custom parameters, call `chem_calculate_rates`, and solve. Python is never involved at runtime. The retained `tools/translate_ts1.py` regenerates `configs/ts1.json` from the supplied TS1 declarations.

The Rosenbrock solver receives a parameter table, a step-size controller function, and factor/solve function pointers. All five reference parameter sets are in `src/parameters.c`; the retained translation script records their provenance. To add a method, fill a `ChemRosenbrock` table with its stages, lower-triangular `a/c`, weights `m/e`, `new_function`, `alpha/gamma`, and estimator order. Assign it to `ChemSolver.parameters`; no integration-loop change is needed.

To replace the controller, assign `solver.controller` and `solver.controller_context`. The callback returns an error-dependent step multiplier. Acceptance/rejection and post-rejection bounds remain explicit in `src/rosenbrock.c`; modifying those policies is localized there. The default is `chem_default_controller`, using the original safety factor, bounds, and estimator order.

Linear algebra defaults to unpivoted in-place Mozart elimination. `solver.factor = chem_lu_doolittle` switches to Doolittle; both use `chem_linear_solve` and the finalized sparse pattern. A new factor/solve pair can replace both callbacks using the same element-major layout and symbolic metadata. Separate L/U APIs are also provided. There is no alternate linear algebra hidden inside the benchmark.

To add a rate type, extend `ChemRateType`, add its parameters to `ChemRate`, initialize defaults in `chem_rate`, and evaluate it in `chem_rate_value`. Update the enum range validation in `chem_model_reaction`, stable configuration type lists, and field/unit normalization in the reader. Generated forcing and Jacobian code requires no change because it consumes computed rate constants. C callbacks (`CHEM_LAMBDA`) support exploratory rates without recompiling the library. Taylor polynomials support the reference limit of 16 coefficients.

Backward Euler has its own `ChemEuler` controls and the reference Newton/reduction sequence. `chem_advance` accepts a generic step callback so experiments can replace an entire integrator. It resumes positive-progress partial solves and aggregates statistics. `chem_advance_many` distributes distinct states over OpenMP workers.

Linear and equilibrium constraints, plus caller-defined residual/Jacobian callbacks, share the same mass-matrix mask. DAE initialization solves for algebraic corrections, uses state-weighted convergence and backtracking, and rolls back failed initialization. Algebraic step changes enter the error estimate. After a solve, only differential variables are clamped nonnegative, matching the reference facade. External process Jacobian entries overlapping constraint dependencies are replaced by constraint entries; unrelated process entries in algebraic rows remain declared and active.

All model storage is owned by the model, all state storage by the state, and states must be destroyed before their model. Reaction and constraint insertion copies caller arrays and labels. Species properties support owned strings, doubles, integers, and booleans; phases supply qualified names and deep copies. Setters return error codes for invalid names and sizes. `chem_state_copy` creates independent writable state and scratch storage.

External components declare their species, parameters, and constraints through the builder APIs, then register a `ChemExternal` callback record with `chem_model_external`. Setup runs once after finalization and can cache sparse and variable indices. Parameter updates run before rate evaluation; forcing and negative-Jacobian callbacks run on each solver evaluation. Multiple components compose in insertion order. Callback contexts are borrowed: callers own their lifetime and provide reentrant callbacks when sharing a model across threads. An external Jacobian can provide a declared sparsity pattern; omission conservatively allocates a dense one. `tests/external_helpers.py` and the aerosol/constraint fixtures demonstrate the contract.

Array utilities independently expose dense/grouped layouts, CSR/CSC sparsity, explicit read-only views, shape-checked callback invocation, and reductions that can exclude padding. Solver kernels always traverse actual cells. Finite-difference verification and sparsity checking live in `src/verify.c` and are independent of the integrators.

The runtime sources include only standard C headers and `omp.h`. SIMD comes from contiguous data, `restrict`, and `omp simd`. Build flags add language, warnings, includes, position-independent code, and OpenMP linkage to the caller's unchanged `OPTFLAGS`. No intrinsics, target attributes, OS APIs, wider floating-point types, fast-math, or tuning flags are used.

# Tier 3 integration and scenario interface

`include/chem_adaptive.h` exposes the opt-in integration API. A `ChemWorkspace` belongs to exactly one state and must be destroyed before that state. It owns an unfactored Jacobian, seven solution-history vectors, and scratch vectors. It preserves the proposed next step across calls. `chem_workspace_reset` discards multistep, matrix, and controller history after a forcing discontinuity, retaining only that step proposal. Call it whenever concentrations, conditions, rates, tolerances, or the method are changed externally; uninterrupted calls with the same ODE may keep the history. These solvers support ordinary chemical ODEs; constraint/DAE models continue to use the inherited API. Both generic reaction evaluation and generated kernels work without integrator changes.

| Scenario method | Formula and implementation |
| --- | --- |
| `micm2` … `micm6` | Exact inherited parameter tables and `chem_solve`, including its original startup, norm, rejection statistics, and end-of-call clamp. `micm3` is the default. |
| `ros2` | Two-stage, second-order L-stable formula; the coefficients also define a second-order W-method with an approximate Jacobian. |
| `ros3` | Three stages, two forcing evaluations, third order with a second-order embedding. |
| `ros4` | Four stages, three forcing evaluations, fourth order with an embedded estimator. |
| `rodas3` | Four stages, three forcing evaluations, stiffly accurate third-order formula. Selected for the operating points. |
| `rodas4` | Six stages, six forcing evaluations, stiffly accurate fourth-order formula. |
| `bdf` | Independently implemented arbitrary-step BDF, variable order up to the requested maximum (1–5), modified Newton iteration and reusable sparse matrices. |

The five new Rosenbrock modes share `src/adaptive.c`. They use the inherited coefficients, not the inherited integration loop. A trial forms `(I/(h gamma) - J) k_i = f(Y_i) + sum(c_ij k_j)/h`, then forms the solution and embedded error from the stage weights. Chemistry is autonomous within an outer interval. There is no derivative of the time-dependent light factor within that interval.

New error control uses `atol_i + rtol * max(abs(y_i),abs(ynew_i))`. `norm=0` is RMS across all species and cells; `norm=1` takes each cell's species RMS and then the maximum across cells; `norm=2` takes the largest scaled component. The latter two avoid dilution by adding easy cells to a group. A trial must have error at most one and finite concentrations. It is rejected if any concentration is below `-positivity * (atol_i + rtol*abs(y_i))`. Smaller negative roundoff/truncation values are projected to zero only after acceptance. This is an explicit numerical policy, not a proof of atom conservation. `clipped_values` and `largest_clip` quantify the projection. Rejected trials retain the old state and shrink the step; nonfinite stages produce a failed error test and a smaller retry. There is no acceptance at a minimum step and no silent fallback to different chemistry.

The proportional controller proposes `safety * error^(-1/p)`, where `p` is the embedded error order from the coefficient table. The PI alternative proposes `safety * error^(-0.7/p) * previous_error^(0.4/p)`. The H211b alternative proposes `safety * (error*previous_error)^(-1/(b*k)) * (h/previous_h)^(-1/b)`; it uses proportional control until an accepted history exists and after rejection. All growth factors are bounded; the next trial after a rejection cannot grow. History is reset when the light factor changes. See `REFERENCES.md` for the mathematical sources. These controller policies are independent of the chemistry.

Jacobian age counts accepted steps. A rejection refreshes a lagged Jacobian before further retrying. Factored matrices may be reused only with an identical step size and the same Jacobian. `step_hold` supplies a symmetric deadband for the step multiplier, which can make that exact equality useful. Lagging `ros3`, `ros4`, or either RODAS formula is an experiment: their full formal order with an arbitrary approximate Jacobian is not claimed. Declared operating points always use a fresh Jacobian. No species-dependent selection is involved.

`src/bdf.c` differentiates a Lagrange interpolant through the new endpoint and the previous `q` states to obtain the BDF coefficients. With dimensionless nodes `x_j=(t_old_j-h)/h`, `d_0=sum(-1/x_j)` and the other derivative weights follow directly from the Lagrange basis. Newton solves the full residual with shift `alpha=d_0/h`. The solution error estimate is the difference from extrapolation through `q+1` old points, divided by `d_0*(1-t_old_q/h)`. Startup uses backward Euler against an explicit tangent predictor. The extra historical point is required before an order can be raised. This implementation uses polynomial history rather than a Nordsieck representation and is not CVODE or NDF.

BDF attempts at most eight chord iterations, stopping when the weighted correction is below 0.03. A stalled iteration first refreshes an old matrix and then reduces the step if necessary. Trial steps are never accepted without nonlinear convergence. A factorization can be reused when its shift differs by at most `reuse_ratio`, its age is below `jacobian_age`, and Newton convergence permits it. Failed Newton iterations invalidate it. Order selection compares neighboring-order divided-difference estimates after enough accepted steps, with an upper bound set by `order`. These fixed numerical policies are visible in one source file; BDF is retained as a tested research alternative, not a recommended operating point.

## Runtime parameters

Every option is passed as a separate `--option key=value`. Invalid, nonfinite, or out-of-range values are rejected. `--threads N` selects workers; `--group N` is an alias for `--option group=N`. Positive sizes need not divide the cell count. The solver uses actual lanes only.

| Option | New-method default | Meaning / range |
| --- | ---: | --- |
| `rtol` | `1e-4` | Positive relative error tolerance |
| `atol` | `1e-18` mol m-3 | Positive uniform absolute error tolerance; the library API also supports species-specific values |
| `group` | build `GROUP_SIZE`, normally 32 | Positive cells per adaptive group |
| `schedule` | `OMP_SCHEDULE`, otherwise static | Driver option: `static`, `dynamic`, or `guided`; an explicit option overrides the environment |
| `schedule_chunk` | runtime default | Driver option: nonnegative integer chunk size; dynamic/guided default to 1 when explicitly selected |
| `h_start` | 1 second | Positive first step proposal; subsequent calls retain their next-step proposal |
| `h_max` | 0 | Nonnegative internal-step cap; zero means the outer interval |
| `max_steps` | 200000 | Positive trial-step resource limit per outer interval |
| `safety` | 0.85 | Positive finite controller multiplier |
| `growth` | 3 | Maximum multiplier, at least 1 |
| `norm` | 1 | 0 group RMS, 1 maximum cell RMS, 2 maximum component |
| `controller` | 0 | Rosenbrock: 0 proportional, 1 PI, 2 H211b |
| `controller_b` | 1 | H211b smoothing parameter, at least 1 |
| `controller_k` | 2 | Positive H211b error exponent parameter |
| `step_hold` | 1 | Rosenbrock deadband: hold `h` when the multiplier is in `[1/step_hold,step_hold]`; 1 disables it |
| `jacobian_age` | 1 accepted step | Positive maximum Jacobian/matrix age |
| `reuse` | 1 | 0 disables matrix reuse, 1 permits it subject to the method's checks |
| `reuse_ratio` | 0.2 | BDF: allowed relative change of matrix shift, in `[0,1)` |
| `order` | 5 | BDF maximum order, integer 1–5 |
| `positivity` | 1 | Nonnegative allowance in units of local concentration scale before rejection |

The shrink floor is 0.1 and rejection always forces a reduction. `h_start`, `h_max`, tolerances and iteration/step limits are numerical controls: deliberately restrictive choices may exhaust the resource limit, which returns a failure rather than an incomplete solution. The inherited modes allow only `rtol`, `atol`, `group`, and the scheduling options; omitted tolerances retain the exact original values (`1e-6`, `1e-3` mol m-3). The default path is not an operating point for the long scenarios.

## Driver and timing

`src/scenario.c` reads the standard `key = value` files, resolves paths relative to the scenario file, reads every cell, and evaluates the piecewise linear photolysis curve at the start of each step. Conditions and nonphotolysis parameters remain fixed. Initial rate calculation is setup. During each timed step it sets every custom label beginning with `j` to its original value times the factor and updates the corresponding user-defined reaction constants from their original, already computed values. This is exact for TS1's linear photolysis rates, including third-body factors; `tests/test_scenario.py` independently recalculates all rates and checks the same solution. All thermal rates, species, and reactions remain present.

The clock includes interpolation, updating light-dependent parameters/rates, parallel group scheduling, integration, and aggregation of operation counts. It excludes parsing, allocation, model construction, initial rate calculation, finite-output validation, and all CSV serialization. There is no warm-up. Each interval is timed separately so snapshots can be written outside the clock, without buffering an entire run. `ts1_final.csv` and `ts1_snapshots.csv` use 17-digit values and the supplied ordering rules, including an unconditional final snapshot when snapshots are requested. A scenario without snapshots produces a header-only snapshot file.

`internal_steps`, forcing, Jacobian, factorization, solve, and Newton counts count group-wide operations, not individual lanes. New rejection counts include startup rejections; inherited counts retain MICM's original convention. `matrix_reuses` counts trials initially using an old matrix; a BDF trial may later refresh it. `clipped_values` counts scalar values; `largest_clip` is in mol m-3. Every outer endpoint is checked for finite, nonnegative output. Any failure gives a nonzero process exit; partial output must not be scored as a completed run.

## Adding a method

For another Rosenbrock table, add its coefficients in `src/parameters.c`, a method identifier in the opt-in API's validation, and a name in `src/scenario.c`. The common loop handles its stages and estimator; no reaction-specific code is needed. For another integrator family, add an implementation beside `src/bdf.c` and dispatch it in `chem_adaptive_solve`, using a workspace owned by the state. New counters belong in `ChemAdaptiveStats`. Keep the inherited `chem_solve` path unchanged.

Add an independent small Python prototype in `prototypes/`, bind any new public fields in `tests/adaptive_api.py`, and test state and operation-count agreement in `tests/test_adaptive.py`. Then add controlled configurations to a JSON file in `scans/` and run `tools/scan.py`, which invokes the unmodified supplied `score.py`, preserves failed attempts, and records native time and operation counts. Timings must be collected without concurrent compiler or benchmark jobs. No Python is loaded by either native executable.
