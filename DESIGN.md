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

The complete generator adds these build parameters. Command-line options use lower-case dashed names, except `SPECIALIZED_STRIDE` maps to `--cell-stride` and `LINEAR_MODE` to `--linear`. The direct CLI independently defaults both cell stride and solver tile to 4; Make derives the stride default from the selected tile/group size as shown below.

| Make parameter | Default | Effect |
| --- | --- | --- |
| `SOLVER_TILE` | 4 | Cells packed into a shared temporary chemistry/linear workspace inside each adaptive group; 0 uses the unblocked three-stage template |
| `SPECIALIZED_STRIDE` | tile size, or group size when tiles are disabled | Emit an additional backend with that literal cell stride/count; 0 disables it. Other strides/counts use the general generated backend |
| `LINEAR_CHUNK` | 1 pivot/row | Number of consecutive factor pivots or solve rows per generated helper |
| `LINEAR_LOOP` | `rows` | `rows` fuses each helper's scalar operations inside a SIMD cell loop; `operations` emits separate loops, `pointers` also declares distinct element pointers, and `dot` accumulates each LU element/solve row in a register |
| `LINEAR_MODE` | `unrolled` | `generic` uses the generic factor/solve callbacks for ablation; chemistry is still generated |
| `FUSE_KERNELS` | 1 | Fuse each chunk's reaction/derivative blocks inside one cell loop; 0 emits a cell loop per block |
| `INITIALIZE_KERNELS` | 1 | Assign each chemistry output on its first contribution and explicitly zero untouched elements; 0 clears the full output before adding terms |
| `VECTOR_ERROR` | 1 | SIMD computation of squared scaled errors, followed by the unchanged serial group sum; 0 scales during the serial sum |
| `FUSE_STAGES` | 1 | Use the generated three-stage pipeline; 0 uses the generic Rosenbrock loop with generated chemistry/LU |
| `SPECIALIZE_RATES` | 1 | Emit fixed rate formulas and custom-parameter indices; 0 calls the generic rate evaluator |

Group size, tile size, generation chunks, worker count, and scheduling are independent. A tile never changes the adaptive group or its norm. Partial final tiles use actual cells only; no padding enters chemistry or error control. Positive sizes are supported subject to ordinary allocation/compiler resource limits; 0 is accepted only for the disabling options above. No value modifies numerical coefficients, tolerances, or precision.

For TS1, a four-cell packed workspace has `(2359 + 547 + 7*209)*4 = 17,476` doubles (139,808 bytes): matrix, rates, concentration, forcing, error, new state, and three stages. A 32-cell state's three unused stage vectors provide 20,064 doubles, so the main benchmark needs no tile allocation. Other groups/tiles use a checked heap workspace when needed. Error terms are copied back to species-major group storage before the serial reduction. Initial forcing is retained across rejected attempts exactly as in the generic method; each attempt recomputes each cell's Jacobian and factorization. Logical group statistics are unchanged.

Positive group and chunk sizes, including partial final groups, use the same formulas. Scheduling never moves a cell between adaptive groups. Changes to grouping can change the roundoff and adaptation history, as in the reference OpenMP tutorial. No tuned value changes tolerances or coefficients.

`tools/configure.py` reads YAML or JSON, resolves `camp-files` relative to each source, validates species/reaction references, expands branches and repeated reactants, and normalizes rate units. Configuration rate coefficients default to molecule/cm3 units, as demonstrated by the tutorial's by-hand conversions; `concentration units: mol m-3` selects coefficients already in SI. Surface properties and custom parameter values are SI. The ordered TS1 sequences retain the exact source reaction arithmetic. See `DECISIONS.md` for the extensions.

The Tier 1 reader emits public-API model construction and optional forcing/Jacobian kernels. `mechanism_layout` calculates the same permutation and symbolic fill as the library so generated constant indices are verified against the finalized model. A generated loader `chem_load_NAME(int specialized)` selects generated or generic evaluation for the same model. Tier 1 rate evaluation, LU, controller, stages, and state storage remain ordinary library functions. Its benchmark still selects generated forcing/Jacobian kernels.

`tools/specialize.py` uses the same normalized configuration and symbolic layout to emit a complete three-stage solver. Rate dispatch, parameter lookups, chemistry indices, derivative jobs, sparse elimination, and substitutions become C literals and straight-line statements inside cell loops. The loader still constructs generic metadata for the public state API and independent verification; this setup is outside timing. The specialized integration template, `tools/rosenbrock3.c.in`, retains the generic adaptation, rejection, clamp, error floor, step bounds, maximum-step checks, and statistics. It unrolls the three stages and fuses vector passes while keeping each scalar addition in order. Coefficients and controller constants come from `chem_rosenbrock_parameters(3)` and `chem_default_controller`, shared with the generic library.

Conditions, arbitrary per-cell air density, custom parameters, concentrations, tolerances, and interval length remain runtime data. Built-in rate coefficients are fixed by configuration; regenerating changes them. Derived branched and surface constants are computed once by the loader with the same C/libm expressions as the generic library. No generated kernel reads another mechanism or uses TS1-specific reaction assumptions. Configurations retain the Tier 1 serializable rate families; C callbacks remain generic-library features because a function pointer has no configuration representation.

The generated solver exports `chem_load_NAME(1)`, `chem_rates_NAME(state)`, `chem_solve_NAME(dt,state,stats)`, `chem_factor_NAME`, and `chem_linear_NAME`. States must come from that mechanism and retain their original dimensions and contiguous cell stride. The rate function must be called after changing conditions/custom parameters; concentration changes do not require it. A model and its code are immutable while states run, and separate states can share them across OpenMP workers. The generic API remains independently available.

For another mechanism, run:

```sh
.conda/bin/python tools/configure.py mechanism.yaml generated/example.c --name example
```

Add that generated source to a Makefile target compiled with `$(CC) $(CFLAGS)`. Call `chem_load_example(1)`, create a state, assign conditions/concentrations/custom parameters, call `chem_calculate_rates`, and solve. Python is never involved at runtime. The retained `tools/translate_ts1.py` regenerates `configs/ts1.json` from the supplied TS1 declarations.

Use `tools/specialize.py mechanism.yaml generated/example_specialized.c --name example_specialized` for the complete specialized variant. `make bench` regenerates TS1 whenever its configuration, reader, generator, template, or specialization options change. `make test-specialized` also regenerates every configuration example and the modified TS1 fixture. Generated files live under `generated/`, are excluded from git, and are never edited. Build flag stamps force recompilation when switching compiler flags or group defaults; generation flag stamps force regeneration when specialization choices change. Compiler-specific object directories allow compiler switching. GNU Make 3.81 uses explicit static pattern targets for specialized outputs and scoped intermediate retention.

The Rosenbrock solver receives a parameter table, a step-size controller function, and factor/solve function pointers. All five reference parameter sets are in `src/parameters.c`; the retained translation script records their provenance. To add a method, fill a `ChemRosenbrock` table with its stages, lower-triangular `a/c`, weights `m/e`, `new_function`, `alpha/gamma`, and estimator order. Assign it to `ChemSolver.parameters`; no integration-loop change is needed.

To replace the controller, assign `solver.controller` and `solver.controller_context`. The callback returns an error-dependent step multiplier. Acceptance/rejection and post-rejection bounds remain explicit in `src/rosenbrock.c`; modifying those policies is localized there. The default is `chem_default_controller`, using the original safety factor, bounds, and estimator order.

Linear algebra defaults to unpivoted in-place Mozart elimination. `solver.factor = chem_lu_doolittle` switches to Doolittle; both use `chem_linear_solve` and the finalized sparse pattern. A new factor/solve pair can replace both callbacks using the same element-major layout and symbolic metadata. Separate L/U APIs are also provided. There is no alternate linear algebra hidden inside the benchmark.

To add a rate type, extend `ChemRateType`, add its parameters to `ChemRate`, initialize defaults in `chem_rate`, and evaluate it in `chem_rate_value`. Update the enum range validation in `chem_model_reaction`, stable configuration type lists, and field/unit normalization in the reader. Tier 1's generated forcing and Jacobian code consumes computed rate constants and requires no change. The complete generator also needs that rate's literal expression in `rate_body` in `tools/specialize.py`; compare it with generic evaluation in the validation configuration. C callbacks (`CHEM_LAMBDA`) support exploratory rates in the generic library. Taylor polynomials support the reference limit of 16 coefficients.

Backward Euler has its own `ChemEuler` controls and the reference Newton/reduction sequence. `chem_advance` accepts a generic step callback so experiments can replace an entire integrator. It resumes positive-progress partial solves and aggregates statistics. `chem_advance_many` distributes distinct states over OpenMP workers.

Linear and equilibrium constraints, plus caller-defined residual/Jacobian callbacks, share the same mass-matrix mask. DAE initialization solves for algebraic corrections, uses state-weighted convergence and backtracking, and rolls back failed initialization. Algebraic step changes enter the error estimate. After a solve, only differential variables are clamped nonnegative, matching the reference facade. External process Jacobian entries overlapping constraint dependencies are replaced by constraint entries; unrelated process entries in algebraic rows remain declared and active.

All model storage is owned by the model, all state storage by the state, and states must be destroyed before their model. Reaction and constraint insertion copies caller arrays and labels. Species properties support owned strings, doubles, integers, and booleans; phases supply qualified names and deep copies. Setters return error codes for invalid names and sizes. `chem_state_copy` creates independent writable state and scratch storage.

External components declare their species, parameters, and constraints through the builder APIs, then register a `ChemExternal` callback record with `chem_model_external`. Setup runs once after finalization and can cache sparse and variable indices. Parameter updates run before rate evaluation; forcing and negative-Jacobian callbacks run on each solver evaluation. Multiple components compose in insertion order. Callback contexts are borrowed: callers own their lifetime and provide reentrant callbacks when sharing a model across threads. An external Jacobian can provide a declared sparsity pattern; omission conservatively allocates a dense one. `tests/external_helpers.py` and the aerosol/constraint fixtures demonstrate the contract.

Array utilities independently expose dense/grouped layouts, CSR/CSC sparsity, explicit read-only views, shape-checked callback invocation, and reductions that can exclude padding. Solver kernels always traverse actual cells. Finite-difference verification and sparsity checking live in `src/verify.c` and are independent of the integrators.

The runtime sources include only standard C headers and `omp.h`. SIMD comes from contiguous data, `restrict`, and `omp simd`. Build flags add language, warnings, includes, position-independent code, and OpenMP linkage to the caller's unchanged `OPTFLAGS`. No intrinsics, target attributes, OS APIs, wider floating-point types, fast-math, or tuning flags are used.
