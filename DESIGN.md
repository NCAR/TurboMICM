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
