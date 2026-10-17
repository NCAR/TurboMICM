# Design

A C11 library (`libtmc.a`, header `include/tmc.h`) that reproduces MICM's CPU
chemistry solver: species, phases, and reactions with all of MICM's
rate-constant types; a solver builder that compiles a mechanism into plans;
states that hold groups of grid cells; Rosenbrock and backward Euler
integrators; sparse LU factorizations; constraints (differential-algebraic
systems) and external models. Python tooling reads MICM's configuration
format and generates, for any mechanism, C tables and straight-line kernels
that the library can run in place of its generic loops. The benchmark
`bin/ts1_bench_cpu` runs TS1 through the same library code the tests
exercise, with TS1's generated kernels attached.

## Architecture

```
tools/translate_ts1.py   ts1_mechanism.hpp (C++) -> configs/ts1.json
tools/mechanism.py       reads YAML/JSON, compiles species order, sparsity, LU fill, fingerprint
tools/generate.py        -> generated/<name>.c: builder calls, initial state, kernels (forcing, Jacobian, LU, solves)

include/tmc.h            the public API
src/species.c            species, properties, phases, systems
src/reaction.c           reactions and their builder
src/rates.c              rate-constant constructors and MICM's formulas
src/parameters.c         the five Rosenbrock parameter sets, backward Euler defaults, the default controller
src/dense.c, sparse.c    vector-ordered dense and block-sparse arrays
src/lu.c                 Markowitz ordering, symbolic fill, Mozart and Doolittle plans, decompositions, solves
src/builder.c            SolverBuilder::Build: names, ordering, jobs, pattern, constraints, external models, fingerprint
src/state.c              states: concentrations, conditions, parameters, tolerances, work arrays
src/process.c            rate constants of a state; generic forcing and Jacobian; complete evaluation with externals and constraints
src/constraints.c        equilibrium and linear constraints
src/external.c           external model dispatch and name maps
src/rosenbrock.c         the Rosenbrock loop, the DAE path, the error norm
src/backward_euler.c     the backward Euler loop
src/util.c               errors, allocation, hashing, version
bench/ts1_bench.c        the benchmark driver; bench/kernel_bench.c times the kernels one by one
tests/                   the ported test suite (tests/README.md)
```

A **solver** (`tmc_solver`) is immutable after `tmc_builder_build` and is
shared by every thread. It holds the species order, the sorted reactions,
the process plans (reactant and product indices, yields, third bodies,
custom-parameter slots), the Jacobian jobs in MICM's order, the Jacobian and
LU patterns with their replay lists, the constraints and their plans, the
external models, the mass-matrix diagonal, and optionally the generated
kernels. A **state** (`tmc_state`) is one group of cells with its own step
size: concentrations, rate constants, custom parameters, conditions,
tolerances, the LU matrix, pivot reciprocals, and the integrator's work
arrays. A benchmark or host model makes one state per group and hands groups
to threads; the library itself opens no parallel region.

### Layout

Every per-cell array is **vector ordered** with width L: cells form blocks of
L lanes and within a block the lane index is fastest, so a block of species
values is `data[(block * n_vars + species) * L + lane]` and a block of matrix
elements `data[(block * nnz + element) * L + lane]`. L = 1 is MICM's standard
(row-major) ordering; L = W = `TMC_WIDTH` (2 by default) is the width the
generated kernels are compiled for. The last block of a state is padded; the
padded lanes hold zeros, compute harmlessly, and are excluded from the error
norm and the convergence tests. A dense `tmc_dense` and a sparse `tmc_sparse`
carry their own L, so the utility arrays serve every layout the tests ask
for while the solver's states use one.

### From a configuration to C

1. `tools/mechanism.py` reads the configuration (resolving `camp-files`),
   converts units as MICM's tutorial programs do by hand, and expands
   branched reactions, quantities, and yields. TS1 comes from
   `tools/translate_ts1.py`, which text-matches the C++ declarations.
2. It compiles the mechanism the way `src/builder.c` does: reactions
   stably sorted by rate-constant type; species in declaration order and
   then MICM's diagonal Markowitz ordering of the Jacobian pattern; the
   Jacobian jobs; the pattern with the diagonal; the symbolic LU fill; and a
   fingerprint (FNV-1a over the species names, reactions, yields, and LU
   pattern) that the C builder recomputes.
3. `tools/generate.py` writes `<name>.c` and `<name>.h`:
   - `<name>_system()` and `<name>_reactions()` build the mechanism through
     the public builder API, so the generated and the hand-built mechanism
     are the same object to the library;
   - `<name>_initial_state()` sets the configuration's initial
     concentrations, parameters, and conditions;
   - `<name>_kernels()` returns a `tmc_kernels` table of four functions with
     the width and the fingerprint. The builder attaches them only when both
     match the solver it just built, and never to a solver with constraints
     or external models, which run the generic kernels.
4. The Makefile runs the generator at build time with `$(PYTHON)` from the
   conda environment and compiles the output with `$(CC)` and `$(OPTFLAGS)`.
   Python never runs when the library or the benchmark runs.

The four kernels act on one state's blocks in lock step (`nb` blocks):

- **forcing**: pass 1 computes every reaction rate k prod(reactants) into a
  per-thread buffer; pass 2 assigns every species its sum of -rate over the
  reactions consuming it and yield rate over those producing it, in MICM's
  accumulation order, four species' sums interleaved.
- **Jacobian**: pass 1 computes every (independent variable, reaction)
  partial k prod(other reactants); pass 2 assigns every LU element its sum
  (+partial for reactant rows, -yield partial for product rows), in MICM's
  job order, and adds alpha to the diagonal last, so the matrix is
  alpha I - J with zeros in the fill positions.
- **LU**: MICM's Mozart elimination in panels of four pivots: within a panel
  the pivots run one by one (reciprocal, column scaling, updates of the
  panel's rows and columns); the trailing elements then take the panel's
  rank-4 update, each loaded and stored once. Every element sees its updates
  in pivot order, so the factors equal the pivot-by-pivot algorithm's bit for
  bit. The pivot reciprocals are stored for the solve.
- **solve**: forward substitution with the unit lower factor, then backward
  substitution. In fast mode (the default) a row with twelve or more terms
  accumulates four partial sums and the back substitution multiplies by the
  stored reciprocal; exact mode (`make KERNEL_MODE=exact`) keeps MICM's
  single sum and division.

The generic kernels in `src/process.c` and `src/lu.c` do the same work from
index lists, for any width, and are what the tests mostly run; the tests also
compare generated and generic results.

### The solvers

`src/rosenbrock.c` is MICM's `AbstractRosenbrockSolver::Solve` with the same
control flow: the forcing and Jacobian at the start of a step, the shifted
matrix alpha M - J factored in place (recomputed after a rejection), the
stages with MICM's reuse of a forcing between stages, the combination into
the new solution and the error estimate, the normalized error accumulated
variable by variable and cell by cell, the bounded step-size factor, and the
acceptance rules. Only the vector arithmetic is restructured: the stage
combinations are single passes of fused multiply-adds in MICM's order. With
constraints, the Newton initialization of the algebraic variables (with
MICM's backtracking line search), the mass-matrix shift and coupling, the
algebraic error estimate, and the differential-only clamp follow MICM's code.
`src/backward_euler.c` is MICM's `AbstractBackwardEuler::Solve`.

### Threading

The library has no parallel region; the benchmark runs `omp parallel for
schedule(runtime)` over groups, dynamic one group at a time by default, for
the untimed warm-up and for the timed steps, and every group's states and
rate constants are made before the timer starts. `_Thread_local` buffers in
the generated kernels hold the reaction rates and Jacobian partials.

## Tuning parameters

| Parameter | Default | Set by | Effect |
|---|---|---|---|
| vector width W | 2 | `make WIDTH=` (`TMC_WIDTH`); `tmc_builder_set_width` per solver | lanes per block; the generated kernels are compiled for the build's W |
| group size | 8 | `make GROUP_SIZE=`; `--group-size` | cells per state, sharing one step size and error norm; blocks run in lock step |
| threads | all cores | `--threads`, `OMP_NUM_THREADS` | OpenMP threads over groups |
| schedule | dynamic, 1 | `--schedule static|dynamic|guided`, `--chunk` | how groups are handed to threads |
| LU panel | 4 | `make GEN_FLAGS="--panel N"` | pivots per generated LU panel |
| kernel mode | fast | `make KERNEL_MODE=exact` | fast: split sums and pivot reciprocals in the solve; exact: MICM's order bit for bit |
| kernel chunk | 16 | `make KERNEL_CHUNK=` | pivots or rows per generated function (compile time only) |
| kernels | generated | `--generic` | the benchmark runs the generic index-list kernels instead |

Every value gives correct results; the gate holds for each (RESULTS.md).
Changing W, the group size, or the kernel mode changes the roundoff and the
step-size history as MICM's own layouts and groupings do, within 1e-13 of
MICM on TS1.

## Making the changes the reviewers will try

- **A Rosenbrock parameter set.** Add a constructor beside the five in
  `src/parameters.c` returning a filled `tmc_rosenbrock_parameters` (stages,
  `a`, `c`, `m`, `e`, `alpha`, `gamma`, `new_function_evaluation`, the
  estimator order) and declare it in `tmc.h`. The loop in `src/rosenbrock.c`
  reads only the struct; up to `TMC_MAX_STAGES` (6) stages are supported.
- **The step-size controller.** Write a `tmc_step_controller` (error and
  parameters in, the factor Hnew/H out) and pass it with
  `tmc_builder_set_step_controller`; `tmc_default_step_controller` is MICM's
  rule. The acceptance and rejection logic stays in `tmc_rosenbrock_solve`
  and is the place to change if the policy itself changes.
- **The linear solver.** The algorithm is `tmc_builder_set_lu_algorithm`
  (`TMC_LU_MOZART` or `TMC_LU_DOOLITTLE`), both built as replay lists in
  `src/lu.c` from the symbolic pattern; a new elimination order is another
  plan builder and executor there, behind the same `tmc_lu_*` entry points.
  The generated LU and solve kernels are emitted by `emit_kernels` in
  `tools/generate.py`; a different algorithm for the kernels changes that
  emitter and nothing in the library, which only calls the function
  pointers in `tmc_kernels`.
- **A rate-constant type.** Add its enumerator to `tmc_rate_type` (keep
  `TMC_RATE_LAMBDA` last, the sort order follows the enum), its parameter
  struct to the union in `tmc_rate_constant`, a constructor and the formula
  in `src/rates.c` (`tmc_rate_constant_calculate`), and, if it reads custom
  parameters, its labels in `tmc_rate_constant_custom_parameter_labels`. For
  configurations, add its keys and unit conversion to `_reaction` in
  `tools/mechanism.py` and its emission to `emit_mechanism` in
  `tools/generate.py`. The kernels are unaffected: they read computed rate
  constants.
