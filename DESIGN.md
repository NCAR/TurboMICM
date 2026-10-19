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
src/*_ext.c, support.c   internals MICM's tests reach (include/tmc_*_ext.h, tmc_support.h): error norm, alpha shift, convergence test, constraint terms, array functions, name maps
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
   then an ordering of the Jacobian pattern, MICM's diagonal Markowitz or
   greedy minimum fill (the TS1 kernels' default; attached kernels impose
   their ordering on the solver they are built into); the
   Jacobian jobs; the pattern with the diagonal; the symbolic LU fill; and a
   fingerprint (FNV-1a over the species names, reactions, yields, and LU
   pattern) that the C builder recomputes.
3. `tools/generate.py` writes `<name>.c` and `<name>.h`:
   - `<name>_system()` and `<name>_reactions()` build the mechanism through
     the public builder API, so the generated and the hand-built mechanism
     are the same object to the library;
   - `<name>_initial_state()` sets the configuration's initial
     concentrations, parameters, and conditions;
   - `<name>_kernels()` returns a `tmc_kernels` table of four kernels and two
     fused variants (the forcing with the Jacobian at the start of a step; the
     forcing at a stage's combination y + a k), with the width and the
     fingerprint. The builder attaches them only when both
     match the solver it just built, and never to a solver with constraints
     or external models, which run the generic kernels.
4. The Makefile runs the generator at build time with `$(PYTHON)` from the
   conda environment and compiles the output with `$(CC)` and `$(OPTFLAGS)`.
   Python never runs when the library or the benchmark runs.

The kernels act on one state's blocks in lock step (`nb` blocks):

- **forcing**: pass 1 computes every reaction rate k prod(reactants) into a
  per-thread buffer; pass 2 assigns every species its sum of -rate over the
  reactions consuming it and yield rate over those producing it, in MICM's
  accumulation order, four species' sums interleaved.
- **Jacobian**: pass 1 computes every (independent variable, reaction)
  partial k prod(other reactants); pass 2 assigns every LU element its sum
  (+partial for reactant rows, -yield partial for product rows), in MICM's
  job order, and adds alpha to the diagonal last, so the matrix is
  alpha I - J with zeros in the fill positions.
- **LU**: MICM's Mozart elimination in panels of 64 pivots, every element
  finished target by target. At its pivot: the diagonal, then the column
  below it (its updates from the panel's earlier pivots in pivot order, then
  the scaling by the reciprocal), then the row to its right; after the
  panel, the trailing elements take the panel's rank-64 update. Each element
  is loaded and stored once per panel and sees its updates in pivot order,
  so the factors equal the pivot-by-pivot algorithm's bit for bit. The pivot
  reciprocals are stored for the solve.
- **solve**: forward substitution with the unit lower factor, then backward
  substitution. In fast mode (the default) a row with twelve or more terms
  accumulates four partial sums and the back substitution multiplies by the
  stored reciprocal; exact mode (`make KERNEL_MODE=exact`) keeps MICM's
  single sum and division.
- **forcing at y + a k**: the second stage's forcing with its combination
  formed block by block into a thread-local buffer (one fused multiply-add
  per element, as the library's pass), so no group-wide pass precedes it.

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
| group size | 32 | `make GROUP_SIZE=`; `--group-size` | cells per state, sharing one step size and error norm; blocks run in lock step, so a group fetches the kernels' code once |
| threads | all cores | `--threads`, `OMP_NUM_THREADS` | OpenMP threads over groups |
| schedule | dynamic, 1 | `--schedule static|dynamic|guided`, `--chunk` | how groups are handed to threads |
| LU panel | 64 | `make GEN_FLAGS="--panel N"` | pivots per generated LU panel; the results do not depend on it |
| solve row budget | 64 | `make GEN_FLAGS="--row-budget N"` | terms per block loop in the generated solves; short rows share a loop |
| kernel mode | fast | `make KERNEL_MODE=exact` | fast: split sums and pivot reciprocals in the solve; exact: MICM's order bit for bit |
| species ordering | minimum fill (Markowitz in exact mode) | `make TS1_ORDERING=markowitz` | the TS1 kernels' pivot order; the library's own default is MICM's Markowitz |
| kernel chunk | 16 | `make KERNEL_CHUNK=` | pivots or rows per generated function (compile time only) |
| kernels | generated | `--generic` | the benchmark runs the generic index-list kernels instead |

Every value gives correct results; the gate holds for each (RESULTS.md).
Changing W, the group size, or the kernel mode changes the roundoff and the
step-size history as MICM's own layouts and groupings do, within 1e-13 of
MICM on TS1.

## Tier 3: the scenario solver and the integrator loops

Tier 3 adds `bin/ts1_scenario` (`bench/ts1_scenario.c`), two integrator
loops beside MICM's, and the pieces they need. The Tier 1 path is untouched:
`--method micm3`, the default, runs MICM's three-stage Rosenbrock in MICM's
loop (`tmc_rosenbrock_solve`) with MICM's parameters and tolerances, and the
MICM test suite and the TS1 gate still exercise exactly that code. Only the
stage formulas were moved into a function both loops call
(`tmc_rosenbrock_stages`, with the same operations in the same order), and
the step's pieces (forcing, factorization, triangular solve) lost `static`.

```
bench/ts1_scenario.c     the scenario driver: cells, f(t), snapshots, methods and options
bench/csv_table.h        the CSV reader both drivers use
src/row.c                the Rosenbrock-Wanner loop with the Tier 3 options
src/extrapolation.c      extrapolated linearly implicit Euler
src/process.c            tmc_solver_update_custom_rate_constants: only the rate constants that read custom parameters
prototypes/              the same methods in Python over the library's pieces (prototypes/README.md)
```

### The scenario driver

    bin/ts1_scenario [--threads N] [--method NAME] [--option KEY=VALUE ...] SCENARIO

- Every cell starts from TS1's InitState, takes its CELLS_CSV row, and has
  its air density recomputed; the rate constants are computed once, before
  the timer. Cells form groups of `group_size` (default 32), each a state
  with its own step size, handed to OpenMP threads with `schedule` and
  `chunk` (default dynamic, one group at a time); each thread runs all the
  steps of its group, recording the snapshots on the way.
- At the start of step k, f(k dt) is interpolated linearly in the photolysis
  file; when it differs from the previous step's value, every photolysis
  parameter (a custom rate parameter whose label begins with "j") is set to
  its CELLS_CSV value times f, and
  `tmc_solver_update_custom_rate_constants` recomputes only the rate
  constants that read custom parameters (TS1's user-defined photolysis,
  heterogeneous, and surface reactions; 165 of 547). Temperature, pressure,
  and the Arrhenius, Troe, and other computed rate constants do not change
  within a scenario, so this is the same result as MICM's full update at a
  fraction of its cost. The rule is mechanism-generic: a reaction's rate
  constant is recomputed when an input it reads changed. The update also
  records the largest relative change of any rate constant, compounded
  since the state's last factorization, for the reuse rule below.
- The timer covers the steps, the photolysis updates, and the snapshot
  copies, and nothing else. The outputs follow `../scenarios/README.md`.
- Unknown methods, options, and option values are errors, so a mistyped
  operating point cannot run with a default.

### Methods

| `--method` | Loop | Formula | Stages | Forcing evaluations | Order (embedded) |
|---|---|---|---|---|---|
| `micm3` (default) | MICM's | MICM's three-stage Rosenbrock | 3 | 2 | 3 (2) |
| `ros2` | ROW | MICM's two-stage Rosenbrock (ROS2, a W-method) | 2 | 2 | 2 (1) |
| `ros3` | ROW | MICM's three-stage Rosenbrock | 3 | 2 | 3 (2) |
| `ros4` | ROW | MICM's four-stage Rosenbrock | 4 | 3 | 4 (3) |
| `rodas3` | ROW | MICM's four-stage DAE Rosenbrock (RODAS3) | 4 | 3 | 3 (2) |
| `rodas4` | ROW | MICM's six-stage DAE Rosenbrock (RODAS4) | 6 | 6 | 4 (3) |
| `extrap` | extrapolation | linearly implicit Euler extrapolated to `order` | | see below | `order` (`order` - 1) |

The ROW loop (`src/row.c`) runs MICM's stage formulas in MICM's transformed
form: one factorization of alpha I - J per step, with alpha = 1 / (H gamma),
and one triangular solve per stage. Each step's forcing and matrix come from
the fused `forcing_jacobian` kernel when the kernels are attached.

The extrapolation loop (`src/extrapolation.c`) takes, for each entry n_j of a
step-number sequence, n_j linearly implicit Euler substeps of H / n_j,
y += (n_j / H I - J)^-1 F(y), with the Jacobian evaluated once per step at
the step's start (the kernel's alpha = 0 call, kept in `jacobian_copy`) and
copied and shifted for each n_j; the first substep of every entry reuses
F(y_0). The rows are combined by Aitken-Neville extrapolation in powers of H,
T_jl = T_j,l-1 + (T_j,l-1 - T_j-1,l-1) / (n_j / n_j-l - 1), the result is the
last row's last entry and the error estimate its difference from the entry
before it. With the harmonic sequence 1, 2, 3, 4 (`order=4`) a step costs 4
factorizations, 10 solves, and 7 forcing evaluations.

### Options

| `--option` | Default | Meaning |
|---|---|---|
| `rtol=R` | MICM's 1e-6 | relative tolerance |
| `atol=A` | MICM's 1e-3 mol m-3 | absolute tolerance, every species |
| `group_size=G` | 32 | cells per state |
| `schedule=S`, `chunk=K` | dynamic, 1 | the OpenMP schedule over groups |
| `generic=1` | 0 | the generic index-list kernels instead of the generated ones |
| `warm_start=0|1` | 0 for `micm3`, else 1 | start each solve from the step the state's previous solve proposed, instead of MICM's 1e-6 of the interval |
| `h_start=S`, `h_max=S` | 0 | the first step and a cap on the step, in seconds (0: MICM's) |
| `norm=rms|cellmax|max` | rms | the error norm: MICM's root mean square over the group; the largest of the cells' root mean squares; the largest scaled error |
| `controller=micm|predictive|h211b` | micm | MICM's bounded safety-factor rule; Gustafsson's predictive controller, the smaller of the elementary and predicted factors as RODAS and RADAU5 take it; Soderlind's H211b filter |
| `clamp=1` | 0 | clamp negative concentrations to zero after every accepted internal step (MICM clamps once, after the solve) |
| `quantize=1` | 0 | round every step to a divisor of the solve interval, so steps tile it and repeat from one interval to the next |
| `reuse_age=N` | 0 | a factorization may serve up to N accepted steps, including across solves, while the step is unchanged and the rate constants have not changed by more than `reuse_change` (default 0.1); valid without loss of order for W-methods (ros2, ros34pw2) and for the linearly implicit Euler substeps of `extrap` (which then keeps one factorization per sequence entry in the state), an approximation for the other tableaus |
| `freeze_band=F` | 0 | keep the step when the controller's factor lies in [1/F, F], which keeps the factorization |
| `max_steps=M` | MICM's 1000 | attempts per solve |
| `safety=D`, `fac_min=Q`, `fac_max=Q` | MICM's 0.9, 0.2, 6 | the controller's safety factor and bounds on the step ratio; Dreger et al. (2025) found the embedded estimate conservative and `safety` 1.3 to 1.5 worth 10 to 30 % of the evaluations |
| `h211_b=B`, `h211_k=K` | 4, the estimator order | H211b's filter parameters; Dreger et al. used b = 1, k = 1.7 |
| `h_start_factor=F` | MICM's 1e-6 | the first step of a solve as a fraction of the interval, when not warm-started |
| `filter_error=1` | 0 | ROW loop: the error estimate is passed through the step's own matrix, (alpha I - J)^-1 alpha e, before the norm (one more solve per step), RADAU5's device against estimates excited by a stiff transient |
| `order=K`, `sequence=harmonic|shifted|even` | 4, harmonic | extrapolation: the number of entries, and the sequence 1,2,3,...; 2,3,4,...; or 2,4,6,... |

Every option value gives a valid run; the operating points in
OPERATING_POINTS.md name the ones the measurements support.

### Adding a method

- **A Rosenbrock or W tableau:** a constructor in `src/parameters.c` filling
  `tmc_rosenbrock_parameters` in MICM's form (the transformed coefficients
  `a`, `c`, `m`, `e`, the diagonal `gamma`, the new-forcing flags, the
  estimator order), a name in `method_parameters` of `bench/ts1_scenario.c`
  and in `proto_tableau` of `prototypes/shim.c`, and a case in
  `prototypes/agreement.py`. The prototype reads the coefficients from the
  library, so they are transcribed once.
- **Another loop:** a `tmc_<name>_solve` in its own file with the signature
  of `tmc_row_solve`, a `TMC_LOOP_*` value, a case in `tmc_solver_solve`
  (`src/builder.c`), any work arrays in `tmc_solver_create_state`, its
  Python mirror in `prototypes/methods.py`, and a case in the agreement test.
  The step-size rules (`tmc_row_initial_step`, `tmc_row_step_factor`,
  `tmc_row_quantize`, `tmc_row_error_norm`) are shared.
- **Another mechanism:** nothing in the loops is TS1's. The driver finds the
  photolysis parameters by their labels; a mechanism whose photolysis labels
  differ needs the prefix changed in one place.

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
