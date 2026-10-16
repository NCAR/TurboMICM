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
| Group scheduling | static contiguous blocks | currently fixed; no speed-specific chunk constant |
| SIMD width | compiler-selected | portable `omp simd`, no hardware-specific width |
| Chemistry evaluator | automatically generated | `--generic` for the generic table evaluator |
| LU method | Mozart in benchmark | library `chem_lu_algorithm` |

There are no target-specific flags, intrinsics, vector extensions, or alternative
floating-point modes. Group size affects work distribution and the common error
norm, not the method or requested tolerances. Every positive group size and
thread count must remain correct; verification is still being expanded.

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
