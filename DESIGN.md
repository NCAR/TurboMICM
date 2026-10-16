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
mechanisms; external-model and constraint support is still under development.

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

`chem_lu_algorithm` selects Mozart or Doolittle operation generation. The current
factor and solve entry points are isolated in `src/linear.c`; a pluggable backend
interface is still to be added for replacing the linear solver.

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
