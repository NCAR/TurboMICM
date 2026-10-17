# Decisions

Choices the prompt leaves open, with the reasoning. Numbers are from this
machine (Apple M5, clang 21); RESULTS.md has the measurements.

## Scope and rules

- **Linear algebra.** The experimenter allowed changes to the linear algebra.
  The default build keeps MICM's algorithm exactly where it costs nothing
  (Mozart elimination order, Markowitz species ordering) and departs in two
  places, both switchable: the generated back substitution multiplies by the
  pivot reciprocals the factorization already computes, and sums with twelve
  or more terms accumulate in four partial sums. `make KERNEL_MODE=exact`
  restores MICM's order bit for bit, and the tests run that mode against
  MICM's reference result (the difference is zero) and the fast mode against
  the gate (largest difference 1e-13 on 10,000 cells).
- **Species ordering.** The library's default is MICM's diagonal Markowitz
  ordering, reproduced exactly, including that MICM computes the Markowitz
  cost in unsigned arithmetic, so a species with no incoming edge wraps
  around and is ordered last, and that the pattern it orders has no forced
  diagonal. The fast TS1 kernels are generated for a greedy minimum-fill
  ordering instead (`TMC_ORDERING_MINIMUM_FILL`, `tools/mechanism.py` and
  `tmc_minimum_fill_reorder` are the same algorithm: the pivot adding the
  fewest elements, ties by the Markowitz cost and then the index), which
  attached kernels impose on the solver they are built into: 2,293 LU
  elements against 2,381, 10 % fewer elimination multiply-adds and 4 % fewer
  in the solves, 3.8 % fewer instructions, and 2.7 % less time at one
  thread, 3.3 % at ten. The results move from 3.8e-14 to 4.5e-14 of MICM
  on the 100-cell input. Exact mode keeps MICM's ordering (`make
  KERNEL_MODE=exact`, or `TS1_ORDERING=markowitz`), so its results still
  match MICM bit for bit; so does the benchmark's generic path.
- **Everything else numerical is MICM's:** the Rosenbrock formulas and all five
  parameter sets, the error norm's accumulation order (variable by variable,
  cell by cell, padded lanes excluded), the step-size controller, the default
  tolerances (relative 1e-6, absolute 1e-3 unless a species sets one), the
  rate-constant formulas and the order of their factors, and the post-solve
  clamp.

## Layout and grouping

- **Vector width 2 (`TMC_WIDTH`, `make WIDTH=`).** Every per-cell array is
  vector ordered: blocks of W cells with the cell fastest. On this machine's
  NEON units a vector instruction covers two doubles, so wider blocks only
  enlarge the generated code (which already exceeds the 192 KB instruction
  cache at W = 4) without doing more work per instruction. W = 4 or 8 suits
  AVX2 or AVX-512 machines; W = 1 is MICM's standard ordering. Every width
  gives correct results; the generated kernels are compiled for one width and
  refuse another.
- **Groups of 32 cells (`GROUP_SIZE`, `--group-size`).** A group is one
  state with one step size and one error norm, as MICM's OpenMP tutorial
  does per thread. The kernels run a group's blocks in lock step: each loop
  body of straight-line code is fetched once per group and reused for every
  block, and the blocks' dependency chains overlap. The generated code
  (500 KB) is far larger than the 192 KB instruction cache, so the fetch is
  what small groups pay for: at one thread, groups of 4 cells run at 2.4
  instructions per cycle, 8 at 2.8, and 32 at 3.7. In the final timing
  session (RESULTS.md) 32-cell groups took 26.4 s against 33.3 s for 16 and
  about 32 s for 8 at one thread, and 4,850 ms against 4,917 (16) and
  5,030 (8) at 10 threads; 64 are slower again at one thread. The error norm
  is a root mean square over the group, so small groups are stricter per
  cell; from 4 cells up TS1 takes MICM's 9 internal steps per 30 s step in
  every group.
- **Scheduling.** The benchmark hands groups to threads with
  `schedule(dynamic, 1)` by default (`--schedule`, `--chunk`): the four
  performance cores and six efficiency cores differ in speed, and 1,250
  groups make the tail short. Static scheduling costs 8 % at 10 threads.

## Generated code

- **What is generated.** For any mechanism, `tools/generate.py` emits the
  species and reactions as calls to the library's builder API, the recorded
  initial state, and four kernels with every species and matrix index a
  literal: forcing, the shifted Jacobian alpha I - J, the in-place LU, and the
  two triangular solves. The library attaches them only if their width and a
  fingerprint of the species order, reactions, yields, and LU pattern match
  the solver it built, so a configuration change cannot run stale kernels.
  Everything else (rate constants, the Rosenbrock loop, the error norm, the
  backward Euler and constraint paths) is library code shared with the
  generic kernels.
- **Gather form.** The forcing computes every reaction rate once, then each
  species sums its terms; the Jacobian computes every (independent variable,
  reaction) partial once, then each matrix element sums its terms and the
  diagonal adds alpha last. Each element is written once, with its terms in
  MICM's accumulation order. The four longest sums are emitted side by side
  (interleaved statements) so their dependency chains overlap.
- **LU in panels, every element finished target by target.** The Mozart
  elimination runs in panels of 64 pivots. Inside a panel, each pivot's
  diagonal, column (then the scaling by the reciprocal), and row take their
  updates from the panel's earlier pivots in one load, a chain of fused
  multiply-adds in pivot order, and one store; the trailing elements then
  take the panel's rank-64 update the same way. Each element still sees its
  updates in pivot order, so the result is bit for bit the pivot-by-pivot
  algorithm's (the 100-cell result is byte for byte the same for every panel
  size). Against panels of four pivots whose own rows and columns were updated
  statement by statement, 64-pivot panels cut the LU from 31,900 to 24,500
  instructions per block (12,600 loads and 5,800 stores to 11,300 and 2,900)
  and the whole run by 9 % at one thread; the full matrix as one panel is
  slower again (`--panel`).
- **Measured rejects.** Variants that measured slower on this machine and were
  dropped: a left-looking (row-wise) LU (long serial sums, 2 times slower), a
  block loop around each pivot instead of each panel (loop overhead), block
  loops through a table of base pointers (defeated the compiler's addressing),
  a constant table for the yields (no change), scattering forcing and
  Jacobian terms into memory as MICM does (store chains on the dense species),
  skewed array addresses against 4 KB aliasing (3 % slower), and block
  offsets the compiler cannot analyze, which remove its loop prefetches
  (2.8 % fewer instructions, 3 % slower: the prefetches hide the L2 latency
  of a block's data), and a generated solve that forms each row's right-hand
  side on the way (f + c0 K0 + c1 K1) instead of the stage combination pass
  and the copy of the initial forcing: bit-identical, but three forward
  substitutions in the code (76 KB more), one more instruction per 150, and
  1.5 % slower at one thread, 3 % at ten; the solution, error estimate,
  and error norm in one variable-major pass (bit-identical, 3 % more
  instructions, the same time at one thread and 2 % slower at ten); not
  writing the LU's fill positions in the Jacobian assembly, which the LU
  then starts from zero (bit-identical, 570 fewer instructions per block,
  but 25 % slower at one thread); and the pivot reciprocals kept in the
  diagonal instead of their own array (bit-identical, 140 fewer
  instructions per block, 5 % slower).

## Configuration format

- YAML is read with PyYAML, whose default schema is YAML 1.1 and reads
  `2.3e2` as a string; the reader adds the YAML 1.2 float forms, so the
  examples' numbers are numbers whether written as `2.3e2` or `2.3e+02`.

- Reactants and products may be CAMP maps (`{"A": {"qty": 2}}`) or ordered
  lists (`["A", "A"]`, `[{"species": "B", "yield": 0.5}]`), an extension that
  keeps the reactant order the C++ declarations use, which fixes the order of
  the floating-point products.
- `"rate constant units": "mol m-3"` on a mechanism marks parameters already in
  MICM's internal units, as TS1's are; without it the reader converts from
  molecule cm-3 the way MICM's tutorial programs do by hand (pre-exponential
  factors scaled by (1e-6 N_A) to the power of the reaction order minus one,
  and falloff parameters likewise).
- A species with `"third body": true` is parameterized on air density, like
  MICM's `SetThirdBody()`; it holds no state variable and multiplies into the
  rate constant.
- `"initial concentration [mol m-3]"` on species, and top-level
  `"initial conditions"` and `"initial rate parameters"`, record TS1's
  `InitState`, so the benchmark's defaults come from the configuration.
- PHOTOLYSIS, EMISSION, and FIRST_ORDER_LOSS become user-defined rate
  constants labelled by their `MUSICA name`, as MICM's by-hand tutorials build
  them; a SURFACE reaction's label is its `MUSICA name` ("C surface" in the
  example, where the by-hand program uses "C").

## Tests

- MICM's cases are ported one function per case id and printed in the
  reference list's order (tests/README.md). Two checks of this port's own run
  after them (`test_cpu --extra`, to stderr, so the 680 lines stay exactly
  the list): the generated TS1 kernels against the generic kernels on the
  reference input, and the experiment's gate on the 100-cell reference input.
- Randomized MICM tests use fixed seeds in C; the checks compare against an
  independent computation, so the values drawn do not change what is tested.
- Internals MICM's tests reach (the error norm, the alpha shift, the backward
  Euler convergence test, constraint terms on their own, array functions, name
  maps) are exposed in `include/tmc_*_ext.h` and `tmc_support.h`, and they
  call the solvers' own functions, so the tests check the code that runs.

## C API

- Species, reactions, rate constants, and phases are value types with
  fixed-size names (128 characters) and property lists (16 per species), so
  they copy without allocation; reactions and phases own small arrays and
  have `_free` functions.
- MICM's template parameters become run-time values: the vector width and
  group size per solver and state, the LU algorithm an enumeration, the
  Rosenbrock parameters a struct (five constructors), the step-size controller
  a function pointer, and the kernels a table of function pointers.
- Errors are returned, never aborted on: functions that can fail take a
  `tmc_error *` and return its category; the codes follow MICM's error.hpp.
  Out-of-memory aborts, since nothing useful can follow.
- MICM's defaults and rules hold where its tests observe them: the builder
  ignores unused species unless told otherwise; a single concentration or
  parameter value is for a one-cell state (the `tmc_state_fill_*` functions,
  which MICM lacks, set every cell, and the generated initial states use
  them); the post-solve clamp is `std::max(y, 0)`, which keeps a NaN visible;
  a reaction has at most four parameterized reactants, whose product is
  accumulated from 1 before it scales the rate constant.
- An external model's context belongs to the caller and is bound to one
  solver when it is built; a builder reused for another configuration needs
  another context. MICM copies the model object instead; C has no copy of an
  opaque context without another callback, and no case needs one.
