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
  instructions per block, 5 % slower); and the LU's trailing update in row
  pairs sharing their loads of the panel's row values (bit-identical, 0.6 %
  more instructions, the same time at one thread); and larger kernel chunks,
  `KERNEL_CHUNK` 64 and 128 against 16 (1.3 % fewer instructions from fewer
  calls per block, the same time at one thread, 2 to 3 % slower at ten). The
  one small pass fusion
  that measured no worse, the second stage's forcing at y + a k formed block
  by block in a thread-local buffer instead of a group-wide pass, stays:
  0.5 % at one thread, bit-identical.

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

# Tier 3

Decisions for the numerical-method tier, in the order they were taken. The
Tier 1 decisions above stand.

## Scope and rules

- **Starting point and record.** This directory is the Tier 1 entry at
  commit 6d2fcde with its history, as the Tier 3 prompt asks. The author of
  this entry also orchestrates the experiment and had seen the other Tier 3
  entries' verification reports (their methods' names, operating points,
  and scores) before this work began, but not their code or documents; the
  methods here are chosen from the literature and from this entry's own
  scans, and the record says so rather than pretend otherwise.
- **The Tier 1 path is frozen.** MICM's loop, parameters, tolerances, and
  tests stay as they are; `--method micm3` runs them. New behavior lives in
  new loops selected by a parameter, with every new parameter zero meaning
  MICM's behavior, so the 680 cases and the TS1 gate test the same code as
  in Tier 1.
- **The literature was read before the methods were written.** Six
  focused literature searches (Rosenbrock and W-methods; accuracy
  requirements and solver comparisons; extrapolation, exponential, and
  SDIRK methods; controllers, norms, reuse, and positivity; alternative
  schemes and precision; the TS1 mechanism and CAM-chem's solver) were run
  by research assistants and then checked source by source; REFERENCES.md
  lists what was used and for what. No code was copied from any source.

## The scenario driver

- **Only the photolysis rate constants are recomputed each step.** In a
  scenario, f(t) changes the photolysis parameters and nothing else, so
  `tmc_solver_update_custom_rate_constants` recomputes the rate constants
  that read custom parameters (165 of TS1's 547) and leaves the Arrhenius,
  Troe, and other computed ones, which depend only on the unchanged
  temperature, pressure, and air density. The results are identical to a
  full update. The rule is generic: recompute what reads an input that
  changed. The full update runs once per group before the timer, as the
  prompt allows ("not setup").
- **A step with an unchanged f(t) skips the update** (at night f = 0 on
  every step). The count of updates is printed (`rate_updates=`).
- **Groups run their whole scenario on one thread.** Each group's steps
  depend on its previous state, so the parallel loop is over groups, each
  thread taking its group through all the steps and recording snapshots
  itself; no barrier per step. The schedule is dynamic, one group at a
  time, as in Tier 1.
- **Mistyped options fail.** An unknown method, option, or value is an
  error, never a silent default, so an operating point runs as declared or
  not at all.

## The loops

- **The ROW loop re-uses MICM's stage formulas.** MICM's transformed
  Rosenbrock form (one factorization of alpha I - J per step, c_ij / H
  terms instead of matrix products) serves every Rosenbrock and W tableau,
  and the generated kernels (fused forcing and Jacobian, forcing at
  y + a k, LU, solves) serve every method unchanged. The stage loop was
  moved into a function MICM's loop also calls, with the operations in the
  same order, so Tier 1's bit-for-bit agreement with MICM holds (checked
  by the tests).
- **The step is carried from one solve to the next** (`warm_start`,
  default on for the new loops). MICM restarts every solve at 1e-6 of the
  interval and needs about eight growth steps to recover, which is why it
  takes nine internal steps per 30-second Tier 1 step whatever the
  tolerance; with 5 to 30 minute scenario steps that restart would dominate
  the cost. The carried step is the controller's proposal, not the last
  step, which the end of the interval may have truncated.
- **Factorization reuse is a rule on the state, not on the method.** A
  factorization is reused when the step is the same, its age is below
  `reuse_age`, and the rate constants have changed by less than
  `reuse_change` since it was made; a rejection or a changed step ends it.
  W-methods keep their order under it (ROS2 is one); for other tableaus it
  is an experiment the scans measure. The rate-change test is cheap (it
  falls out of the photolysis update) and catches dawn and dusk, when the
  photolysis rates change by large fractions from one step to the next.
- **Steps that tile the interval** (`quantize`) round the step to
  dt / n. With a carried step this makes the step sequence repeat from one
  scenario step to the next, which is what lets a factorization live across
  scenario steps, and it removes the truncated last step that would
  otherwise force a new factorization every interval.
- **Extrapolation evaluates the Jacobian once per step and shifts a copy
  for each substep size.** The alternative, re-assembling alpha I - J for
  each n_j, would cost a Jacobian assembly (438 ns per cell) per entry
  instead of a copy and a diagonal add.
- **Clamping inside the loop is an option, off by default.** MICM clamps
  once after the solve; clamping after each accepted step changes the
  trajectory slightly and is measured, not assumed.

## What the scans settled (in the order they ran)

- **S1, the first step:** the carried step with quantization is the
  default for the new loops (35 % less time than MICM's restart at equal
  accuracy; SCANS.md). The scan also showed that at rtol 1e-3 the loop
  rejects a quarter of its attempts once the step is carried, because
  MICM's controller lets the step grow six-fold after every success; the
  controller scan addresses it, and a filtered error estimate
  (`filter_error`, one extra solve, RADAU5's first-step device) was added
  as an arm for it.

- **S2, the tableau:** RODAS3 is the base method for both levels, as the
  1997 benchmarks found for one to three digits: at the tolerances that
  reach each level on all six disclosed scenarios it is the fastest, with
  ROS3 and ROS4 close. Extrapolated linearly implicit Euler is kept as the
  method for the accurate end (from about 3.5 digits it is the fastest, and
  it reaches four digits where RODAS3 would need rtol 1e-4 and many steps),
  not for the ranking levels: at loose tolerances it takes one large step
  per interval whose error its estimate does not see on the 30-minute and
  polluted scenarios. ROS2, RODAS4, and ROS34PW2 are kept as options and
  used in no operating point.
- **S6a, factorization reuse:** not pursued for the operating points. By
  day the photolysis rates change by more than a few percent every interval,
  so a kept factorization lasts an interval or two; the extrapolation gains
  18 % and loses up to half a digit, RODAS3 (not a W-method) loses its
  order, and ROS34PW2 with reuse remains slower than RODAS3 without. The
  mechanism stays in the code, tested, for rates that change more slowly.
- **S4, the controller:** H211b with b = 4 and k = p at safety 0.9 is the
  candidate for RODAS3, because it is the one controller above the
  elementary controller's work-precision curve (7 to 17 % less time at
  equal worst-species digits, from a smoother step sequence with half the
  rejections and a third more steps per tolerance). The predictive
  controller only trades rejections for smaller steps; Dreger's (1, 1.7)
  setting acts like the elementary controller on this problem; safety
  above 1 and a larger growth limit gain nothing; the filtered error
  estimate lands on the same curve as a looser tolerance and costs a solve
  per step. S12 confirmed it under the final conventions (quantized step,
  absolute tolerance 1e-20), and it is the controller of the balanced point;
  the elementary controller stays for extrapolation and, by S13b, for the
  level-2 point.
- **S5, the norm:** MICM's group RMS norm stays. The per-cell and maximum
  norms are the same curve at a rescaled tolerance (by 2.5 and 20), and
  controlling the worst component directly does not help the species that
  sets the score.
- **S3, the absolute tolerance:** 1e-18 for the level-2 points, where the
  gain saturates; for level 3 the absolute tolerance is a lever (at rtol
  5e-4, 1e-20 and 1e-22 each add 0.2 digit for the worst species for 10 to
  15 % more time, cheaper than the same digits through rtol), because the
  species that set the score live near the 1e-14 floor the score counts
  from. The level-3 points come from the joint grid S12.
- **S7, clamping:** off, as in MICM; every run passes the gate without it,
  and it changes nothing measurable.
- **S9, grouping:** 32-cell groups and dynamic scheduling, as in Tier 1;
  8 to 32 are equal within noise at 10 threads, 64 and 128 cost 15 to
  30 %, static scheduling 5 to 10 %.
- **S8, extrapolation sequences:** the shifted sequence (2, 3, 4, ...
  substeps) replaces the harmonic one for the extrapolation points: it
  skips the single-substep entry whose error the estimate misses, gains 0.2
  to 0.4 digit for the worst species at equal tolerance, and dominates on
  the frontier. Orders above 4 are not used (the step cannot grow into
  them within an interval), and the even sequence at orders 5 and 6 gave
  non-finite output on the 30-minute scenario.
- **S9 at one thread:** 32-cell groups confirmed (64 equal, 16 costs 40 %
  more, 8 up to 90 % more).
- **S12, the joint grid:** the absolute tolerance is 1e-20 for every
  operating point but the accurate end's alternative at 1e-22. At rtol
  1e-2 the phenol floor of 1e-18 (2.06 digits on diurnal-10min) lifts to
  2.5 with the elementary controller and 3.0 with H211b; at level 3 the
  tighter absolute tolerance buys digits for the worst species more cheaply
  than the relative one.
- **S13, the confirmation on all six scenarios:** the operating points come
  from this scan, under ACCURACY.md's rule (declared for a level only with
  0.2 digit of margin for the worst species on every disclosed scenario).
  Level 3 goes to extrapolation of order 3 with the shifted sequence at
  rtol 1e-3 and atol 1e-20 (3.61 digits, 2.49 for the worst species, in
  2.65 s over the six scenarios), ahead of RODAS3 with H211b at 1e-3 / 1e-18
  (3.29 and 2.32 in 2.58 s) on margin and of the elementary controller at
  5e-4 / 1e-20 (3.39 and 2.35 in 2.95 s) on time; the recommended point is
  order 4 at 5e-4 / 1e-22 (3.62 and 2.69 in 3.04 s), and the accurate point
  order 4 at 2e-4 / 1e-20 (3.87 and 3.00 in 3.45 s). RODAS3 with H211b at
  2e-3 / 1e-20 is declared for level 2 as the balanced point: it reaches
  three digits with 0.17 of margin on all six (3.17 and 2.17 in 2.21 s),
  which the ranking credits at level 3 if it holds on the undisclosed
  scenarios, and it is the fastest configuration above 2.5 digits. The
  loose RODAS3 points at rtol 1e-2 reach 2.04 digits (BrO 1.15) on
  diurnal-5min whatever the controller, the contest's floor; they are
  reported, not declared. The level-2 point is the elementary controller
  at 5e-3 / 1e-20 (2.42 and 1.29 in 1.16 s): S13b showed every loose H211b
  setting at safety 1.0 to 1.2 meeting the same 2.04-digit floor on
  diurnal-5min.
- **The order of the remaining scans** followed from the cost picture: the
  fastest configurations take one to five RODAS3 steps per interval, and
  the step count is set by the worst species (isoprene) reaching L - 1
  digits while the overall error is a digit better. So: the controller and
  the safety factor (S4), then the error norm, whose maximum form controls
  the worst species directly (S5), the absolute tolerance (S3), the group
  size and scheduling (S9), clamping (S7), the extrapolation's sequences
  (S8), and the generic kernels (S10); precision (S11) last, if the
  algorithmic levers leave it worth the generator work.

## Prototypes

- **The prototypes use the library's linear algebra.** A Python mirror of
  each loop over the library's forcing, Jacobian, and sparse LU (through a
  small ctypes shim) agrees with C to roundoff and counts the same
  operations, which is what the prompt's agreement test asks; a prototype
  with its own dense LU would differ at the level of the factorization's
  rounding on a stiff matrix and prove less. The shim and the shared
  library are built apart from the benchmark (`make proto`, objects with
  -fPIC in build/proto), so the benchmark's build is unchanged.
