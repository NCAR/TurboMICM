# Extension exercises across the nine Tier 1 entries

You are helping the second judge of TurboMICM's Tier 1 assess **extensibility**,
one of three secondary criteria (readability, maintainability, extensibility).
The Tier 1 design says the reviewers "try to add a Rosenbrock parameter set,
replace the step-size controller, swap the linear solver, and add a
rate-constant type"; the judge adds a fifth: use another mechanism. The SSE (the
scientific software engineer who runs the experiment) chose that the exercises
be done **on paper**: you work out, from the code as it stands, exactly what
each change takes. You do the same exercise for all nine entries, so that one
standard applies to all.

## The contest, in short

Each entry is one OpenAI model at one reasoning effort that ported MICM's CPU
chemistry solver (C++) to ISO C11 with libm and OpenMP, kept MICM's numerical
methods (Rosenbrock with MICM's parameter sets and controller, backward Euler,
sparse LU), ported MICM's 680 test cases, reads MICM's YAML/JSON
configurations with a Python reader and generator, and runs the TS1 benchmark
(210 species, 547 reactions). Rules: `~/NSF-NCAR/TurboMICM.wiki/Tier-1-Design.md`.

## The entries (read-only)

Each is a git worktree at the entry's final commit:
`~/NSF-NCAR/tier1-work/verify/<run>/clang/<run>/`, for `<run>` in
sol-high, sol-xhigh, sol-max, astra-high, astra-xhigh, astra-max, luna-max,
luna-xhigh, luna-high. Their README.md and DESIGN.md say how each is meant to
be extended. A description of each entry's architecture, with citations, is in
`~/NSF-NCAR/TurboMICM.wiki/Tier-1-Implementations.md` (use it to
find your way; check everything against the code). MICM is at
`~/NSF-NCAR/ODE/benchmarks/micm/` for comparison.

Read-only: do not edit, build, or run anything in these directories. Light
commands (`grep`, `sed -n`, `wc`, `git show`, reading files) are fine; no
builds, no tests, no benchmark.

## The exercises

- **E1. Add a Rosenbrock parameter set.** A new four-stage method with its own
  coefficients (a, c, m, e, alpha, gamma), order, and new-function-evaluation
  flags, selectable at run time beside MICM's sets; for definiteness, RODAS3
  (Sandu et al., 1997). Where are parameter sets defined and selected? Is
  there a maximum stage count? Does the integrator that TS1 uses loop over any
  stage count, or does it hard-code three stages? Would generated code or the
  benchmark need changes? What test would cover the new set?
- **E2. Replace the step-size controller.** Replace MICM's controller (step
  factor 0.9·err^(−1/(order+1)) clamped to [0.2, 6], with MICM's rules after
  rejections) by a PI controller, such as Gustafsson's or Söderlind's H211b,
  which needs the previous step's error and step size. Is there a controller
  hook, and what does it receive? Where are the error norm and the
  accept/reject decision? Is the controller logic in more than one code path
  (grouped and single-cell paths, Rosenbrock and backward Euler)? Where would
  the controller's state live between steps and between solver calls?
- **E3. Swap the linear solver.** Replace the sparse in-place LU in the
  Rosenbrock stages by another solver (for example a dense LU, or an external
  library's sparse solver), chosen at run time, without editing the
  integrator. Is there a linear-solver interface (function pointers, an enum),
  and what does the solver receive (matrix storage, sparsity, group width)?
  How many call sites does the integrator have? Do generated Jacobian kernels
  or the species ordering assume the current storage? What would the tests
  check?
- **E4. Add a rate-constant type.** A new type with its own parameters that
  depends on temperature and pressure, for definiteness
  k = A·exp(C/T)·(T/300)^B·(1 + D·P/101325), added end to end: MICM-format
  configuration, Python reader (and unit conversion), generator output, C data
  structures, the per-cell rate evaluation, and a test. Count the places that
  must change, in Python and in C. Does a callback type for user rates exist,
  and would it serve instead?
- **E5. Use another mechanism.** Run a new mechanism in MICM's configuration
  format (one of MICM's tutorial configurations, or a successor of TS1 such as
  TS2) through the entry: what a user runs, what must be edited by hand,
  whether a generated-kernel path supports it, whether the tests run the
  generator on any mechanism besides TS1, and any defect that would give wrong
  results (one is known: in luna-high, `tools/read_config.py:83-86` gives each
  species the index of its input position and `tools/generate_c.py:100-137`
  emits it beside a species table in Markowitz order, so reactions name the
  wrong species when that order differs from the input's; check it, and look
  for similar hazards in every entry).

## For each entry and exercise

- **Edit sites:** every file and function that must change, with
  `path:start-end`, and the count of distinct sites.
- **Extension point:** the existing run-time hook used, if any (table entry,
  enum value, function pointer, struct field), with its citation.
- **Hazards:** hidden limits (array sizes, maximum stages), logic in two
  places (Python and C, or two integration paths), generated code to change,
  coverage the tests would lack.
- **Documentation:** whether DESIGN.md or README.md says how to do it, and
  whether that matches the code.
- **Patch sketch:** at most 12 lines of pseudo-diff for the essential change.
- **Ease, 1 to 5:**
  - 5: one place, through an existing run-time extension point, documented or
    obvious;
  - 4: a few edits in one module, no hidden coupling;
  - 3: several files, or a hidden limit or coupling (a fixed maximum stage
    count, logic duplicated between Python and C or between code paths);
  - 2: invasive: the path TS1 uses hard-codes the current choice, or
    generated code must change too;
  - 1: needs restructuring.

When you have done one exercise for all nine entries, put the nine ratings
side by side, re-read the evidence, and make them consistent: the same kind
of change gets the same rating everywhere. Then order the nine for that
exercise, ties allowed.

## What to hand back

A report at the path your task names, in this order: for each exercise, a
summary table (entry, ease, edit-site count, extension point, main hazard),
then the per-entry details as listed above. Facts and citations; no words
about the agents or models beyond the run names. Your final message: only the
summary tables, plus any defect you found and how sure you are of it.
