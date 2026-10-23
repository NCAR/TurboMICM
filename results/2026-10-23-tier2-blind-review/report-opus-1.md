# Blind review of nine Tier 2 entries: the generators' readability, maintainability and extensibility

I read each entry's README, DESIGN, SPECIALIZATION and DECISIONS files, its generator sources, its generator tests and its Makefile. I built and ran nothing, so every exercise below was done on paper. Line references point into each entry's own files. A few problems affect almost every entry; they are listed once, under "Shared caveats" in the comparison, and each exercise (c) below refers to them.

---

## Entry A

**Summary.** `tools/specialize.py` (337 lines) reuses the Tier 1 reader in `tools/configure.py` and builds a `Plan`: the Markowitz species order (MICM's unsigned-cost tie rule), filled sparse rows, packed indices, and Mozart or Doolittle operation lists. An `Emitter` groups scalar "units" into `#pragma omp simd` cell loops and splits them into functions in a separate `_kernels.c` file, with constant-stride and runtime-stride versions. The generated code fixes the rate formulas (`rate_expression`), the forcing and Jacobian terms (in the generic library's accumulation order), the diagonal shift, the LU and both triangular sweeps. A hand-written template, `tools/specialized_solver.h`, does the three-stage integration; it keeps the library's coefficient table, ordered error norm and step controller. `make test-specialized` builds `tests/specialized.c` for six mechanisms: TS1, the three tutorials, an all-rate-family fixture and an empty mechanism (`tests/generator_cases.py`). Each run compares rates, forcing, Jacobian, LU, the linear solve, and five adaptive steps (including step and rejection counts) with the generic library, at 1e-10 relative with a 1e-30 floor. The cell counts are 1, 3, 7, 32 and 33, each at default and tight tolerances.

**Scores**
- **Readability 6.** The layout is clear, and `DESIGN.md` ("Specialized compiler") and `SPECIALIZATION.md` describe it accurately; `Plan.lu_operations` reads like the algorithm. The code itself is terse:
  - Rate families are bare integers 0–9 (`tools/specialize.py:180-225`); their names exist only in `RATE_NAMES` in `configure.py`.
  - Parameters are positional (`p[offset+2]`), and the branched-rate constants carry no comment.
  - `Emitter.kernel` builds call arguments with string tricks.
  
  The generated C has no reaction comments, but a JSON sidecar records the species order and the counts.
- **Maintainability 8.** Every generated kernel and full trajectories are compared on all six mechanisms (`tests/specialized.c:42-89`). Make regenerates when the configuration, reader, generator or options change, using `build/generator-flags` and `build/flags` stamps (`Makefile:56-111`). The reader raises `ValueError` for bad species, types, `qty` values or non-finite numbers. It silently ignores unknown keys, though, and merges custom parameters that share a MUSICA name (`tools/configure.py:97-100,157-158`). The documentation matches the code: the "8,439 Mozart operations" figure equals 209 + 1,244 + 6,986, the counts other entries give.
- **Extensibility 7.** Editing a reaction only touches `configs/ts1.json`, and the `bin/test_specialized_%` pattern rule makes adding a mechanism to the suite a two-line change. A new rate type touches five files with positional slots, and the guide in `DESIGN.md:232-241` never mentions `specialize.py`. LU variant, block size, fusion and fixed width are already options.

**Exercises**
- (a) Add a reaction to TS1 — **easy**
  - Append an object to the reactions in `configs/ts1.json`: a `reactants` list and a `products` list of `{species, yield}`. The mechanism declares `rate units: mol m-3`, so a cm³ molecule⁻¹ s⁻¹ constant must be multiplied by 6.022e17 for each extra reactant.
  - `make bench`, then `make test-specialized` (runs `bin/test_specialized_ts1`).
  - Care: a new user or photolysis label becomes a custom parameter that defaults to 0.
- (b) Change a rate parameter — **easy**
  - Edit `A`, `C` and so on in `configs/ts1.json`. Reactions have no names, so find the right one by grepping for its reactants.
  - `make bench test-specialized`.
- (c) Another mechanism — **moderate**
  - `.conda/bin/python tools/specialize.py ../mechanisms/chapman/config.json build/generated/chapman_specialized.c --name chapman_specialized` (the README gives this form).
  - Add a generation rule (copy the Robertson one, with `$(wildcard ../mechanisms/chapman/*)` as prerequisites) and add `chapman` to `SPECIALIZED_EXAMPLES`; `make test-specialized` then builds `bin/test_specialized_chapman` from the pattern rules.
  - No documentation describes this recipe. In Carbon Bond 5, the photolysis and the emission named `NO2` (and those named `ALD2`) end up sharing one parameter. The shared caveats also apply.

**Strengths**
1. Thorough kernel-level and trajectory tests on six mechanisms, including every rate family and an empty mechanism.
2. A compact generator with a clean split into shared reader, Plan, Emitter and rate emitters, plus a readable integration template that reuses the library's parameter table and controller.
3. Exact rebuilds through option and flag stamps; a JSON sidecar maps generated indices back to species.

**Weaknesses**
1. Terse code: numeric rate codes, positional parameter slots, uncommented constants.
2. The reader silently ignores unknown keys and merges parameters that share a label.
3. No recipe for adding a mechanism to the tests; the rate-type guide leaves out the generator; generated code reaches into library internals (`internal.h`).

---

## Entry B

**Summary.** `tools/specialize.py` (374 lines) imports the validated reader `tools/configure.py`. Its `bind` step applies MICM's stable rate-type order, drops parameterized species, computes the reference Markowitz order and fill, and stores explicit derivative and pivot schedules in a `Plan` dataclass. `rate_body` emits each rate family with named, fixed parameters. `Emitter.function` and `dispatch` batch work into SIMD cell loops inside bounded helper functions, with fixed-stride and dynamic-stride versions called through `src/specialized.c`. The output defines a `ck_specialized` descriptor (`include/chem_specialized.h`) holding update, forcing, Jacobian, factor, solve and integrate entry points; integration comes from the template `src/specialized_solver.h`. `make test-specialized` runs four Python structural tests (`tools/test_specialize.py`, one of which edits TS1 and regenerates) and 70 numerical cases in `tests/test_specialized.c`. The cases cover seven generated variants (TS1, three tutorials, an all-rate fixture, a zero-order fixture, and an "alternate" build with every optional transformation and odd widths), at widths 1, 3, 32, 37 and 64 and two tolerance regimes. Each case compares rates, forcing, Jacobian, LU, the solve, and six steps with every statistic, at 1e-10.

**Scores**
- **Readability 7.** Named rate types and parameters (`tools/specialize.py:126-177`), the dataclass plan, and docstrings explaining why helpers keep runtime loop bounds make the generator easy to follow. `DESIGN.md` describes it accurately: its parameter table matches the Makefile and argparse. On the minus side, `emit` is one long function, and the template hard-codes the Rosenbrock coefficients as unlabeled literals (`src/specialized_solver.h:50-104`).
- **Maintainability 9.** Testing is the strongest of the nine: every kernel, trajectories and full statistics on seven variants, including one that combines all non-default options (`Makefile:145-147`), plus ten reader bad-input tests (`tools/test_configuration.py`). Make stamps the compiler settings and generator options and lists every dependency. Minor issues: the long option list is copied into six rules, and wrapped reader errors ("Invalid configuration: 'X'") give little context.
- **Extensibility 8.** Configuration edits regenerate everything, and the README documents how to generate and link another mechanism through the descriptor. Adding a mechanism to the suite takes three edits in two files. A new rate type touches `configure.py` (TYPES, DEFAULTS, units), `specialize.py` (`PRIORITY`, `rate_body`) and the C library; DESIGN's guide leaves out the generator steps.

**Exercises**
- (a) **easy:** append to `mechanisms/ts1.json` (top-level `concentration units: mol m-3`); `make bench test-specialized`.
- (b) **easy:** edit the value; same commands. `test_edit_reaction_and_parameter` already shows the generator picks such edits up.
- (c) **moderate**
  - `.conda/bin/python tools/specialize.py ../mechanisms/chapman/config.json generated/sp_chapman.c --symbol sp_chapman` (README).
  - Add a `generated/sp_chapman.c` rule with the option flags, and add `build/sp_chapman.o` to the prerequisites and link line of `bin/test_specialized`.
  - Add `sp_chapman` to the `extern` list and to `cases[]` in `tests/test_specialized.c`; `make test-specialized`.
  - Care: Carbon Bond 5's two `NO2` user parameters get separate slots under the same name. The shared caveats also apply.

**Strengths**
1. The most complete verification, including a test of all optional transformations together and of editing then regenerating.
2. A readable generator with named formulas and explicit schedules, behind a clean descriptor API.
3. Precise documentation and exact rebuild stamps.

**Weaknesses**
1. The template repeats MICM's coefficients as literals instead of using the library's table.
2. Option lists are duplicated in the Makefile, and adding a mechanism is undocumented.
3. Labels have no PHOTO./EMIS. prefix, and the rate-type guide omits `specialize.py`.

---

## Entry C

**Summary.** The generator is split into small modules: `tools/configure.py` (reader); `tools/symbolic.py` (Markowitz order emulating MICM's unsigned cost, fill, Jacobian jobs; shared with the Tier 1 emitter `generate.py`); `tools/specialize_rates.py`; `tools/split_kernels.py`; `tools/tiled.py`; and the template `tools/ros3.c.in`. `specialize.py` emits two copies of every kernel. One is a flat set with Mozart LU, exposed for inspection. The other is a packed four-cell-tile set with Doolittle LU and reciprocal pivots, used by the solver; `tiled.py` produces it by rewriting the emitted lines and the template text. Rates, forcing, Jacobian, LU and the triangular sweeps are fixed-index straight-line code, and the three stages use the template's literal coefficients. `tests/specialized/check.c` runs seven mechanisms (TS1, C's own Chapman, three tutorials, an edge fixture and an empty mechanism) at 1, 3, 7, 32 and 33 cells and two tolerance regimes. It compares rates, forcing, Jacobian, the flat factorization, the solve, and six steps with step counts and a mid-run state copy, at 1e-10 per value plus a relative infinity-norm check.

**Scores**
- **Readability 6.** The modules are small and their docstrings state each pass's purpose; rate formulas use named keys. The solver path, however, is produced by string surgery:
  - `fixed_kernels` regex-renames functions and inserts `const size_t cells = 4;` after every `{`.
  - `transform_solver` replaces specific lines of the template (`tools/tiled.py:11-91`).
  - `split_kernels` re-parses its own output by indentation.
  
  To see what actually runs, a reader has to apply these rewrites in their head.
- **Maintainability 8.** The tests are strong, and Make tracks every generator module and the template (`Makefile:142`) along with option and compiler stamps. But the text replacements would fail silently if the template changed. The kernel test also checks the flat Mozart factor (`check.c:67`), not the packed Doolittle LU that the timed solver uses; that path is tested only through trajectories. There are no bad-input tests for the generator.
- **Extensibility 7.** Configuration edits regenerate, and the README documents generating another solver and the order of calls to use it. Adding a mechanism to the suite needs a rule, an `SPECIAL_OBJECTS` entry and two array entries in `check.c`. A new rate type touches `rate_initializer`, `formula` and `src/rates.c`. New layout specializations would have to be written as more text transforms.

**Exercises**
- (a) **easy:** edit `configs/ts1.json` and copy the per-reaction `"concentration units": "mol m-3"` key (leaving it out silently applies the molecule-cm⁻³ conversion); `make bench test-specialized`.
- (b) **easy:** edit the value; same commands.
- (c) **moderate**
  - `.conda/bin/python tools/specialize.py ../mechanisms/chapman/config.json generated/chapman2_specialized.c --symbol chapman2_specialized`.
  - Add a rule like the Robertson one and add `chapman2` to `SPECIAL_OBJECTS`.
  - Add the extern and entries in `solvers[]` and `names[]` of `check.c`; `make test-specialized`. The shared caveats apply.

**Strengths**
1. Modular generator with a symbolic layer shared with Tier 1.
2. Strong tests on seven mechanisms, with both per-value and norm checks.
3. Complete Make dependencies and documented options and experiments.

**Weaknesses**
1. Tiling implemented by fragile text rewriting.
2. The timed LU kernel is never compared directly, and every kernel exists twice in the generated file.
3. The per-reaction units key is a trap when adding a reaction; unknown reaction keys are ignored.

---

## Entry D

**Summary.** `tools/specialize.py` (406 lines) reuses the validation in `tools/configure.py` and runs five passes, described in `DESIGN.md:195-216`:
1. stable rate order and custom indices;
2. MICM Markowitz order and Mozart fill (`Plan`);
3. ordered forcing and derivative operations;
4. emission of rate formulas, pivots and substitutions as SIMD loops in bounded helpers;
5. instantiation of the shared integrator source `src/rosenbrock_step.inc`.

Kernels use a fixed stride (default 32) and fall back to the generic library for any other group size. An exact-input cache skips rate evaluation when the conditions and parameters are byte-identical. `tests/test_specialized.c` covers TS1, three tutorials, an extra-rate fixture and a mutated fixture (adds species and reactions, changes a coefficient) at 1, 7, 16, 32 and 33 cells and two tolerance regimes. It checks rates and cache invalidation, LU and triangular solves on a synthetic matrix, singular-pivot detection, and twelve steps with statistics, at a pure 1e-10 relative bound. The generator writes depfiles listing every included configuration file.

**Scores**
- **Readability 7.** There is one function per pass, and comments explain the subtle numerics: the unsigned Markowitz cost, the first-occurrence rule, and why the singular-pivot check is a SIMD reduction (`tools/specialize.py:44-46,216,262-273`). Two things are harder to follow: the integrator is instantiated by redefining library names (`#define chem_forcing forcing`, `:320-335`), and several emitters build long one-line C strings.
- **Maintainability 8.** Depfiles plus settings and compiler stamps make regeneration exact (`Makefile:117-165`), and the reader has 49 rejection tests. However, at the default width 32 every other cell count runs generic kernels (`specialize.py:186,226,291`), so only the 32-cell runs actually test generated code; forcing and Jacobian are never compared directly. `DESIGN.md:191-193` still says generated C "supplies mechanism data, not specialized arithmetic".
- **Extensibility 8.** The README gives exact commands to generate and compile another mechanism, and adding it to the suite means one list, one rule and `DECLARE`/`CASE`. DESIGN states every place a new rate type must touch, including `rate_statements`. The tuning options pass through `GEN_OPTIONS`.

**Exercises**
- (a) **easy:** edit `mechanisms/ts1.json`; `make bench test-specialized`.
- (b) **easy:** edit the value; same commands.
- (c) **moderate**
  - `.conda/bin/python tools/specialize.py ../mechanisms/chapman/config.json build/specialized/chapman.c --symbol chapman --width 32 --depfile build/specialized/chapman.d` and `make build/specialized/chapman.o` (README).
  - Add `chapman` to `SPECIALIZED_NAMES` with a generation rule, and add `DECLARE(chapman)`/`CASE(chapman)` to the test; `make test-specialized`.
  - Care: use `GEN_WIDTH=0` to test the generated kernels at sizes other than 32. The shared caveats apply.

**Strengths**
1. Clear per-pass generator with the subtle numerics commented.
2. One Rosenbrock source shared by the generic and generated solvers, and exact depfile-based regeneration.
3. Strong reader validation, a mutation fixture, and documented commands for new mechanisms.

**Weaknesses**
1. The default test exercises generated kernels only at 32 cells.
2. The macro-override instantiation of the integrator is obscure.
3. A stale sentence in DESIGN.md.

---

## Entry E

**Summary.** `tools/specialize.py` (350 lines) normalizes through `tools/configure.py` (reorder and unit conversion) and builds the filled CSR pattern. It emits direct reaction and Jacobian code, rate formulas and an SIMD error-component kernel. For the LU it emits either unrolled code or, by default, a compact table interpreter (`tools/compact_linear.py`); `tools/compact_chemistry.py` optionally shares kernels between reactions of the same shape. Kernels have a fixed width of 32 and fall back to generic kernels at other sizes, and the shared `src/rosenbrock_impl.h` is instantiated through hooks. `tests/test_specialized.py` (using ctypes) checks nine variants at 1, 7, 32 and 33 cells. It checks every sparse-plan coordinate, rate constants at 1e-13, six steps at both tolerance regimes with all statistics, and the other Rosenbrock tables.

**Scores**
- **Readability 6.** Functions are well separated, and DESIGN's two option tables match argparse. Rate formulas use positional slots whose meaning is defined only in `configure.py` (`specialize.py:182-233`), and the default LU is a generated interpreter over static tables, which is less direct to read than straight-line code.
- **Maintainability 7.** The plan-coordinate and statistics checks are strict, and Make tracks every module plus option and compiler stamps (`Makefile:100-157`). But forcing, Jacobian and LU are never compared kernel by kernel, and generated kernels run only at 32 cells (7 in two variants). The non-TS1 harness also starts every species at 0 except one named `A` (`tests/test_specialized.py:65`). Only the include-cycle case has a bad-input test.
- **Extensibility 6.** Generating another mechanism is a one-line command, but the test looks only in `configs/` and starts other mechanisms from all-zero states. Carbon Bond 5 is rejected outright ("duplicate custom rate label NO2", `configure.py:152-157`), and the rate-type guide omits `rate_formula`.

**Exercises**
- (a) **easy:** edit `configs/ts1.json`; `make bench test-specialized`.
- (b) **easy:** edit the value; same commands.
- (c) **moderate**
  - `.conda/bin/python tools/specialize.py ../mechanisms/chapman/config.json build/chapman_specialized.c --prefix chapman --fixed-width 32`.
  - Copy the configuration under `configs/`, add a rule and an entry in `SPECIALIZED_CONFIGS`, and add an entry to `CONFIGS` in the test.
  - Fix the harness's initial values, or the trajectories stay at zero.
  - Carbon Bond 5 cannot be run as given.

**Strengths**
1. Well-factored emitter with documented alternatives.
2. Strict comparison of the sparse plan and of statistics across nine variants.
3. Explicit errors for duplicate labels and missing surface properties.

**Weaknesses**
1. Positional rate parameters.
2. Narrow kernel coverage, and trivially zero states for non-TS1 mechanisms in the harness.
3. Rejects Carbon Bond 5.

---

## Entry F

**Summary.** `tools/generate_c.py` (651 lines) reads either F's own normalized JSON (used for TS1) or a CAMP configuration. It uses its own Markowitz tie rule (2,360 TS1 entries) and no-pivot fill. With `--specialized` it emits per-cell scalar functions for rates, custom-rate updates, combined forcing and Jacobian, an LU using `fma` with pivot checks, and triangular solves. These go into an optional `chem_specialized_kernels` table on the mechanism, which the generic solver calls when the sparse pattern matches; the integrator, error norm and controller stay in the library. `--without-*-kernel` flags reproduce the ablations. The test (`tests/test_specialized.c`) runs TS1, Robertson and the two tutorials with 4 cells, three 30-s steps and default tolerances, and compares final states only, at 1e-10.

**Scores**
- **Readability 5.** One file mixes reading, ordering, symbolic work and emission; the docstrings help, and the emitted C has descriptive temporaries. Rate parameters are positional, unit fixes are ad hoc, and DESIGN's generator section is short.
- **Maintainability 4.**
  - Tests: one cell count, three steps, final state only, no kernel or statistics checks.
  - Build: no flag stamp — the README asks for `make clean` between compilers (`README.md:21-29`); the tutorial rules omit their included files (`Makefile:160-167`); `make env` is not idempotent.
  - The CAMP reader converts only Troe (with a fixed exponent) and ternary k_inf, never Arrhenius A, and expects `scaling_factor`, not MICM's `scaling factor` (`generate_c.py:25-29,78-94`). The tests cannot notice, because the generic path uses the same reader.
- **Extensibility 4.** A new rate type stays local to `KINDS` and `emit_rate_body`. But another CAMP mechanism is misread unless the reader is fixed, its harness uses zero initial concentrations, and each binary holds a single mechanism.

**Exercises**
- (a) **easy:** append a reaction to `mechanisms/ts1.json` in F's undocumented schema (`parameters` in SI units, products as `{name, yield}`); `make bench test-specialized`.
- (b) **easy:** edit `parameters`; same commands.
- (c) **hard**
  - `.conda/bin/python tools/generate_c.py ../mechanisms/chapman/config.json build/specialized_chapman.c --specialized`.
  - Add generation, object and test-binary rules, listing the included files, and add the binary to `test-specialized`.
  - Add `initial concentration` to the species, or the comparison is vacuous.
  - Fix the reader before trusting the result: bimolecular Arrhenius A stays in cm³ molecule⁻¹ s⁻¹ units, and Carbon Bond 5's ×9 and ×0.0036 scaling factors are dropped.

**Strengths**
1. Simple optional kernel table, chosen at run time, with measured ablations.
2. Readable scalar generated code from a single generator.
3. Tests cover all four configured mechanisms.

**Weaknesses**
1. The CAMP reader misreads general mechanisms.
2. Weak tests.
3. Fragile build: no stamps, missing include dependencies.

---

## Entry G

**Summary.** `tools/specialize.py` (395 lines) uses the reader and symbolic layout in `tools/configure.py`, which also emits the Tier 1 loader and Make depfiles. It emits fixed rate formulas, fused chemistry kernels that assign on first write, and an unrolled Mozart LU with substitutions, in four loop styles. It then fills either `tools/rosenbrock3.c.in` or, by default, `tools/rosenbrock3_blocked.c.in`. The blocked template runs all three stages in four-cell tiles inside each 32-cell adaptive group and keeps the group's serial error sum. Coefficients and the controller come from the library. `tests/specialized.py` checks TS1, three tutorials, an all-family validation mechanism, an edge fixture, and a TS1 variant from `tools/specialization_fixture.py` that adds a reaction and changes a coefficient. At 1, 2, 4, 7, 32 and 33 cells it compares rates, forcing, Jacobian, LU, solves and three intervals with T, P, air density and parameters changing between them, at 1e-10 relative, with expected zeros matched exactly.

**Scores**
- **Readability 6.** The docstrings and named `rate_body` help, and the templates are plain C. But four linear-loop modes, fixed and dynamic strides, and optional tiles make `emit_linear` and `specialize` dense: the latter takes 12 arguments and chains `.replace` calls (`specialize.py:296-362`). DESIGN.md is written in long paragraphs.
- **Maintainability 9.** Its regeneration handling is the most complete: depfiles for recursive includes (`configure.py:108-119`, `Makefile:187`), per-compiler build directories, and flag and option stamps. Every kernel is tested at every cell count, including an automatic add-reaction/change-coefficient variant. Minor issues: DESIGN calls its order the "reference" Markowitz, yet it counts the diagonal (2,359 entries against the 2,381 reported for MICM's order), and it merges parameters that share a label (`specialize.py:18-32`).
- **Extensibility 8.** Editing TS1 is already exercised automatically, the README documents generating a new solver, and DESIGN names every place a new rate type must touch. Adding a mechanism to the suite takes about five Makefile edits plus the test's names tuple.

**Exercises**
- (a) **easy:** edit `configs/ts1.json` (`reactant sequence`/`product sequence` lists); `make bench test-specialized`. `ts1_variant` already automates a version of this.
- (b) **easy:** edit the value; same commands.
- (c) **moderate**
  - `.conda/bin/python tools/specialize.py ../mechanisms/chapman/config.json generated/chapman_specialized.c --name chapman_specialized`.
  - In the Makefile: `SPEC_NAMES`, `CONFIG_chapman`, a `generated/chapman.c` rule, and `$(BUILD)/chapman.o` added to `libspecialized`.
  - Add `'chapman'` to `main()` in `tests/specialized.py`. The shared caveats apply.

**Strengths**
1. A built-in test of editing a reaction and a coefficient, with every kernel compared at all cell counts.
2. Exact, compiler-aware regeneration.
3. Every option and extension point is documented.

**Weaknesses**
1. Many options make the generator hard to read.
2. Dense documentation, and a slightly inaccurate "reference ordering" claim.
3. Same-label parameters are merged, and adding a mechanism takes many Makefile edits.

---

## Entry H

**Summary.** `tools/generate.py` (843 lines) reads H's normalized JSON or a CAMP file and picks its own minimum-fill ordering. Into one file it writes both the generic descriptor and a specialized solver: rate loops, forcing and Jacobian SIMD loops, an LU split into one function per pivot, table-driven triangular solves, and a full three-stage integrator written as C text inside Python (`_emit_integrator`). The integrator transposes the cell-major state, allocates its workspace on every call, and carries its own copy of the coefficients. The test (`tests/test_specialized.c`) runs TS1 and two tutorials on three identical cells for one 1-ms interval and compares the final state with the generic path.

**Scores**
- **Readability 4.** The rate code is readable, but normalization is done twice (`solver_mechanism`, `emit`), about 150 lines of solver sit in Python strings (`generate.py:541-690`), and dead code remains (lines 790-791).
- **Maintainability 3.** The test is very weak (`tests/test_specialized.c:23-52`). The Makefile has no flag stamp and a non-idempotent `env`, and the generic benchmark compiles the specialized file. The CAMP reader does not convert Arrhenius A (`generate.py:131-136`), follows includes only one level, and forces initial concentrations to 0 (line 71).
- **Extensibility 3.** DESIGN names `_rate_code` for new rate types. But each binary holds one mechanism (fixed `chem_generated_*` symbols), no command is documented for other mechanisms, and changing the integrator means editing C inside Python strings.

**Exercises**
- (a) **easy:** edit `config/ts1.json` (types such as `ArrheniusRateConstant`, parameters in SI units); `make bench test-specialized`. The test barely exercises the change.
- (b) **easy:** edit the value; same commands.
- (c) **hard:** `.conda/bin/python tools/generate.py ../mechanisms/chapman/config.json build/generated/chapman.c`; add object and test rules. The harness's states are all zero because CAMP input cannot set initial concentrations, the Arrhenius conversion is missing, and a third body is recognized only through `parameterized: air_density`.

**Strengths**
1. Readable rate emitters with named parameters.
2. Species-major generated loops.
3. Candid measurements in SPECIALIZATION.md.

**Weaknesses**
1. Integrator embedded in Python strings, and duplicated normalization.
2. Near-vacuous test.
3. Reader gaps, and one mechanism per binary.

---

## Entry I

**Summary.** `tools/generate_c.py` (380 lines, most of it one 280-line `main()`) reads I's normalized JSON; CAMP files first go through `tools/read_config.py`. It emits the mechanism tables and, under `#ifdef CHEM_SPECIALIZED`, four per-cell scalar functions: forcing, a negative Jacobian in a dense n×n layout, a fixed-index LU over the Markowitz-ordered species, and triangular solves. Rate constants, the integrator and the controller stay generic. The test runs both TS1 benchmarks on the 100-cell reference input (100 × 30 s) and compares them with `compare_final.py` at 1e-10. For two tutorials read from `../micm`, it compares rates, forcing and Jacobian for a single cell at 1e-13, plus one 1e-6-s solve.

**Scores**
- **Readability 4.** The code reads linearly but is not decomposed: everything is emitted inline in `main()`, with positional slots (`values[7]` holds the custom index, line 201). DESIGN's Tier 2 section is short.
- **Maintainability 4.** The TS1 end-to-end check is meaningful. The tutorial checks are tiny, and they depend on configuration files outside the repository. The Makefile recompiles everything in one command per binary, has no stamps, hard-codes the libomp path, and has a non-idempotent `env`. The reader skips the Arrhenius conversion and the `scaling factor` key (`read_config.py:24-79`).
- **Extensibility 3.** Only TS1 is wired into Make, generating another mechanism takes two steps plus hand compilation, and rate constants are not specialized. New specializations would have to go inside the monolithic `main()`.

**Exercises**
- (a) **easy:** edit `mechanisms/ts1.json` (`{name, coefficient}` terms); `make bench test-specialized`.
- (b) **easy:** edit the value; same commands.
- (c) **hard:** `read_config.py ../mechanisms/chapman/config.json build/chapman.json`, then `generate_c.py build/chapman.json build/chapman.c`, then compile with `-DCHEM_SPECIALIZED`, the library sources and `tests/specialized_compare.c` (the test script only looks in `../micm`). The Arrhenius A values come out wrong unless the reader is fixed.

**Strengths**
1. Meaningful TS1 end-to-end comparison.
2. Clean `CHEM_SPECIALIZED` switch and a small table of specialized operations.
3. Small generated source (about 1.1 MB).

**Weaknesses**
1. Monolithic generator.
2. Reader gaps, and tests that depend on `../micm`.
3. Minimal build discipline, and only TS1 is supported in practice.

---

## Comparison

| Rank | Entry | Readability | Maintainability | Extensibility | Sum |
|---|---|---:|---:|---:|---:|
| 1 | B | 7 | 9 | 8 | 24 |
| 2 | G | 6 | 9 | 8 | 23 |
| 3 | D | 7 | 8 | 8 | 23 |
| 4 | A | 6 | 8 | 7 | 21 |
| 5 | C | 6 | 8 | 7 | 21 |
| 6 | E | 6 | 7 | 6 | 19 |
| 7 | F | 5 | 4 | 4 | 13 |
| 8 | I | 4 | 4 | 3 | 11 |
| 9 | H | 4 | 3 | 3 | 10 |

**Why this order**
- **B is first.** It combines a readable generator (named formulas, explicit schedules) with the most complete verification and a clean descriptor API.
- **G and D tie at 23.** G has the stronger verification and rebuild handling; D has the cleaner generator and more precise extension notes. I put G ahead because its suite checks every generated kernel at every cell count and automatically tests adding a reaction and changing a coefficient, while D's default suite exercises generated kernels only at 32 cells.
- **A and C tie at 21.** I put A ahead because its emission is explicit and the LU the timed path actually uses is compared directly. C builds its timed path by text rewriting and checks it only through trajectories.
- **E (19)** is well structured but held back by positional rate parameters, narrow kernel coverage, and refusing Carbon Bond 5.
- **A gap separates F (13), I (11) and H (10).** All three have single-file or monolithic generators and weak tests. Their CAMP readers also skip MICM's molecule-cm⁻³ conversion for Arrhenius rates, so a new mechanism is misread identically by the generic and generated paths, which their tests cannot catch. F is ahead for its clean kernel table and four-mechanism test; I is ahead of H for its meaningful TS1 end-to-end test.

**Shared caveats** (they apply to every exercise (c) above)
- No entry reads MICM's `tracer type: THIRD_BODY`, so Chapman's M becomes an integrated species everywhere.
- No entry reads `Ea`. It is 0 throughout Carbon Bond 5, so nothing changes here, but a nonzero value would be silently ignored.
- No entry prefixes photolysis and emission labels the way MICM does (`PHOTO.`, `EMIS.`), so Carbon Bond 5's `NO2` and `ALD2` labels collide. A, F, G, H and I silently merge them into one parameter; E refuses the mechanism; B, C and D keep two slots under the same name.
- No generator annotates the generated C with species or reaction names (A's JSON sidecar comes closest), and none documents an end-to-end recipe for adding a mechanism to `make test-specialized`.

```json
{"scores": {"A": [6, 8, 7], "B": [7, 9, 8], "C": [6, 8, 7], "D": [7, 8, 8], "E": [6, 7, 6],
            "F": [5, 4, 4], "G": [6, 9, 8], "H": [4, 3, 3], "I": [4, 4, 3]},
 "exercises": {"A": ["easy", "easy", "moderate"], "B": ["easy", "easy", "moderate"], "C": ["easy", "easy", "moderate"], "D": ["easy", "easy", "moderate"], "E": ["easy", "easy", "moderate"], "F": ["easy", "easy", "hard"], "G": ["easy", "easy", "moderate"], "H": ["easy", "easy", "hard"], "I": ["easy", "easy", "hard"]},
 "rank": ["B", "G", "D", "A", "C", "E", "F", "I", "H"]}
```
