# Blind review of the nine Tier 2 generators

**How I worked.** I read only files in this folder and did not build or run anything. Every claim about generated code comes from reading the generators, templates and tests. I could not check claims about MICM's own reader against MICM's source, which is not in this folder.

**Difficulty scale for the exercises:**
- **easy:** edit the configuration, run `make`, and the generated solver plus `make test-specialized` work with no other code edits.
- **moderate:** build files or test files also need editing, or the change needs real care.
- **hard:** the generator or library must change.

For every entry, exercises (a) and (b) change the generic TS1 binary as well, because both binaries read the same configuration. Any test that compares TS1 with MICM's results will then fail, and that failure is expected.

**Reader gaps that matter for exercise (c):**
- **No entry handles everything in Chapman and Carbon Bond 5.** None of the nine readers recognizes `"tracer type": "THIRD_BODY"`, so Chapman's M becomes an ordinary variable species, and none reads `"Ea"`.
- **Labels collide.** None namespaces user-defined labels the way MICM does with `PHOTO.`/`EMIS.`. Carbon Bond 5's photolysis and emission parameters named `NO2` and `ALD2` therefore collide. Most entries merge them silently; E rejects the file with an error.
- **Unit conversion splits the field.** A–E and G convert pre-exponential factors by reaction order, which fits all three tutorial fixtures. F, H and I hard-code the factors that fit the tutorial's particular Troe and ternary reactions and leave Arrhenius A unconverted.
- **No generated C points back to the chemistry.** None of the nine annotates its output with reaction equations or species names.

---

## Entry I

**Summary.**
- **Generator and input:** `tools/generate_c.py` (380 lines, almost all in one `main()`) reads a project-specific normalized JSON. TS1 is the committed `mechanisms/ts1.json`, made once by `tools/import_ts1.py` using regular expressions over MICM's C++ header; CAMP configurations go through `tools/read_config.py`.
- **What it specializes:** it emits the static tables used by both the generic and specialized binaries. Under `#ifdef CHEM_SPECIALIZED` it adds straight-line forcing, negative Jacobian, a natural-order LU on the Markowitz-reordered species, and triangular solves.
- **Storage:** all of these index a dense n×n row-major array with constant offsets (`src/solver.c:109–110` allocates n² doubles).
- **What stays generic:** rate constants, the Rosenbrock controller and the per-cell loop.
- **Tests:** `tests/test_specialized.py` compares the full generic and specialized TS1 runs (100 cells, 100 steps) at 1e-10 with `compare_final.py`. For the two MICM tutorial configurations, `tests/specialized_compare.c` compares forcing and Jacobian at 1e-13 and one 1e-6 s solve.

**Scores**
- **Readability: 4.** `main()` mixes parameter packing, pattern and fill construction, and C emission in f-strings, with two comments (`generate_c.py:93–376`); only Markowitz ordering and string quoting are factored out. DESIGN.md is mostly the Tier 1 API, with a short, accurate Tier 2 paragraph. SPECIALIZATION.md has two measured rows. The generated code (`rate_17`, `f[42]`) gives no way back to the chemistry.
- **Maintainability: 5.** The end-to-end TS1 test is real, and `Makefile:37–38` regenerates `generated/ts1.c` from the JSON and the generator. But the reader silently mishandles MICM inputs:
  - `"scaling factor"` passes through `read_config.py`, while `generate_c.py:205` reads `scaling_factor`, so CB5's 0.0036 and 9 become 1.
  - Arrhenius A is never converted, and Troe uses fixed exponents (`read_config.py:75–79`).
  - An unknown type raises a bare `KeyError` at `generate_c.py:99`, before the friendly error at `:150`.
  - There are no compiler-flag stamps.
- **Extensibility: 4.** Editing TS1 data is easy. There is no documented route for other mechanisms: the commands exist only in `tests/test_specialized.py:47–55`, with tutorial paths hard-coded. A new specialization or rate kind means editing the monolithic `main()` plus the C enum and evaluator (`DESIGN.md:152–154`).

**Exercises**
- **(a) Add a reaction — easy**
  - Append `{"type": "arrhenius", "reactants": [{"name": …, "coefficient": 1.0}], "products": […], "parameters": {"A": …, "C": …}}` to `mechanisms/ts1.json`. A is in mol m⁻³ units; nothing converts it on this path.
  - Run `make bench`, then `make test-specialized`. Ordering and fill update automatically.
  - Care: a misspelled species gives a raw `KeyError`.
- **(b) Change a rate parameter — easy**
  - Edit `"A"`/`"C"` in that reaction's `parameters`, then run `make bench test-specialized`.
- **(c) Another mechanism — moderate**
  - Run `.conda/bin/python tools/read_config.py ../mechanisms/chapman/config.json /tmp/ch.json`, then `tools/generate_c.py /tmp/ch.json /tmp/ch.c`.
  - Compile `/tmp/ch.c` with `tests/specialized_compare.c` and the library sources using `-DCHEM_SPECIALIZED`, copying the command at `test_specialized.py:51–54` (or edit that script's hard-coded loop), and run it.
  - Expect for Chapman: M treated as ordinary, because the reader looks for a `"third body"` key (`read_config.py:36`).
  - Expect for CB5: scaling factors dropped; bimolecular A values about 6e17 off compared with an order-based reader; unimolecular Troe given bimolecular factors; NO2/ALD2 labels merged (`generate_c.py:113–124`). The generic-vs-generated check still passes and hides all of this.

**Strengths**
1. The TS1 test compares complete 100-cell, 100-step runs of the two binaries.
2. The specialization hook is small and clean (`chem_specialized_ops`, four function pointers).
3. The generator needs only the standard library, and make regenerates its output.

**Weaknesses**
1. A monolithic, almost uncommented generator, with thin generator documentation.
2. The CAMP reader silently diverges from MICM's format (units, scaling factor, label namespaces, third body).
3. Rate evaluation is unspecialized, and the matrix is dense n×n.

## Entry H

**Summary.**
- **Generator and input:** `tools/generate.py` (843 lines) reads H's normalized JSON (`config/ts1.json`, from `tools/translate_ts1.py`) or a CAMP configuration (`from_camp`). It emits the generic descriptor plus a complete specialized solver.
- **What it specializes:** per-reaction rate formulas; cell-vectorized forcing and Jacobian on a species-major group layout; one generated function per LU pivot; a table-driven triangular solve; an error norm.
- **Integrator:** the whole three-stage Rosenbrock integrator is written out from Python string literals (`_emit_integrator`, `generate.py:541–690`).
- **Ordering:** a minimum-actual-fill heuristic in `solver_mechanism`.
- **Tests:** `tests/test_specialized.c` compares the generated integrator with the generic one at 1e-10 for TS1 and the two tutorial configurations, but on three identical cells over a single 1 ms interval.

**Scores**
- **Readability: 5.** The functions are sensibly split (`from_camp`, `solver_mechanism`, `_rate_code`, `specialized_source`, `emit`), and DESIGN.md's Tier 2 section explains ordering and layout. But about 150 lines of integrator live as Python strings, `solver_mechanism` returns a dict of a dozen parallel structures, and comments are sparse. SPECIALIZATION.md is candid but has no isolated measurements.
- **Maintainability: 4.** `Makefile:49–51` regenerates TS1 when the configuration, generator or header changes. Against that:
  - Identical cells cannot catch cell-stride errors, and 1 ms barely exercises the controller.
  - CAMP species start at 0 (`generate.py:71`), so the no-user tutorial's state stays identically zero during the test.
  - The string-emitted integrator duplicates the generic controller in `src/chem.c` and must be kept in step by hand.
  - Objects carry no compiler stamp, so README's `make cpu CC=clang` followed by `make cpu CC=gcc-15` silently reuses clang objects.
- **Extensibility: 5.** TS1 edits are easy, and `DESIGN.md:171–174` lists the four places a new rate type touches. A new mechanism needs new Makefile rules and a test harness that seeds nonzero concentrations. Unit conversion is a tutorial special case (`generate.py:127–136`).

**Exercises**
- **(a) Add a reaction — easy**
  - Add `{"type": "ArrheniusRateConstant", "reactants": [names], "products": [{"species", "yield"}], "label": null, "parameters": {A…E}}` (SI units) to `config/ts1.json`.
  - Run `make bench test-specialized`.
  - Breaks: `tests/test_rates.c:217–218,242` pins 547 reactions, 168 parameters and 142 user rates, so `make test-cpu` fails until updated.
- **(b) Change a rate parameter — easy**
  - Edit `parameters` in `config/ts1.json`, then run `make bench`.
- **(c) Another mechanism — moderate**
  - Run `tools/generate.py ../mechanisms/chapman/config.json build/generated/chapman.c`; it accepts `camp-files`.
  - Copy the `NO_USER_GEN` / `specialized_no_user` / `test_specialized_no_user` rules (`Makefile:80–82,268–283`) and add the new binary to `test-specialized`.
  - Edit `tests/test_specialized.c` to set nonzero concentrations and photolysis rates; otherwise the comparison is vacuous.
  - M is not a third body (the reader expects `"parameterized": "air_density"`). CB5's labels merge, and its Arrhenius A values go unconverted.

**Strengths**
1. One generator specializes everything, including rates and the integrator, on a cell-vectorized layout.
2. Clear errors for unknown species and types and for non-integer quantities; `"scaling factor"` is read correctly.
3. Make dependencies cover the configuration, generator and header.

**Weaknesses**
1. The integrator is duplicated as Python strings.
2. Tests are close to vacuous for configuration mechanisms and weak for TS1.
3. No compiler stamp (contradicting README's build sequence), and unit conversion only fits the tutorial.

## Entry G

**Summary.**
- **Two generators, one reader:** `tools/configure.py` is the Tier 1 path. It emits a loader that builds the model through the public API, plus fixed-index forcing and Jacobian kernels. `tools/specialize.py` is the Tier 2 path.
- **What the Tier 2 path specializes:** rate formulas; unrolled Mozart LU and substitutions with selectable loop backends; a complete three-stage integrator from `tools/rosenbrock3_blocked.c.in`, which packs four-cell tiles inside each 32-cell adaptive group.
- **Shared constants:** the templates take their coefficients from the library's `chem_rosenbrock_parameters(3)` and `chem_default_controller`.
- **Reader:** handles recursive includes, order-based units with a `concentration units` extension, and writes make depfiles.
- **Tests:** `tests/specialized.py` (via ctypes) compares rates, forcing, Jacobian, LU, solves and three integrated intervals at `1e-10·|generic|`, and expected zeros must match exactly. It covers TS1, three tutorials, an all-rate validation mechanism, an edge fixture, and a TS1 variant that make builds by changing an A and adding a reaction (`tools/specialization_fixture.py`), each at 1, 2, 4, 7, 32 and 33 cells.

**Scores**
- **Readability: 6.** The pieces are named and have docstrings: `custom_indices`, `rate_body`, `linear_plan` (operations as tuples), `emit_linear`. DESIGN.md tabulates every parameter. But the code is very dense (for example `configure.py:156–216`), and four LU loop backends × fixed/dynamic stride × two templates hide the default path. The prose documentation is hard going.
- **Maintainability: 8.** The best build hygiene of the nine: compiler-flag and generation-option stamps, plus depfiles for transitive includes (`Makefile:81–147,187`). The tests are broad and include the build-made TS1 mutation, and coefficients have a single source. Minuses: the two templates repeat the control flow of `src/rosenbrock.c`; the tests never compare solver statistics; `tests/tutorials.c:276` pins TS1 at 210 species and 547 reactions.
- **Extensibility: 7.** The build already performs exercises (a)/(b) on its variant. Both generators' command lines are documented, and so are the rate-type steps including `rate_body` (`DESIGN.md:64`). Adding a mechanism to the tests takes about five coordinated Makefile edits plus `tests/specialized.py:129`.

**Exercises**
- **(a) Add a reaction — easy**
  - Add a CAMP-form reaction to the TS1 MECHANISM in `configs/ts1.json`. The file-level `concentration units: mol m-3` applies, so A is in SI units.
  - Run `make bench test-specialized`.
  - The "TS1 dimensions" check (`tests/tutorials.c:276`) fails in `make test-cpu` and needs updating.
- **(b) Change a rate parameter — easy**
  - Edit A/B/C, then run `make bench test-specialized`.
- **(c) Another mechanism — moderate**
  - Makefile: add `example_chapman` to `SPEC_NAMES`, set `CONFIG_example_chapman := ../mechanisms/chapman/config.json`, add a `generated/example_chapman.c` rule (calling `configure.py`), and add `$(BUILD)/example_chapman.o` to `GENERATED_OBJECTS`.
  - Add the name to the tuple at `tests/specialized.py:129`, then run `make test-specialized`.
  - The harness seeds 1e-8 concentrations, so the comparison is meaningful.
  - M is ordinary. CB5's photolysis and emission NO2 share one parameter, and Ea is ignored; units, EMISSION and scaling factor are handled.

**Strengths**
1. The most thorough generator tests, including a build-generated TS1 mutation.
2. Make regenerates on any change to the configuration, its includes, the generator, or the options.
3. The integrator templates take their coefficients from the library; the options are documented with measurements.

**Weaknesses**
1. Dense code, and many options obscure what the default emits.
2. The integrator templates duplicate the generic control flow.
3. Labels are not namespaced; THIRD_BODY and Ea are not recognized; the TS1 size is pinned in a Tier 1 test.

## Entry F

**Summary.**
- **Generator and input:** `tools/generate_c.py` (651 lines) reads F's normalized JSON (`mechanisms/ts1.json`, from `tools/translate_ts1.py`) or a CAMP configuration, and always emits the immutable tables for the generic library.
- **What `--specialized` adds:** per-reaction rate formulas, a custom-rate update, a combined forcing/negative-Jacobian kernel, a fixed LU with explicit `fma` updates, and fixed triangular solves.
- **How it plugs in:** these go into a `chem_specialized_kernels` table, which the generic solver calls when the sparsity pattern matches. The kernels work per cell, and the controller stays in the library.
- **Ablations:** the `--without-*-kernel` flags build the variants reported in SPECIALIZATION.md.
- **Tests:** `tests/test_specialized.c` integrates 4 cells with different T and P over three 30 s intervals, scaling custom parameters each interval, and compares generated with generic at 1e-10 for TS1, Robertson and both tutorials.

**Scores**
- **Readability: 7.** A single, well-factored file with docstrings (`read_config`, `markowitz_order`, `sparse_program`, `emit_rate_body`, `emit_specialized_kernels`, `compile_model`). Its emitted C uses named temporaries (`pivot_k`, `multiplier_k_row`, `derivative_i_a`) and is the easiest of the nine to trace by hand. The Tier 2 documentation and the ablation table are concise and accurate. The reader embeds tutorial-specific unit factors under a comment that overstates their generality (`generate_c.py:88–94`).
- **Maintainability: 5.** The optional kernel table keeps the generator and solver cleanly separated, and the ablation targets are reproducible (`Makefile:117–144`). But:
  - CAMP species start at 0 and the harness perturbs only nonzero values (`test_specialized.c:74–76`), so the Robertson and no-user tests are vacuous.
  - `"scaling factor"` is silently dropped (`KINDS` reads `scaling_factor`).
  - Troe factors are fixed, and Arrhenius A is not converted.
  - There is no compiler stamp (README says to clean between compilers).
- **Extensibility: 6.** Adding a kernel or rate formula stays inside one Python file plus `kinetics.c`. But `make mechanism CONFIG=…` builds only the generic library, so comparing a new mechanism needs three new Makefile rules and harness changes.

**Exercises**
- **(a) Add a reaction — easy**
  - Add `{"type": "ARRHENIUS", "reactants": [names], "products": [{"name", "yield"}], "parameters": {…}}` (SI units) to `mechanisms/ts1.json`.
  - Run `make bench test-specialized`.
- **(b) Change a rate parameter — easy**
  - Edit the parameter, then run `make bench test-specialized`.
- **(c) Another mechanism — moderate**
  - Add a rule `$(BUILD)/specialized_chapman.c: ../mechanisms/chapman/config.json tools/generate_c.py` (calling `--specialized`), a `$(BIN)/test_specialized_chapman` link rule, and add it to `test-specialized`.
  - Edit the harness to seed nonzero concentrations; Chapman supplies none, so the test would otherwise compare zeros.
  - M is ordinary. CB5 loses its scaling factors and Arrhenius conversion, gets bimolecular factors on its unimolecular Troe reactions, and merges NO2/ALD2.

**Strengths**
1. A compact, readable generator whose generated code is traceable.
2. A clean optional-kernel interface with fallback to the generic path.
3. An honest per-kernel ablation table.

**Weaknesses**
1. Tests of non-TS1 mechanisms are largely vacuous.
2. The reader diverges from MICM's format (units, scaling factor).
3. Every new mechanism needs hand-written Makefile rules, and there is no compiler stamp.

## Entry E

**Summary.**
- **Reader:** `tools/configure.py` normalizes CAMP configurations (includes with cycle detection, order-based units with a `rate units` extension) and makes duplicate custom labels an error.
- **What `tools/specialize.py` emits:** fixed-index forcing and Jacobian, an alpha shift, the LU, triangular solves, rate formulas and a SIMD error norm. The default LU is a compact operation table interpreted in short SIMD loops (`compact_linear.py`); `--unrolled-linear` gives straight-line code.
- **Integrator:** it instantiates the library's own algorithm, `src/rosenbrock_impl.h`, through macros, so the generic and specialized solvers share one integrator.
- **Tests:** `tests/test_specialized.py` runs 9 configurations at 1, 7, 32 and 33 cells: TS1, three tutorials, parameterized-species, all-feature and inert fixtures, and two builds with alternative generator options. It checks:
  - every Jacobian coordinate against the C plan;
  - rates under probed parameters, at 1e-13;
  - six steps at default and tight tolerances with T, P, density and parameters changing, at 1e-10 with every statistic equal;
  - the other Rosenbrock methods.

**Scores**
- **Readability: 7.** An `Emitter` class and small per-kernel emitters, with comments on the delicate points (repeated reactants). The alternative emitters sit in their own modules, and DESIGN.md:157–245 documents every control. The default LU is a table plus a switch, which is less legible than straight-line code; the documentation is dense.
- **Maintainability: 7.** One integrator source, strict tests and flag stamps (`Makefile:105–108`). But:
  - `"scaling factor"` is silently ignored (`configure.py:54` reads `scaling_factor`).
  - The unknown-species error doesn't name the species (`configure.py:111–112`).
  - `tests/test_configs.py:97` asserts that `configs/ts1.json` equals the translation of MICM's header. DESIGN.md:19–20 documents this, but it contradicts DECISIONS.md:120–123's "edit the configuration".
- **Extensibility: 6.** The specialization options are well factored. But any TS1 edit trips the pin; the harness initializes only a species named `A` for non-TS1 mechanisms (`test_specialized.py:65`); and CB5 is rejected for duplicate labels. The rate-type steps (`DESIGN.md:144–148`) omit `specialize.rate_formula`.

**Exercises**
- **(a) Add a reaction — easy**
  - Add the reaction to the MECHANISM in `configs/ts1.json`; `rate units: mol m-3` applies, so A is in SI units.
  - Run `make bench test-specialized`.
  - `make test-cpu` then reports "supplemental FAIL" from `verify_reader` (`test_configs.py:97–98`), and that test must be edited.
- **(b) Change a rate parameter — easy**
  - Same edit and commands; the same pin, plus a MICM-derived TS1 regression, fail.
- **(c) Another mechanism — moderate**
  - Add a `build/chapman_specialized.c` rule (`specialize.py … --prefix chapman`), add `chapman` to `SPECIALIZED_CONFIGS`, and add `'chapman': '../../mechanisms/chapman/config.json'` to `CONFIGS`.
  - Change the harness initialization; otherwise every Chapman concentration is 0.
  - Run `make test-specialized`.
  - M is ordinary. CB5 fails with `duplicate custom rate label NO2` (`configure.py:155–156`), and its scaling factors would silently become 1.

**Strengths**
1. Generated kernels plug into the library's single integrator.
2. Very strict tests: exact statistics, all methods, plan coordinates, alternative option builds.
3. Duplicate labels fail loudly instead of merging.

**Weaknesses**
1. `"scaling factor"` is silently dropped.
2. The TS1 file is pinned to MICM's translation, and the harness initialization is tied to a species named `A`.
3. The default LU is a table interpreter, harder to read as emitted code.

## Entry D

**Summary.**
- **Reader:** `tools/configure.py` is a strict CAMP reader (order-based units, `rate units`, `parameterization: air density`, `initial state`) that emits a C builder.
- **Generator structure:** `tools/specialize.py` builds a `Plan` (stable rate-type sort, Markowitz ordering, fill, Mozart operations) and an `Emitter`.
- **What it emits:** fixed-stride (default 32) and variable-stride kernels for rates (with an exact-input rate cache), forcing, Jacobian, fused-pivot LU and solves.
- **Integrator:** it then instantiates the library's `src/rosenbrock_step.inc` through macros, so one stage-and-controller implementation serves both solvers. Depfiles cover transitive includes.
- **Tests:** `tests/test_specialized.c` covers TS1, three tutorials, an extra-rates fixture and a make-generated mutated fixture (new species, repeated reactants, density multiplier, emission, shared label). It runs at 1, 7, 16, 32 and 33 cells with default and tight tolerances, and checks:
  - rates, including cache hits and each invalidation input;
  - LU and solves on a synthetic matrix, and singular-pivot detection;
  - 12 integration steps with statistics compared.

**Scores**
- **Readability: 7.** DESIGN.md:195–216 lists the generator's five passes and the Plan/Emitter split, and the code follows it. Comments explain the MICM-matching subtleties (`specialize.py:44–46,215–216`). The style is compact, the rate cache is a long block of C inside Python strings (`specialize.py:167–199`), and `DESIGN.md:192–193` still says generated C "supplies mechanism data, not specialized arithmetic".
- **Maintainability: 9.** A single integrator source, the most systematic tests, settings stamps and depfiles (`Makefile:125–165`), and validation of numbers, units, labels and initial-state names. Minor: `tests/config_equivalence.c:86` pins TS1's counts, and same-named labels are shared silently.
- **Extensibility: 7.** `README.md:108–121` documents generating and wiring another mechanism, and `DESIGN.md:213–215` names every place a new rate type touches. Adding a mechanism to the tests is three edits.

**Exercises**
- **(a) Add a reaction — easy**
  - Edit `mechanisms/ts1.json`; mechanism-level `rate units: mol m-3` applies.
  - Run `make bench test-specialized`.
  - `make test-cpu`'s count check (`config_equivalence.c:86`) and the MICM-derived `tools/test_ts1.py` regression fail.
- **(b) Change a rate parameter — easy**
  - Edit the parameter and run the same commands; only the MICM regression fails.
- **(c) Another mechanism — moderate**
  - Run `.conda/bin/python tools/specialize.py ../mechanisms/chapman/config.json build/specialized/chapman.c --symbol chapman --width 32 --depfile build/specialized/chapman.d`.
  - Add a rule and `chapman` to `SPECIALIZED_NAMES`, and add `DECLARE(chapman)`/`CASE(chapman)` to `tests/test_specialized.c`.
  - Run `make test-specialized`; the harness seeds nonzero values.
  - M is ordinary. In CB5, the two NO2 parameters are set together by name, and Ea is ignored.

**Strengths**
1. A shared integrator, with the generator emitting kernels only.
2. The most systematic tests: cache, singularity, statistics, mutation.
3. A strict reader, plus depfiles and stamps.

**Weaknesses**
1. Compact code with large C-in-Python blocks.
2. Stale sentence in DESIGN.md.
3. Labels are not namespaced, THIRD_BODY is not recognized, and TS1's counts are pinned.

## Entry C

**Summary.**
- **Reader:** `tools/configure.py` is strict: it rejects unknown document keys and uses a per-reaction `concentration units` key.
- **Symbolic layer:** `tools/symbolic.py` holds a `Plan` dataclass for ordering, fill and jobs.
- **Emitters:** `tools/specialize.py` emits forcing, Jacobian, LU (Doolittle with reciprocal pivots by default) and solves in batched SIMD loops; `specialize_rates.py` emits rate formulas.
- **Post-passes that rewrite emitted C:** `split_kernels.py` splits functions by line count behind volatile function pointers; `tiled.py` renames kernels and string-substitutes the `ros3.c.in` template to run packed four-cell tiles.
- **Tests:** `tests/specialized/check.c` compares rates, forcing, Jacobian, factor (on a diagonally dominant matrix), solve and six integration steps (including a night step and a mid-run state copy). It covers TS1, C's own Chapman, Robertson, both tutorials, an edge fixture and an empty mechanism, at 1, 3, 7, 32 and 33 cells and two tolerance regimes, element-wise at 1e-10 plus an infinity-norm check.

**Scores**
- **Readability: 6.** The modules are split by concern and have docstrings. But the C that is finally compiled exists only after `split_kernels.py` and `tiled.py` have rewritten text with regexes and exact-whitespace replacements (`tiled.py:80–90`). `configure.py` is dense one-liners, and `ros3.c.in` hard-codes the coefficients.
- **Maintainability: 7.** Broad, careful tests (including Chapman and an empty mechanism), compiler/specialization/driver stamps (`Makefile:96–153`), and strict validation. Risks: `tiled.py`'s `str.replace` silently does nothing if `ros3.c.in`'s formatting changes, which only the tests would catch; `ros3.c.in` duplicates the method constants kept in `src/parameters.c`.
- **Extensibility: 7.** TS1 edits are configuration-only, `README.md:132–144` documents generating another solver, and the tests need three edits. The per-reaction units key is a trap: a new TS1 reaction without it is converted as molecule cm⁻³.

**Exercises**
- **(a) Add a reaction — easy, with care**
  - Add the reaction with `"concentration units": "mol m-3"` to `configs/ts1.json`.
  - Run `make bench test-specialized`.
- **(b) Change a rate parameter — easy**
  - Edit and rebuild.
- **(c) Another mechanism — moderate**
  - Add a rule `generated/chapman_camp_specialized.c` (calling `specialize.py ../mechanisms/chapman/config.json $@ --symbol chapman_camp_specialized $(SPECIAL_FLAGS)`); the symbol must not clash with the existing `chapman_specialized`.
  - Add the object to `SPECIAL_OBJECTS`, and add the extern and `solvers[]`/`names[]` entries in `check.c`. Run `make test-specialized`.
  - M is ordinary. CB5's duplicate labels get separate slots, but name lookup returns the first (`src/model.c:34`); Ea is ignored.

**Strengths**
1. A strict reader with useful messages.
2. Varied, careful tests, including an empty mechanism and a Chapman case.
3. Specialization stages are modular and each was measured.

**Weaknesses**
1. Fragile text-rewriting passes over generated C.
2. The template duplicates coefficients and controller constants.
3. The per-reaction units key, plus the label and THIRD_BODY gaps.

## Entry B

**Summary.**
- **Reader:** `tools/configure.py` raises `ConfigurationError` and wraps `KeyError`/`TypeError` into readable messages. It uses a document-level `concentration units` key and order-based conversion. It also rescales user-defined scaling factors by reaction order, so an EMISSION is divided by 1e-6·N_A.
- **Generator structure:** `tools/specialize.py` binds a `Plan` (rate-type priority, Markowitz ordering, fill, derivative and pivot schedules) and emits `rate_body` formulas. An `Emitter` splits kernels into bounded helpers with fixed- and dynamic-stride versions, dispatched through `src/specialized.c`.
- **Integrator:** `src/specialized_solver.h` is the three-stage template.
- **Tests:** `tests/test_specialized.c` checks 7 builds (TS1, three tutorials, an all-rate fixture, a zero-order fixture, and the all-rate fixture rebuilt with every non-default option) at widths 1, 3, 32, 37 and 64 and two tolerances. It compares rates, forcing, Jacobian, factorization, solve and six trajectory steps with every statistic equal, plus dt = 0 and dt < 0 edge cases. `tools/test_specialize.py` adds structural and edit-a-reaction checks; `tools/test_configuration.py` exercises the reader's errors.

**Scores**
- **Readability: 8.** The clearest generator of the nine. Docstrings state the reasons (`specialize.py:189–235` on helper splitting and runtime loop bounds), names follow the chemistry, and `DESIGN.md:174–238` walks through each stage and parameter. The second-translation-unit dispatch takes a moment to follow, and the generated C has no chemistry annotations.
- **Maintainability: 8.** Excellent tests and the best input validation (`configure.py:173–179`), plus settings stamps. Minuses:
  - `specialized_solver.h` hard-codes MICM's coefficients and controller constants beside `src/parameters.c`.
  - Include files must be listed per mechanism (`DESIGN.md:235–237`).
  - The user-rate unit convention needed compensating factors in the tests (`tests/test_configuration.c:29–38`).
- **Extensibility: 7.** `README.md:140–153` documents generating another solver. Adding one to the tests needs a rule repeating all the flags, the object list, and `cases[]`. The rate-type steps (`DESIGN.md:93–99`) omit `PRIORITY` and `rate_body` in `specialize.py`.

**Exercises**
- **(a) Add a reaction — easy**
  - Add the reaction to `mechanisms/ts1.json`; the document's `concentration units: mol m-3` applies.
  - Run `make bench test-specialized`.
- **(b) Change a rate parameter — easy**
  - Edit and rebuild.
- **(c) Another mechanism — moderate**
  - Run `tools/specialize.py ../mechanisms/chapman/config.json generated/sp_chapman.c --symbol sp_chapman`.
  - Add a Makefile rule listing the three Chapman files; add `build/sp_chapman.o` to `bin/test_specialized`'s prerequisites and link line; add the `extern` and `&sp_chapman` to `cases[]`. Run `make test-specialized`.
  - M is ordinary. CB5's emission rates are reinterpreted as molecule cm⁻³ s⁻¹ (divided by 6.022e17), the duplicate labels resolve to the first slot by name, and Ea is ignored.

**Strengths**
1. The most readable and best-explained generator.
2. Robust, tested input validation.
3. Strict tests with exact statistics, alternative-option builds and edge cases.

**Weaknesses**
1. The integrator template duplicates the method constants.
2. An unusual unit convention for user-defined rates.
3. Per-mechanism include lists, and the label/THIRD_BODY gaps.

## Entry A

**Summary.**
- **Reader:** `tools/configure.py` normalizes CAMP configurations (mechanism-level `rate units`, order-based conversion, air-density parameterization, initial state) into a dataclass and emits constructor calls.
- **Generator structure:** `tools/specialize.py` builds a `Plan` (Markowitz ordering, fill, Mozart or Doolittle operations as generators of statement lists) and an `Emitter`.
- **What it emits:** rates, forcing, Jacobian, shift, factor and solve kernels, written into a separate `_kernels.c` file in 32-unit blocks with fixed- and runtime-stride versions, plus a header and a JSON metadata file.
- **Integrator:** `tools/specialized_solver.h` is the readable three-stage template; it takes coefficients, controller and error norm from the library's parameter table.
- **Tests:** `tests/specialized.c` is compiled once per mechanism through `SPEC_HEADER`/`SPEC_NAME` macros. It compares rates, forcing, Jacobian, LU, solve and five 30 s solves (tight mode forces rejections with `h_start = 30`) against the generic library for TS1, three tutorials, an all-rate fixture and an empty mechanism, at 1, 3, 7, 32 and 33 cells and two tolerances.

**Scores**
- **Readability: 7.** The Plan/Emitter/`rate_expression` separation is clear, and `DESIGN.md:125–166` explains the outputs and layout. But `rate_expression` dispatches on bare integers (`kind == 5  # branched`, `specialize.py:180–227`) that must be decoded from `RATE_NAMES` in `configure.py`. `Emitter.kernel` builds declarations, definitions and dispatch in one dense method.
- **Maintainability: 8.** A generic per-mechanism test harness (`Makefile:106–111`), flag and generator-option stamps, and an integrator that reads the library's parameter table. Minuses: duplicate labels silently share one parameter (`configure.py:97–100`); missing keys raise raw `KeyError`; branched-rate constants are computed in Python while the generic library computes them in C.
- **Extensibility: 8.** Adding a mechanism to the generator tests is a two-line Makefile edit, and TS1 edits are configuration-only. The rate-type steps (`DESIGN.md:232–240`) don't mention `rate_expression`, and its numeric codes make that edit error-prone.

**Exercises**
- **(a) Add a reaction — easy**
  - Add the reaction to `configs/ts1.json`; the MECHANISM's `rate units: mol m-3` applies.
  - Run `make bench test-specialized`.
- **(b) Change a rate parameter — easy**
  - Edit and rebuild.
- **(c) Another mechanism — easy**
  - Add `build/generated/chapman_specialized.c: $(wildcard ../mechanisms/chapman/*) tools/configure.py tools/specialize.py build/generator-flags .conda/bin/python` with recipe `$(PYTHON) tools/specialize.py ../mechanisms/chapman/config.json $@ --name chapman_specialized $(GENERATOR_OPTIONS)`.
  - Add `chapman` to `SPECIALIZED_EXAMPLES` and run `make test-specialized`. Pattern rules build and run the harness, which replaces zero concentrations with 1e-10.
  - M is ordinary. CB5's NO2/ALD2 share one parameter, and Ea is ignored.

**Strengths**
1. The easiest path to generating and checking a new mechanism.
2. The integrator template uses the library's parameter table and controller.
3. A clear Plan/Emitter split, with JSON metadata recording the species order.

**Weaknesses**
1. Rate kinds are dispatched by bare integers.
2. Labels merge silently, and some missing-key errors are raw tracebacks.
3. The rate-type extension documentation omits the generator.

---

## Comparison

| Rank | Entry | Readability | Maintainability | Extensibility | Sum |
|---:|---|---:|---:|---:|---:|
| 1 | B | 8 | 8 | 7 | 23 |
| 2 | D | 7 | 9 | 7 | 23 |
| 3 | A | 7 | 8 | 8 | 23 |
| 4 | G | 6 | 8 | 7 | 21 |
| 5 | E | 7 | 7 | 6 | 20 |
| 6 | C | 6 | 7 | 7 | 20 |
| 7 | F | 7 | 5 | 6 | 18 |
| 8 | H | 5 | 4 | 5 | 14 |
| 9 | I | 4 | 5 | 4 | 13 |

**Ties.** B, D and A tie at 23, and E and C tie at 20. I broke ties by the criteria in REVIEW.md's order: readability, then maintainability, then extensibility.
- B ranks first on readability.
- D and A are level on readability; D ranks above A on maintainability, because D's generated kernels plug into the generic solver's own integrator source and its tests are the most systematic.
- E ranks above C on readability.

**Reasons for the order.**
- **The top three share the same design:**
  - a validating CAMP reader with order-based units;
  - a generator split into a symbolic plan and an emitter;
  - a template-based integrator;
  - tests that compare every kernel and multi-step trajectories (with solver statistics) across many cell counts and both tolerance regimes;
  - make that regenerates on any change.
- **What separates the top three:** B's code and documentation explain themselves best; D has one integrator source and the most checks; A makes trying a new mechanism almost free.
- **G** matches them on build hygiene and tests (it even mutates TS1 during the build), but its density and many backends make the generator the hardest of the strong group to follow, and its templates repeat the control flow.
- **E** has the cleanest integrator sharing, but silently drops `"scaling factor"`, rejects CB5, and pins TS1 to MICM's translation.
- **C**'s tests and reader are strong, but its text-rewriting post-passes and hard-coded template are fragile.
- **F** is simple and readable, but its non-TS1 tests are largely vacuous and its reader does not follow MICM's units or scaling factor.
- **H and I** trail. H's tests run three identical cells for 1 ms, and its integrator is duplicated in Python strings. I's generator is a single uncommented function, and both readers silently mistranslate CAMP inputs.

**Shared weaknesses.** None of the nine recognizes `"tracer type": "THIRD_BODY"` or `"Ea"`, none namespaces photolysis/emission labels, and none annotates generated C with the reactions it implements. In every entry the generic library shares the reader with the generator, so a generic-vs-generated comparison cannot catch a reader error.

```json
{"scores": {"A": [7, 8, 8], "B": [8, 8, 7], "C": [6, 7, 7], "D": [7, 9, 7], "E": [7, 7, 6],
            "F": [7, 5, 6], "G": [6, 8, 7], "H": [5, 4, 5], "I": [4, 5, 4]},
 "exercises": {"A": ["easy", "easy", "easy"], "B": ["easy", "easy", "moderate"], "C": ["easy", "easy", "moderate"], "D": ["easy", "easy", "moderate"], "E": ["easy", "easy", "moderate"], "F": ["easy", "easy", "moderate"], "G": ["easy", "easy", "moderate"], "H": ["easy", "easy", "moderate"], "I": ["easy", "easy", "moderate"]},
 "rank": ["B", "D", "A", "G", "E", "C", "F", "H", "I"]}
```
