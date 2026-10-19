# Blind review of the four Tier 3 entries

I read each entry's documents, scan CSVs, plots, prototypes and new C sources, and I built and ran nothing. Where I say a number was "checked", I recomputed it by hand from the entry's own CSVs. I could not check coefficient tables against their published sources, or what the cited papers actually say.

---

## Entry A

**Summary.** A leaves the Tier 1 solver untouched and adds a separate adaptive engine (`include/chem_integrator.h`, `src/integrator.c`). The engine runs MICM's ROS2, ROS3, ROS4, RODAS3 and RODAS4 tables, the RODAS3.1 coefficients from the KPP documentation, and linearly implicit Euler extrapolation of fixed order 2 or 3. Its options are I/PI/H211b control, three error norms, positivity-aware rejection, warm step starts, Jacobian age with an LU deadband, and exact reuse of the initial forcing on retries. The five points are RODAS3.1 at rtol 2e-3 and 1e-3 (level 2), RODAS3.1 at rtol 3e-4 with atol 1e-16 (level 3), RODAS3.1 at rtol 3e-4 with atol 1e-18 (recommended), and six-stage RODAS4 at rtol 1e-4 (level 4); all use 32-cell groups and guided scheduling. The recommended point's minima are 3.705 overall / 2.734 worst-species digits on the disclosed scenarios, and 3.923 / 3.194 on three 48-hour cases A generated itself.

**Scores**

1. **Use of the literature: 8.**
   - `REFERENCES.md` lists nine well-chosen sources: Sandu et al. parts I and II, the KPP Rosenbrock and RODAS3.1 pages, Hairer–Lubich, SUNDIALS CVODE, Söderlind (2002) and Dreger et al. (2025).
   - Each entry says what was taken and admits access limits frankly ("Only the abstract was accessible"). The formulas in `DESIGN.md` (PI exponents, H211b equation 6, reuse rules) are tied to these sources.
   - The accuracy target rests on Sandu (1997) alone, plus one extra digit declared as a design choice.
2. **Rigor of the scans: 8.**
   - `scans/` covers a family × three-tolerance grid and 13 one-at-a-time controller, norm, start, safety, growth and positivity variants (`controller.csv`). It also covers Jacobian age × deadband (`reuse.csv`), H211b, absolute floors and group size, all on full 1,000-cell scenarios, with binary hashes and provisional rows labelled.
   - Every `SCANS.md` number I checked reproduces. These include the family table (the corrected ROS2 rows sum to 23,207 ms), 175,758 vs 313,504 attempts for warm vs cold starts, and 156,360 factorizations with 2.236 species digits at Jacobian age 4 / deadband 0.2.
   - Coverage stays within linearly implicit one-step schemes. `SCANS.md` also calls warm starts a clear improvement without mentioning that cold starts were more accurate (minimum 3.866 vs 3.537 digits).
3. **Soundness of the accuracy argument: 8.**
   - `ACCURACY.md` sets a budget of 0.1 % overall and 1 % per species, and distinguishes the pooled metric from Sandu's worst-species aggregate. It shows the recommended point about five times inside that budget.
   - Every point is checked with GCC, by a 10,000-cell gate (accuracy stated for the prefix only), and on three 48-hour cases. Those cases' references converge to 10.08–10.35 digits (`scans/stress-reference-convergence.csv`), and all five points pass them.
   - Limits are explicit. Only the justification of the target itself is thin.
4. **Readability: 6.**
   - `DESIGN.md` documents every option, the controller formulas and the reuse rules. The code has one statement per line, designated initializers and named workspace arrays.
   - Against that, identifiers are single letters (`v`, `o`, `p`, `s`, `w`), methods are magic integers (31, 43, 12/13, with `method - 10` giving the extrapolation order), and comments are sparse.
   - `prototypes/check.py` is written as hard-to-read semicolon-joined one-liners.
5. **Maintainability: 8.**
   - `make test-methods` runs:
     - 48 strict C/prototype comparisons with exact operation counts, and 64 forcing-cache identity checks;
     - analytic order and invariant tests (`tests/test_new_methods.c`);
     - driver tests, including parity with the supplied MICM runner and rejection of invalid options (`tools/test_scenario.py`).
   - A build stamp forces a rebuild when the compiler or flags change.
   - The documentation matched the code everywhere I checked: option defaults and ranges, controller formulas, the stale-Jacobian retry, and default group shrinking.
6. **Extensibility: 7.**
   - Adding a Rosenbrock tableau takes a case in `src/parameters.c` and a name in `method_id`. `DESIGN.md` gives that procedure, including the required prototype, test and scans.
   - A new method family needs another branch in `ck_integrator_advance`, keyed on magic integers, within a fixed stage capacity.

**Strengths**
- Rigorous, verifiable controlled scans run at full scale, with careful provenance.
- The strongest support for its operating points:
  - the largest recommended margin (about 0.7 digits) and an extra level-4 point;
  - all points pass the independent 48-hour cases;
  - MICM baselines are matched in accuracy at both levels (`micm-level2.csv`, `micm-tuned.csv`).
- The best literature trail, and an exact-work optimization that the data confirms: reusing the retry forcing saves exactly one call per rejection (RODAS4 tier1: 31,974 → 31,706 calls with 268 rejections).

**Weaknesses**
- Narrow family coverage: no Newton-based implicit, BDF, QSSA/Patankar, splitting or precision experiments.
- Terse C, magic method numbers, and a nearly unreadable prototype script.
- The target rests on one source, and the warm-start trade-off between cost and accuracy is reported one-sidedly.

---

## Entry B

**Summary.** B adds `chem_adaptive_solve` (`src/adaptive.c`) with three method families:
- MICM's Rosenbrock tables in a new adaptive loop;
- linearly implicit Euler extrapolation of order 2–4 (frozen Jacobian, Neville tableau, optional variable order that accepts early);
- stiffly accurate SDIRK2/3 with simplified Newton.

Its options are P/PI control, Jacobian age, a step-retention window for factor reuse, stage clipping and three norms. The driver adds rate caching by dependency, four-step blocking and dynamic scheduling. The points are RODAS3 at rtol 0.01 and 0.003 (level 2) and fourth-order extrapolation at rtol 0.004, 0.002 (recommended) and 0.001 (level 3), all with atol 1e-18 and 32-cell groups. The recommended point's minima are 3.448 / 2.550 digits on the disclosed scenarios and 3.929 / 2.848 on three 48-hour stress cases.

**Scores**

1. **Use of the literature: 7.**
   - Six sources, each with its use stated: Sandu I/II, Verwer et al. (1999), Alexander (1977) with Butcher's tutorial for SDIRK, and Hairer–Lubich for extrapolation. Verwer's discussion of ROS2 with approximate Jacobians and of positivity motivated the reuse and clipping experiments.
   - For the methods built, the gaps are notable. The PI controller in `DESIGN.md` has no source, and no controller literature is cited at all.
   - `REFERENCES.md` admits the extrapolation overview was used only "in the final documentation review".
2. **Rigor of the scans: 9.**
   - This is the broadest controlled programme of the four. It covers:
     - the method family, and structure (Jacobian age, P/PI, warm start, clipping, retention window, each across three methods);
     - tolerance curves including SDIRK2/3, fixed vs automatic order, absolute floors and norms, and layout;
     - rate caching and blocking, which give bit-identical outputs;
     - full 1,000-cell work–precision curves.
   - An extrapolation candidate at rtol 0.006 failed a stress case (2.838 digits). That failure changed the selection and is kept in the record.
   - The numbers I checked match: RODAS3 at Jacobian age 4 drops from 2.013 to 1.714 species digits; cold extrapolation restarts raise time from 0.99 to 3.64 s; the atol sweep runs from 3.588 to 3.238 digits.
   - Most screening used 64-cell prefixes on four threads, which B acknowledges.
3. **Soundness of the accuracy argument: 8.**
   - `ACCURACY.md` motivates 3/2 digits from Sandu and keeps transport, mechanism and forcing-discretization errors separate.
   - It quantifies:
     - sensitivity to atol and norm;
     - conservation under projection (a stiff test that accounts for the excess to 1e-10);
     - the largest individual error (0.998 %);
     - accuracy of the key oxidants.
   - It reports the 0.035-digit species margin of the "fast" point honestly.
   - It misses that "accurate" and "recommended" score identically on diurnal-5min: 3.448 / 2.550 in `operating_clang1000.csv`, despite different step counts, and `error64.csv` shows the same plateau across atol and norms. The recommended point therefore adds no margin on the scenario that limits it.
4. **Readability: 6.**
   - `DESIGN.md` gives the most precise description of the formulas of any entry: the Neville recurrence, the SDIRK tableau and embedded weights, the Newton stopping rule, and when caches are invalidated.
   - The C is dense. Lines hold several statements, and `chem_adaptive_defaults` is a 17-value positional initializer.
   - Workspace slots are reached through per-method magic offsets (`w->vectors+9*n`), and SDIRK's stage count is carried in `extrap_order`.
5. **Maintainability: 8.**
   - `make test-tier3` runs:
     - 44 prototype comparisons with exact operation counts;
     - convergence-order checks against a matrix exponential (`prototypes/convergence.py`);
     - tests of invalid input, backend failure and stiff projection, and of rate dependencies;
     - scenario tests against the MICM runner, generated vs generic kernels, and cache/blocking variants (`tools/test_scenario.py`).
   - Binary and output hashes are kept for every run, and `DESIGN.md`'s option table matches the parser.
   - The new solver's C tests sit in `prototypes/`, which is an odd place for them.
6. **Extensibility: 6.**
   - A new tableau is just a table entry.
   - A new family must be threaded through the branches of `chem_adaptive_solve` and through twelve hand-allocated workspace slots whose layout differs by method.
   - New rate types must also be registered with the rate cache.

**Strengths**
- The widest exploration of method families, leading to a choice that was not obvious but is well supported: order-4 extrapolation beats RODAS3 above three digits (`precision_final1000.csv`, `results/work_precision_1000.png`).
- Exemplary provenance and candour:
  - the failed candidate is kept;
  - contaminated scans are disclosed;
  - a score is reused only when the output hash is identical.
- A deep accuracy analysis, and a C MICM baseline run through the same driver. That baseline separates the integrator's gain (2.38×) from implementation effects (7.51× against the C++ MICM).

**Weaknesses**
- Dense, fragile solver code: positional defaults and magic slot offsets.
- The literature is thin for the number of methods built.
- Thin margins at the cheap points, an unexplained 5-minute plateau, and screening done on prefixes.

---

## Entry C

**Summary.** C adds `ch_integrator` (`src/integrator.c`). It runs MICM's ROS2, ROS3, ROS4, RODAS3 and RODAS4 tables, RODAS3.1 from the KPP documentation, and a two-stage L-stable SDIRK with chord iterations. Its options are I/PI/H211b control, three norms, Jacobian lagging with LU holding, a bounded projection charged to the error estimate, and an internal step cap `h_max`. All four points are RODAS3 at rtol 3e-3:
- economy: atol 1e-16, no cap (level 2);
- with atol 1e-18, a step cap of 300 s (fast3), 200 s (balanced, recommended) or 100 s (conservative).

Across the six disclosed and three self-generated 48-hour cases, balanced reaches at least 3.433 overall and 2.490 species digits.

**Scores**

1. **Use of the literature: 7.**
   - Five papers plus a web description, each with its use: Sandu II, the KPP RODAS3.1 table, Alexander (1977) for SDIRK2, Dreger et al. (2025) for H211b and for the tolerances a global model uses, and Gustafsson (1991) for PI.
   - Using Dreger's operational tolerances to frame the target is apt. The base is small, though, and no source addresses why a step cap should help.
2. **Rigor of the scans: 8.**
   - 960 runs in eleven numbered series. Every structural row changes one choice relative to a named parent, all on full 1,000-cell scenarios, with hashes recorded and failures kept.
   - The headline finding is well supported. A 100 s cap gives 3.867 / 3.065 digits in 7.6 s summed (`06-stepcaps.csv`). For comparison, uncapped RODAS3 at rtol 2e-4 gives 3.671 / 2.616 in 11.0 s, and RODAS3.1 gives 3.802 / 2.853 in 9.5 s (`04-workprecision.csv`).
   - A sign error in SDIRK's projected error was found, regression-tested in both the C and the prototype, and rescanned with before/after hashes.
   - The exploratory figure is an undifferentiated scatter in which capped and uncapped runs cannot be told apart.
3. **Soundness of the accuracy argument: 7.**
   - `ACCURACY.md` sets 3/2 digits from Sandu and Dreger, leaving room for splitting and transport errors. It converts atol and the threshold into molecules cm-3.
   - It reports minima over the nine cases with honest caveats: fast3's species margin is only 0.057 digits, and conservative does not reach four digits.
   - The recommended margin (0.433 / 0.490) is the smallest of the four entries, and the reason the step cap helps is never analysed.
   - `RESULTS.md` has no MICM baseline matched in accuracy. It offers only the failing defaults and a single run at rtol 1e-5 that reaches more than 5 digits.
4. **Readability: 7.**
   - The trial routines (`rosenbrock_trial`, `sdirk_trial`) are separate from one shared control loop. State fields are named, and the subtle points carry comments.
   - `prototypes/adaptive.py` is the most readable prototype of the four.
   - On the minus side: an 18-value positional default initializer, magic method IDs (43, 431, 102), and the scenario runner placed in `tools/`.
5. **Maintainability: 8.**
   - `make test-cpu` runs the 680 Tier 1 cases through a ctypes harness, plus `test-methods`:
     - 147 C/prototype comparisons with exact counters, across methods, policies and changing rates;
     - stiff analytic tests, formal-order tests, and tests of overflow, Newton failure and DAE rejection;
     - the projection-sign regression and scenario contract tests (`tests/test_integrator.py`, `tests/test_scenario.py`).
   - `validation/` keeps the gate records, and `DESIGN.md`'s option table matches the defaults in the code.
6. **Extensibility: 7.**
   - `DESIGN.md` states a clear contract: either add a Rosenbrock table, or add a trial routine that writes `new_y` and `error` without touching `y`. The shared loop then supplies control, positivity and acceptance.
   - Selection is still by magic integers inside `ch_integrator_step`.

**Strengths**
- An original structural finding with good evidence: on these scenarios, internal step caps beat tighter tolerances for RODAS3.
- Exemplary transparency in fixing and re-validating a bug.
- A clean verification design: a general, readable prototype, plus analytic, order and failure tests.

**Weaknesses**
- No MICM comparison at equal accuracy.
- An unreadable exploratory figure, and no explanation of why the caps work.
- A narrow set of points (one method, one rtol), the smallest recommended margin, and a thin literature base.

---

## Entry D

**Summary.** D adds `chem_adaptive_solve`. MICM's ROS2, ROS3, ROS4, RODAS3 and RODAS4 tables run in a new loop (`src/adaptive.c`). Beside them is an independently derived BDF with variable step and order (`src/bdf.c`): orders 1–5, coefficients from differentiating a Lagrange interpolant, modified Newton, and guarded matrix reuse. Options are I/PI/H211b control, three norms, step holding and Jacobian age, including a W-method experiment with ROS2. All five points are RODAS3 with atol 1e-18, 16-cell groups and dynamic scheduling, at rtol 0.01 and 0.005 (level 2) and 5e-4, 3e-4 (recommended) and 1e-4 (level 3). The recommended point's disclosed minima are 3.571 / 2.554 digits, and it passes three 48-hour cases.

**Scores**

1. **Use of the literature: 7.**
   - Six sources with stated uses: Sandu II, Verwer et al. (1999), the CVODE documentation for BDF, Dreger et al. (2025) for H211b, and Bates et al. (2024).
   - D uses Bates et al. to argue that chemical effects below one percent matter when choosing a default tolerance. Among the entries, that is the use of the literature most grounded in atmospheric science; I could not check the paper's content.
   - The PI controller and the BDF error estimator have no specific source.
2. **Rigor of the scans: 7.**
   - The series cover:
     - five Rosenbrock formulas and BDF orders 1–5;
     - tolerances, atol, norms and the starting step;
     - PI, six H211b settings, safety and growth;
     - Jacobian age, step holding and ROS2 lagging;
     - grouping and scheduling at 10,000 cells.
   - The values I checked match: the family table, the Jacobian-age minima (3.611 → 3.218), and 2,521 / 13,034 vs 2,196 / 12,093 ms in `large-structure.csv`.
   - `SCANS.md` is mostly prose, and its one table reports accepted steps rather than cost. The CSV shows the cost: about 2.0 s summed for RODAS3 vs 6.3 s for BDF5.
   - Screening used a 128-cell prefix, which leaves most of the ten workers idle.
3. **Soundness of the accuracy argument: 8.**
   - `ACCURACY.md` justifies 3/2 digits from four sources, including why 1 % is too loose as a default.
   - It reports:
     - per-point RMS bounds;
     - accuracy of the key oxidants, with individual maxima (`scans/throughput-validation-species.csv`);
     - candid projection figures (a clip of 1.99e-14 at "balanced", above the scoring threshold).
   - The "controlled" MICM baseline at rtol 3e-4 is tighter than level 3 needs, so the 3.2–11.6× ratios overstate the gain. `RESULTS.md` states this caveat.
4. **Readability: 7.**
   - The Tier 3 code is compact and well partitioned:
     - Rosenbrock in `adaptive.c`;
     - BDF in `bdf.c`;
     - norm and positivity helpers shared between them.
   - Comments explain the BDF coefficients and error estimator, and defaults are set by name.
   - Lines are often dense, methods are integers (2–6, 10), and the prototype mirrors the C closely.
5. **Maintainability: 7.**
   - `make test-tier3` runs:
     - 45 prototype comparisons with exact counts;
     - a stiff multi-lane transfer test;
     - tests of invalid parameters, workspace ownership and forcing discontinuities;
     - a scenario test that recomputes every rate independently to check the driver's photolysis update.
   - Each compiler gets its own build directory.
   - There is no analytic convergence-order test, so an error shared by the C and the mirrored prototype would pass. The Makefile is also `.SILENT`, which hides commands.
6. **Extensibility: 7.**
   - Each family lives in its own file behind a one-line dispatch, `chem_advance` accepts a generic step callback, and `DESIGN.md` explains how to add a table or a family.
   - The driver updates photolysis by scaling precomputed rate constants of user-type reactions labelled "j…". That is exact for TS1, but another mechanism would have to satisfy the same assumption.

**Strengths**
- The only entry to implement and measure a multistep Newton family (variable-order BDF), with a sound explanation of why it loses when forcing restarts the integration every step.
- A confirmed W-method experiment: ROS2 at Jacobian age 16 kept its 3.35 digits while cutting summed time by about 30 %.
- The best-grounded accuracy target, useful species-level diagnostics, and compact, well-structured code.

**Weaknesses**
- `SCANS.md` under-reports the measured quantities, and screening was done on a small prefix.
- Verification has no convergence-order test. Positivity also accepts negatives up to the local error scale without charging them to the error estimate.
- The MICM comparisons are not matched in accuracy, and `RESULTS.md`'s plots show only the chosen RODAS3 points, not the family comparison.

---

## Comparison

**Research (criteria 1–3): B, A, C, D.**
- **B and A tie at 24**, and the order could reasonably flip. I put B first because its scans cover more of the method space: three families, order-4 extrapolation, variable order, and Newton-based SDIRK. Its selection visibly responded to evidence (the rtol 0.006 stress failure), and its analysis of accuracy limits is the deepest.
- A is very close behind. It has the better literature trail, screening at full 1,000-cell scale for every family, and the strongest operating-point margins, including a level-4 point.
- **C and D tie at 22.** C ranks third because its scans are more rigorous: fully documented, run at full scale, and organized around named parents. They also yielded a genuine structural finding (step caps) and showed exemplary transparency about a bug. D's BDF work and its literature-based target are real strengths, but its `SCANS.md` is thin and its screening was on a 128-cell prefix.
- No entry used Python prototypes to explore methods on TS1. In all four, the prototypes mirror the C on 2–3-species problems to verify it.

**Code (criteria 4–6): C, A, D, B.** The spread is small (20–22).
- **C (22)** is first: a clean split between trial routines and one shared control loop, the most readable prototype, and comprehensive tests.
- **A and D tie at 21.** A is placed ahead on maintainability: analytic order tests, driver parity with MICM, and documentation that matches the code. D is more readable, with each family in its own file, but its verification is thinner.
- **B (20)** is last. Its documentation and tests are excellent, but its solver code is the densest and the hardest to extend safely.

```json
{"scores": {"A": [8, 8, 8, 6, 8, 7], "B": [7, 9, 8, 6, 8, 6], "C": [7, 8, 7, 7, 8, 7], "D": [7, 7, 8, 7, 7, 7]},
 "research_rank": ["B", "A", "C", "D"], "code_rank": ["C", "A", "D", "B"]}
```