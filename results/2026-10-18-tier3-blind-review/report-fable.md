All verification is complete. Every quantitative claim I checked in the four SCANS.md and RESULTS.md files reproduced from the committed CSVs, including C's 200-second cap rows just now. Here is the full report.

# Blind review of four Tier 3 entries

Reviewed in the order D, C, B, A. For each entry I read the deliverable documents, the scan CSVs behind the headline tables, the prototypes and their agreement records, the new-method C sources, the scenario drivers, the Makefiles and the test scripts. I recomputed the sums, minima and medians quoted in each SCANS.md and RESULTS.md from the CSVs; every number I checked reproduced. Nothing was built or run. I could not independently confirm timing reproducibility, the GCC runs, the 10,000-cell runs, or the prototype agreement numbers beyond the committed records.

## Entry D

**Summary.** D reruns MICM's five Rosenbrock tables in a new adaptive loop with a per-cell maximum RMS norm, positivity-aware rejection with post-acceptance projection, three step controllers, Jacobian ageing and factor reuse with a step deadband, and adds an independently written variable-order BDF with modified Newton iteration. All five declared points use RODAS3 with a uniform absolute tolerance of 1e-18, 16-cell groups and dynamic scheduling, differing only in relative tolerance. The recommended point is RODAS3 at rtol 3e-4 and targets three digits overall and two per species; its disclosed minima are 3.571 and 2.554 digits. Two looser points serve the two-digit tier and one tighter point gives extra precision.

| Criterion | Score |
|---|---|
| Use of the literature | 7 |
| Rigor of the scans | 8 |
| Soundness of the accuracy argument | 9 |
| Readability | 6 |
| Maintainability | 8 |
| Extensibility | 7 |

**Literature, 7.** `REFERENCES.md` lists six sources and states what each supplied, and `ACCURACY.md` ties its target to Sandu et al. 1997, Dreger et al. 2025 and Bates et al. 2024. The Bates citation is a genuinely good justification of why sub-percent chemistry error matters. The PI controller formula in `DESIGN.md` has no source, and the BDF's error estimator and order selection rest only on the CVODE documentation. Coverage is thin but accurate.

**Scans, 8.** `scans/README.md` and the JSON configurations give a controlled series: family and order in `families-prefix128.csv`, then norm, start step, Jacobian age, step hold and group size one at a time in `structure-prefix128.csv`, controller parameters in `controller-prefix128.csv`, scheduling at 10,000 cells in `large-structure.csv`. I recomputed the family minima and step sums, the Jacobian-age series and the large-structure times, and all match `SCANS.md`. Failed runs and an overlapped batch are kept and labelled. Screening used 128-cell prefixes, which the entry admits do not occupy the machine, and the prototypes in `prototypes/integrators.py` only verify the C code rather than screen accuracy.

**Accuracy argument, 9.** `ACCURACY.md` sets the target, converts digits to RMS bounds, reports per-species diagnostics for major oxidants from `scans/throughput-validation-species.csv`, quantifies projection, and states what is and is not claimed. The minima table matches `scans/throughput-validation.csv`, the 48-hour cases in `scans/stress-final.csv` have references converged to more than ten digits in `scans/stress-reference-generation.csv`, and the narrow margin of the fast point is disclosed.

**Readability, 6.** `src/adaptive.c` holds the whole loop, norm, positivity test and three controllers in 159 lines with several statements per line and one-letter names, and `src/bdf.c` and `src/scenario.c` follow the same style. Comments appear only at a few policy points. `DESIGN.md` is thorough enough that a reader can follow the method, but the code itself is hard to scan.

**Maintainability, 8.** The Makefile uses compiler-specific build directories and a `test-tier3` target that runs `tests/test_adaptive.py` and `tests/test_scenario.py`. Parameters are validated and failures return status codes, every scan carries a metadata file with source hashes, and the option table in `DESIGN.md` matches the parser in `src/scenario.c`. The prototype comparison covers 45 configurations with exact operation counts.

**Extensibility, 7.** `DESIGN.md` gives explicit steps for a new table or family, and they match the code. A new table still requires editing the hard-coded method range in `src/adaptive.c` as well as `src/parameters.c` and `src/scenario.c`. Controllers are integer-coded inside the loop and the new path hard-codes the Mozart factorization, so the callback hooks of the legacy solver do not reach it. Other mechanisms are supported through `tools/configure.py`.

**Strengths.** Honest, well-quantified accuracy policy with species-level diagnostics. Reproducible scan infrastructure with per-run hashes and retained failures. A tested BDF implementation reported as a negative result.

**Weaknesses.** Narrowest method coverage: only MICM's tables plus BDF, with no new tableau, extrapolation or SDIRK. Very dense C that a maintainer will struggle with. Prefix-128 screening and prototypes that verify but do not explore.

## Entry C

**Summary.** C adds a new adaptive driver with seven formulas: the five MICM Rosenbrock tables, the published RODAS3.1 table transcribed from the KPP documentation, and an SDIRK2 with chord iteration. It studies controllers, norms, restart policy, Jacobian lagging with an LU holding band, absolute floors, projection charged to the error estimate, and internal step caps. All four declared points use RODAS3 at rtol 0.003 with 32-cell groups; they differ by absolute tolerance and by a cap on the internal step of 300, 200 or 100 seconds. The recommended point is the 200-second cap with minima of 3.433 and 2.490 digits over nine cases.

| Criterion | Score |
|---|---|
| Use of the literature | 8 |
| Rigor of the scans | 9 |
| Soundness of the accuracy argument | 9 |
| Readability | 8 |
| Maintainability | 9 |
| Extensibility | 8 |

**Literature, 8.** `REFERENCES.md` has five sources, each with a precise statement of use: Sandu for rationale and the lagged-Jacobian warning, the KPP documentation for the RODAS3.1 table, Alexander 1977 for SDIRK, Dreger for H211b, Gustafsson for PI. `DESIGN.md` and `ACCURACY.md` tie their claims to these. The set is small and there is no BDF or extrapolation literature, but nothing cited is misused.

**Scans, 9.** `SCANS.md` names the parent of every structural change and I reproduced every number I checked from `scans/01-families.csv`, `03-structure.csv`, `06-stepcaps.csv`, `07-caps-tolerance.csv`, `09-groups.csv` and `10-copy-elision.csv`, including counter sums for lagged Jacobians. Failed legacy runs are kept with their stderr, each row carries source and binary hashes, and a final review found a sign error in SDIRK's projected error that was fixed, regression-tested in `tests/test_integrator.py` and re-scanned. The prototype comparison in `prototypes/adaptive.py` covers 147 cases with identical counters, with order and failure tests besides.

**Accuracy argument, 9.** `ACCURACY.md` justifies the three-digit target, reports margins over six disclosed and three independent cases that match `scans/final-1000.csv` and `scans/final-stress.csv`, quantifies projection, and labels prefix-only scoring at 10,000 cells. The reference convergence in `scans/stress-reference.csv` is above ten digits. The narrow species margin of fast3 is stated plainly. The step cap is justified only empirically, and the entry says so.

**Readability, 8.** `src/integrator.c` is compact but decomposed into a norm, a Rosenbrock trial, an SDIRK trial and a step loop, with comments that explain non-obvious choices such as the sign convention of the embedded error. `tools/ts1_scenario.c` is straightforward. Method ids such as 43, 431 and 102 are magic numbers, though the `DESIGN.md` table lists them.

**Maintainability, 9.** Tests cover prototype equivalence on two cells with changing rates across three intervals, stiff analytic chains, convergence orders, overflow and Newton failure, the projection regression, and a scenario contract that includes equivalence with the Tier 1 benchmark in `tests/test_scenario.py`. The `validation/` directory and `COMPLETION.md` record hashes for inputs, binaries and outputs. Options are validated in one place and the `DESIGN.md` table matches the parser.

**Extensibility, 8.** A new Rosenbrock table is one case in the factory plus a runner name, and `ch_integrator_create` accepts it automatically through its stage count. A new family needs a trial routine and a dispatch, with the SDIRK special case appearing in three places. The linear solver is replaceable through the public hook and the new path uses it; the controller callback exists only on the legacy path.

**Strengths.** Exemplary provenance and verification, including hash-matched equivalence checks and a regression for a bug found in review. The step-cap study is a clear, measured improvement in global error per unit cost. The most complete test suite of the four.

**Weaknesses.** Moderate method breadth, with no BDF or extrapolation and SDIRK only at order two. All points share one relative tolerance, so the declared family is driven by caps rather than tolerance. Screening timings are single runs and one batch overlapped a compilation, both disclosed.

## Entry B

**Summary.** B adds adaptive Rosenbrock with the five MICM tables, linearly implicit Euler extrapolation of orders two to four with an optional adaptive order, and SDIRK of orders two and three with simplified Newton. It scans controllers, Jacobian age, a step-retention window, stage clipping, norms, absolute floors, group layout, constant-rate caching and step blocking. Two declared points use RODAS3 for the two-digit tier and three use order-four extrapolation; the recommended point is extrapolation at rtol 0.002 with minima of 3.448 and 2.550 digits on the disclosed cases. Rate caching, four-step blocking and dynamic scheduling are on by default and were shown to be bit-identical to the plain driver.

| Criterion | Score |
|---|---|
| Use of the literature | 8 |
| Rigor of the scans | 9 |
| Soundness of the accuracy argument | 9 |
| Readability | 7 |
| Maintainability | 8 |
| Extensibility | 7 |

**Literature, 8.** `REFERENCES.md` cites both Sandu papers, Verwer et al. 1999 for ROS2 and positivity, Alexander with Butcher's tutorial for the SDIRK family and its third-order root, and Hairer and Lubich for extrapolation, stating honestly that only the abstract of the 1988 paper was read. Each claim in `DESIGN.md` points to a source. The PI exponents have no citation.

**Scans, 9.** The scan table in `SCANS.md` covers the widest set of choices of the four entries, and I reproduced the restart cost from `scans/structure64.csv`, the layout sums from `layout256.csv`, the cache and blocking sums from `rate_cache64.csv` and `blocking1000.csv`, and the floor and norm minima from `error64.csv`. Failures are kept, two provenance defects are disclosed and the affected sweep was rerun. The 44 prototype comparisons in `results/prototype_agreement.csv` match counts exactly, and `prototypes/convergence.py` checks orders against a matrix exponential. The stress cases exposed a failing candidate at rtol 0.006 in `scans/stress_gcc.csv`, which was replaced. Screening used 64-cell prefixes.

**Accuracy argument, 9.** `ACCURACY.md` justifies the target, reports a controlled floor scan, charges projection to the error norm and accounts for it exactly on a stiff test in `prototypes/edgecases.c`, and gives per-species minima and the largest individual error. The stress references converge to more than ten digits in `scans/stress_references.csv`, and the entry says the stress cases informed selection rather than calling them holdouts.

**Readability, 7.** `src/adaptive.c` is organized into step routines with comments on policy, but lines carry several statements and one field, `extrap_order`, doubles as the SDIRK stage count. `src/scenario.c` is dense but followable. The formula section of `DESIGN.md` is clear and complete.

**Maintainability, 8.** The Makefile stamps the compiler configuration and `test-tier3` runs the prototypes, C edge cases, rate-dependency tests and `tools/test_scenario.py`, which checks parity with the MICM runner and bit-identical output across cache and blocking settings. Validation transcripts are kept under `results/`, the scan tool refuses to overwrite a CSV, and the option table matches the parser.

**Extensibility, 7.** A new table touches the enum, the parameter factory, two name arrays and an enum bound check. A new family requires branches in the central loop for the shift, the factorization and the trial. The linear backend hook is honoured by the new path, as the edge-case test shows. Generated kernels and a generic evaluator serve other mechanisms.

**Strengths.** Broadest method exploration and a non-obvious recommendation supported by its own data. Driver engineering verified to leave results bit-identical. A stress test that caught and corrected a candidate.

**Weaknesses.** Early provenance defects, though disclosed. Screening on 64-cell prefixes. Dense code with an overloaded option name.

## Entry A

**Summary.** A builds a separate adaptive engine with the five MICM tables, the KPP RODAS3.1 table, and Euler extrapolation of orders two and three, with I, PI and H211b control, three norms, warm starts, positivity-aware rejection, exact reuse of the initial forcing on retries, and Jacobian ageing with a deadband. Four declared points use RODAS3.1 and one uses six-stage RODAS4, all with 32-cell groups and guided scheduling. The recommended point is RODAS3.1 at rtol 3e-4 and atol 1e-18, with disclosed minima of 3.705 and 2.734 digits and 3.923 and 3.194 on three independent 48-hour cases.

| Criterion | Score |
|---|---|
| Use of the literature | 9 |
| Rigor of the scans | 8 |
| Soundness of the accuracy argument | 8 |
| Readability | 8 |
| Maintainability | 8 |
| Extensibility | 7 |

**Literature, 9.** `REFERENCES.md` lists nine sources with page-level use: Sandu for the experiment design and metric caveat, the KPP documentation for stage equations and the RODAS3.1 table, Hairer and Lubich for extrapolation, SUNDIALS for the Jacobian versus factorization distinction, Söderlind for the PI exponents, Dreger for H211b. It is candid about abstract-only access and `ACCURACY.md` separates the supplied pooled metric from the paper's worst-species aggregate.

**Scans, 8.** Scans run on all 1,000 cells. I reproduced the controller findings from `scans/controller.csv`, the lagging trade-off from `scans/reuse.csv`, the RODAS3.1 minima from `scans/rodas31.csv` and `final-tune.csv`, the family table from `scans/family.csv`, and the scheduling study in `scans/scale.csv`. The first family and baseline rows were contaminated by a shared output directory, which `scans/README.md` discloses and `family-corrected.csv` supersedes. The prototype in `prototypes/check.py` passes 48 strict comparisons and 64 cache checks, with H211b sensitivity kept as a separate diagnostic. Extrapolation was tested only at orders two and three and with uncached forcing, so its rejection rests on weaker evidence than B's order-four result would suggest.

**Accuracy argument, 8.** `ACCURACY.md` derives the target from Sandu with the extra digit stated as a design choice, converts the margins to RMS error, and the minima match `scans/operating-points.csv`. The 48-hour cases span NOx factors from 0.02 to 50 with references converged beyond ten digits in `scans/stress-reference-convergence.csv`. There are no per-species diagnostics beyond the worst species, and the 10,000-cell check is a single run scored on the prefix.

**Readability, 8.** The code is formatted by the committed `.clang-format` and `src/integrator.c` is split into clear functions with comments at the policy points. `src/scenario.c` is long but linear, and `include/chem_integrator.h` documents every option field.

**Maintainability, 8.** A configuration stamp forces rebuilds on compiler changes and `make test-methods` runs the prototype, the scenario contract and `tests/test_new_methods.c`. `AUDIT.md` records requirement-by-requirement evidence and sanitizer runs. Scan suites are hard-coded in `tools/scan.py` rather than kept as data, and `src/scenario.c` ignores the return of the parameter initializer for the extrapolation ids.

**Extensibility, 7.** A new table needs only the parameter initializer and a runner name, since the integrator accepts any id the initializer accepts. A new family needs the two extrapolation ids checked in four places. The new path hard-codes Mozart factorization and integer-coded controllers, so the legacy callback hooks do not reach it, and the extrapolation table is fixed at three columns.

**Strengths.** The best use of the literature. Clean, consistently formatted code. Adoption of RODAS3.1 with independently verified order, and a careful roundoff-sensitivity analysis of the aggressive controller.

**Weaknesses.** Incomplete exploration of the alternative family. Early scan contamination. Solver hooks that stop at the legacy engine.

## Comparison

**Research ranking: C, B, A, D.** C and B tie on points. I place C first because the rigor of its experimental process is the strongest of the four: every structural comparison has a named parent, every row carries hashes, the equivalence of driver optimizations was proven by output hashes, and a real defect was found, fixed and regression-tested. B is second for the broadest coverage of method families and structural choices and for a stress test that changed its selection, held back by prefix screening and early provenance slips. A is third: its literature use is the best and its scans are careful and on full inputs, but the alternative family was under-explored and its first scans were contaminated. D is fourth with a clean, honest study whose method coverage is the narrowest and whose code-side prototypes only verify.

**Code ranking: C, A, B, D.** C leads on maintainability and extensibility, with the most complete tests, explicit validation records, and a factory that accepts new tables without touching the loop. A is second for the most readable code and solid tests, with hooks that do not reach the new engine. B is third with well-organized but dense code and a central loop that new families must edit in several places. D is fourth: its documentation and tests are good, but the C sources are the hardest to read of the four and the new path bypasses the legacy callbacks.

```json
{"scores": {"A": [9, 8, 8, 8, 8, 7], "B": [8, 9, 9, 7, 8, 7], "C": [8, 9, 9, 8, 9, 8], "D": [7, 8, 9, 6, 8, 7]},
 "research_rank": ["C", "B", "A", "D"], "code_rank": ["C", "A", "B", "D"]}
```