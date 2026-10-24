# Tier 1 Judging: Judge 2, Claude Fable 5.1

*The second of three judgments of Tier 1's secondary criteria, readability,
maintainability and extensibility, by Claude Fable 5.1 at max effort, on
2026-10-24. Judge 1 was Claude Opus 5.5 ([Tier-1-Judging-Opus]); Judge 3 is
the SSE. The criteria, the exercises and the scale are on
[[Tier 1 Judging|Tier-1-Judging]]; each entrant's design is on its
`Entrant-<run>` page. Nothing in this report bears on the primary ranking,
TS1 wall time, or on the correctness gate, which all nine pass.*

## How this judgment was made

**In the open, like Judge 1.** The session that judged is the one that
orchestrates the experiment. It knew each entry's model, effort and speed.

**What the judge had and had not read.** Before scoring, the judge had read
the project's instructions and the [[Tier 1 Judging|Tier-1-Judging]] page,
both of which show Judge 1's totals and ranks, and the README of Judge 1's
results folder, which lists its procedure, its scores and its findings for
the SSE. The judge had not read Judge 1's report page, its helpers' notes,
exercise reports or test traces, or the "Extension points" sections of the
entrant pages, which summarize Judge 1's exercises; those were read only after
the scores below were fixed and committed, for the comparison at the end of
this page. The scores were formed from the code.

**Procedure.**

1. **The judge read every entry.** For each of the nine, at its final commit
   in a read-only worktree: README.md, DESIGN.md, DECISIONS.md, the Makefile,
   the public header, the benchmark driver, the integrator and its parameter
   sets, the step control, the sparse LU and its plan, the rate evaluation,
   the forcing and Jacobian, the configuration reader and generator, and the
   test runner, plus the generated TS1 source where the entry has one. Notes
   with citations are in `notes/<run>.md` of the results folder.
2. **Five helpers, fresh sessions, the same briefs as Judge 1's.** Three
   worked the exercises (E1 and E2; E3; E4 and E5) and two traced the twelve
   sampled test cases (cases 1 to 6; 7 to 12) in all nine entries, each to the
   briefs in `brief/`, which are Judge 1's with "first judge" changed to
   "second judge" and an added instruction not to read Judge 1's materials.
   The helpers did not know the judge's own ratings. For E5 only, the helper
   copied each entry's Python tools into a scratch folder and ran them on
   MICM's tutorial configurations, on the same configurations with the species
   listed in reverse order, and on a units probe, with a checker of its own
   (`exercises/e5-check/`); nothing in C was built or run, by the helpers or by
   the judge.
3. **Reconciliation.** The judge compared the helpers' ratings and findings
   with its own notes and re-read the code behind every finding that moves a
   score; where a helper's standard differed from the judge's first reading
   the judge chose one standard and applied it to all nine. The differences
   are noted under each exercise.
4. **Scores and report** by the judge: `scores.json` and this page. Compiler
   warnings come from the verification records, which note them only for the
   three Luna entries; nothing was rebuilt for this judgment.

**The rubric** is Judge 1's, stated on [[Tier 1 Judging|Tier-1-Judging]]:
readability and maintainability 1 to 10 each; extensibility twice the mean of
five exercise ratings of 1 to 5, kept to one decimal; equal weights, a total
of 30; ties go to extensibility, then maintainability.
## Scores

Totals out of 30; extensibility is twice the mean of the five exercise ratings
(E1 to E5 in brackets).

| Rank | Entrant | Readability | Maintainability | Extensibility (E1–E5) | Total of 30 |
|---:|---|---:|---:|---|---:|
| 1 | [astra-high](Entrant-astra-high) | 8 | 8 | 8.8 (5, 4, 5, 4, 4) | 24.8 |
| 2 | [sol-xhigh](Entrant-sol-xhigh) | 7 | 8 | 8.4 (5, 4, 5, 3, 4) | 23.4 |
| 3 | [astra-xhigh](Entrant-astra-xhigh) | 7 | 7 | 9.2 (5, 5, 5, 4, 4) | 23.2 |
| 4 | [sol-high](Entrant-sol-high) | 7 | 8 | 8.0 (5, 4, 4, 3, 4) | 23.0 |
| 5 | [astra-max](Entrant-astra-max) | 6 | 7 | 8.4 (5, 4, 5, 3, 4) | 21.4 |
| 6 | [sol-max](Entrant-sol-max) | 6 | 6 | 8.0 (5, 4, 4, 3, 4) | 20.0 |
| 7 | [luna-max](Entrant-luna-max) | 5 | 4 | 7.2 (4, 3, 5, 3, 3) | 16.2 |
| 8 | [luna-xhigh](Entrant-luna-xhigh) | 4 | 5 | 6.4 (5, 2, 3, 3, 3) | 15.4 |
| 9 | [luna-high](Entrant-luna-high) | 3 | 2 | 4.4 (2, 2, 2, 3, 2) | 9.4 |

Three groups. The six Sol and Astra entries are one generic design each, with
the same five parameter sets, a sparse LU with a precomputed plan, a hook for
the linear solver and, in six of them, a hook for the controller; they differ
in how much of that is readable and tested, and in whether generated kernels
duplicate the ordering logic. Ranks 2 to 4 are within 0.4 of a point and the
judge would not defend their order against a different weighting of the same
facts; the gap between rank 1 and rank 2 and the gap between rank 6 and rank 7
are real. The three Luna entries port MICM's cases less faithfully, accept a
NaN step, and convert configuration units with fixed exponents; luna-high
prints its test results from a list and generates wrong reactions for any
mechanism whose species order its generator changes.
## Readability

Can a scientific software engineer follow the solver and match it to MICM? The
judge read each entry's benchmark driver, integrator, step control, sparse LU,
rate evaluation, forcing and Jacobian, reader and generator, and test runner,
in that order, and asked how far each is from MICM's `rosenbrock.inl` and
`lu_decomposition_mozart_in_place.inl`. The common counts (lines, function
sizes, comment share) are those of [[Tier 1 Judging|Tier-1-Judging]].

**astra-high, 8.** The rate constant is a tagged union of per-type structs
with MICM's field names (`include/chem.h:28-58`), so `k0_A`, `kinf_C`, `Fc`
and `N` appear in the formulas (`src/rates.c:64-72`) as they do in MICM. The
integrator is 80 lines with MICM's variable names and structure, the clamp
inside it with a comment on why (`src/rosenbrock.c:146-225`). It has the
highest comment share of the nine (3.9 %) and a DESIGN.md with a file map and
a paragraph per extension. Minor issues: several statements per line in
places (`src/rosenbrock.c:203-207`), and the TS1 path's forcing, sparse
Jacobian and dense Jacobian sit in three files, the sparse derivative plan
inside the solver object (`src/rosenbrock.c:26-81`).

**sol-high, 7.** The integrator mirrors MICM line for line in 104 lines, with
`new_function`, `reject_last`, `reject_more` and the `c/h` coefficients where
MICM has them (`src/rosenbrock.c:28-131`), and comments mark each place where
a correspondence with MICM matters (`src/rosenbrock.c:18,98`,
`src/mechanism.c:147-149`). DECISIONS.md (114 lines) records every format
extension and test-port choice. Two real costs: the rate parameters are
positional, `ChemRate.p[24]` with a comment key in the header
(`include/chem.h:20-25`), so a formula reads `p[0] * exp(p[2] / t) *
pow(t / p[3], p[1])` (`src/rates.c:27-30`); and the style packs several
statements on one line throughout (`src/rosenbrock.c:110`,
`src/mechanism.c:174-176`).

**sol-xhigh, 7.** One statement per line, Allman braces, named rate fields
(`include/chem.h:54-67`, `src/rates.c:62-80`), and an AUDIT.md that maps the
prompt's requirements to files. Against it: the lowest comment share of the
six Sol and Astra entries (1.0 %; the 201-line solver has two comments,
`src/solver.c:1-3,28`); `process.c` holds the forcing, the Jacobian, both LU
factorizations, the solves and the dispatch, which neither its name nor
DESIGN.md says (`src/process.c:108-300`); and `ck_model_create_with_order`
runs 297 lines with the same three-line cleanup repeated four times
(`src/model.c:248-544`, `:320-324,333-336,346-349,360-363`).

**astra-xhigh, 7.** The public structs show every array the solver uses
(`include/chem.h:107-175`), the integrator is 87 lines with MICM's names
(`src/rosenbrock.c:34-120`), the rate union has named fields
(`include/chem.h:27-55`), and DESIGN.md (174 lines) gives a file map, a
tuning table and a paragraph per extension point. Against it: a compact
style of `if (...) { a; b; }` lines (`src/rosenbrock.c:38,56-57,101-110`), an
LU plan replayed as a switch on four numbered op kinds of which only the
first is commented (`src/linear.c:53-67,78-103`), and a 156-line
`chem_model_create` (`src/model.c:53-209`). MICM's two DA sets are named
`CHEM_RODA3` and `CHEM_RODA4` (`include/chem.h:213`).

**astra-max, 6.** DESIGN.md is the fullest of the nine (211 lines: ownership,
layout, methods, a tuning table, an extension section per exercise with the
parameter slot table). The code is denser than the document: `ch_rate.p[32]`
positional slots, the branch selected by `p[4] == 0.0` and prepared values in
`p[5]`, `p[6]` (`include/chem.h:35-41`, `src/rates.c:32-36,67`); two-space
indentation with several statements per line (`src/rosenbrock.c:99-103`,
`src/state.c:12-18`); negative statuses returned beside nonnegative builder
ids (`include/chem.h:17-22,75-86`); and a rate-type order kept as a literal
table by enum position (`src/inspect.c:18-22`). The 86-line integrator itself
reads well (`src/rosenbrock.c:34-119`).

**sol-max, 6.** The 79-line integrator is clear (`src/rosenbrock.c:114-192`)
and the loader-plus-kernels generation is explained in DESIGN.md, but the
style is the densest of the nine (`src/rosenbrock.c:9,37,73,165-170`,
`src/model.c:34-35,47`), `struct ChemModel` and `struct ChemState` are public
with every internal field, which the benchmark and the generated loader write
directly (`include/chem.h:94-122,142-149`; `src/bench.c:127-130,168`;
`tools/configure.py:202,233`), two rate fields are overloaded ("B stores
k_r", "A stores the scaling factor", `src/rates.c:68,70`) with a rename map
in the generator to match (`tools/configure.py:203`), and README.md (41
lines) and DESIGN.md (49 lines) are short.

**luna-max, 5.** One statement per line and long names, but 0.3 % of lines
carry a comment (the lowest share of the nine), the rate parameters are
`p[16]` positional with no key in C (`src/kinetics.c:70`; the key is Python's
`KINDS`, `tools/generate_c.py:16-30`), nearly every function returns -1 for
every failure (`src/solver.c:945-1001`), and the 154-line integrator runs a
loop over the group's cells inside the stage loop, cell by cell, with the
factorization copying the matrix (`src/solver.c:644-797`, `:326-356`), a
structure that does not match MICM's. DESIGN.md's four-item extension list is
accurate.

**luna-xhigh, 4.** `src/chem.c` is 1,744 lines holding rates, parameter sets,
Markowitz, the symbolic build, forcing, Jacobian, two LU factorizations, the
solves, the norms, the ODE and DAE integrators and the finite-difference
utilities. The rate slots are irregular, Troe in `p[5..12]`, JPL19's reference
temperatures at `p[18], p[19]`, branched at `p[13..16]`, k_r at `p[17]`
(`src/chem.c:52-151`), with no key in C. The 207-line integrator also does the
symbolic build, the allocation and two transposes of the state
(`src/chem.c:1061-1267`), and the library prints diagnostics to stderr
(`:1127,1132,1157,1207`). The names are descriptive and the documents honest.

**luna-high, 3.** Two Rosenbrock integrators: a generic single-cell one and
the TS1 group one with the three-stage coefficients and the controller
constants written inline (`src/solver.c:489-690,767-945`; `:796-804,908`).
Two rate evaluators that overlap (`src/mechanism.c:140-225,243-295`). A
727-line header mixing matrix views, CSR views, species records, constraints,
dense LU and two solver APIs. Dense `n x n` Jacobian storage with lists of the
nonzero columns beside it (`src/solver.c:90-91,254-303`). Comments on 0.4 % of
lines. DESIGN.md says the factorization "can be replaced behind a
linear-solver interface" (DESIGN.md, Extension points), which does not exist.
## Maintainability

Can the code be changed safely? The judge weighed the tests first (are MICM's
cases ported case by case, and does each id stand for a check that runs for
it?), then the build, the separation of library, benchmark and tooling, logic
kept in two places, error handling, documents against code, and dead code.
The twelve-case sample of [[Tier 1 Judging|Tier-1-Judging]] was traced in all
nine entries by two helpers (reports T1 and T2 in the results folder); the
judge re-read the code behind every grade that moves a score. Compiler
warnings were taken from the verification records, which note them for the
three Luna entries only.

**sol-high, 8.** Each of the 680 ids is looked up by name: a Python callable
keyed by `Suite.Test` or by the full id, else a C check by name in
`bin/test_rates --list`, else `FAIL case port pending`
(`tests/run_cases.py:41-64`). The twelve sampled cases grade 11 F and 1 P (the surface error asserted at rate construction rather than store construction).
After the 680 lines the runner verifies the reader, compares
configuration-built mechanisms with hand-built ones, and reruns TS1's stiff
cell and the 100-cell fixture at group sizes 1, 7, 32 and 128 against MICM's
reference to 1e-10 (`tests/supplemental.py:10-36`); the generator runs on
four configurations in the build. Statuses carry a reason, allocation failures
unwind (`src/mechanism.c:165-168`), and the integrator inspects every
candidate value so that a NaN state cannot be accepted through the floored
norm (`src/rosenbrock.c:98-101`). Against it: the diagonal Markowitz ordering
is implemented twice, in Python for generation and in C as an API
(`tools/configure.py:167-205`, `src/linear.c:150-180`); and a rate type
appended to the enum is rejected by a range check and reads past a family
table unless both are edited (`src/mechanism.c:30,213-221`).

**sol-xhigh, 8.** The runner is in C: each id is split into suite and test
and offered to 23 dispatchers that return pass, fail, N/A or "not handled",
and an id no dispatcher knows prints `FAIL port not implemented`
(`tests/runner.c:15-89`). The sample grades 9 F and 3 P. Ten reader unit
tests run first (`tools/test_configuration.py`), the three tutorial
configurations are generated in the build and compared with hand-built
descriptors, the ordering lives once, in C at run time, and failure paths
free everything (`src/model.c:33-86,541-543`). The NaN case propagates to
`CK_NAN` (`src/solver.c:38,134-138`). Against it: the 297-line model
constructor and its repeated cleanup, and a `(unsigned)r->type > CK_LAMBDA`
check that a type appended after `CK_LAMBDA` must move (`src/rates.c:47`).

**astra-high, 8.** Every check is invoked through a macro that prints the
id's own PASS or FAIL line; the Python driver rejects duplicate or unknown
ids, adds the eight tutorial comparisons, and prints `FAIL case port pending`
for any id without a report (`tools/run_tests.py:18-34`). No N/A: 680 PASS
with both compilers. The sample grades 12 F. The reader has
forty-nine validation checks, the generator runs on four configurations and
the results are compared with hand-built mechanisms (`Makefile`, `bin/test_config`,
`bin/test_config_extra`), the ordering lives once in C, and a 19-value status
enum names each failure (`include/chem.h:12-16`). Against it: the test binary's
25-file compile line is written twice in the Makefile, each group builds its
own copy of the Jacobian plan (`src/rosenbrock.c:26-81`), and the benchmark
recomputes the rate constants inside the timed loop (`src/ts1_bench.c:154`),
which the tests do not notice but a maintainer should.

**astra-xhigh, 7.** The strongest test plumbing of the nine: `TEST(prefix,
Suite, Name)` macros are scanned into a registry that rejects an id absent
from the case list or registered twice (`tools/test_registry.py:8-17`), and
the runner prints `FAIL not yet ported` for the rest (`tests/main.c:13-24`);
no N/A; `tools/check_ts1.py` runs nine thread, group, schedule and evaluator
combinations against the 100-cell reference. The sample grades
12 F. The weakness is structural: the generated kernels' flat
indices come from a Python copy of the Markowitz ordering and of the fill
(`tools/generate.py:13-34,61-75`), the C model computes its own
(`src/linear.c:128-155`, `src/model.c:122-137`), and nothing compares the two
when the kernels are installed (`src/model.c:201-202`). The shipped TS1 agrees
with MICM, and the three tutorials and Chapman are generated with kernels
in the build, but a new mechanism on which the two orderings differed would
run wrong silently. Also: the post-solve clamp is a separate call that the
caller must remember (`src/bench.c:194,205`).

**astra-max, 7.** Ids register through `@case` and `@group` decorators that
assert the id is in the list and unregistered; the runner prints `FAIL case
port not yet implemented` otherwise (`tests/run_cases.py:18-34`); a C program
run before the cases installs Doolittle through the linear-solver hook and a
counting controller, checks the call counts against the statistics and the
solution against the default (`tests/native_cases.c:33-65`), the only test
of the extension hooks as such among the nine. The sample grades
9 F and 3 P. The ordering lives once in C; the three examples are
generated and compared with hand-built mechanisms. Against it: the positional
slots couple Python's `FIELDS` lists and `RATE_NAMES` order to the C enum and
`p[]` layout (`tools/configure.py:17-33,168`, `include/chem.h:28-41`), a
literal order table by enum position (`src/inspect.c:18-22`), a per-call
allocation in `ch_update_rates` (`src/state.c:135`), and a model left
half-built after a failed finalization (DESIGN.md says to rebuild it).

**sol-max, 6.** Ids register by full id and the runner fails closed
(`tests/run.py:35-57`); the generator runs on three examples, a validation
configuration and TS1, and README says generic and generated kernels are
compared. The sample grades 10 F and 2 P. Against it: the generated
kernels' indices come from a Python copy of the ordering and the fill, checked
at load only by the nonzero count (`tools/configure.py:118-157,230`), while
DESIGN.md says the indices "are verified against the finalized model"; the
reaction order is sorted in Python and again in C by two parallel lists that
must stay aligned (`tools/configure.py:10-12,101`, `include/chem.h:20-24`,
`src/model.c:186-191`); the public structs let the benchmark and the loader
write internals; the Mozart factorization mallocs a buffer on every call and
on failure writes NaN into the diagonal (`src/linear.c:17-18`).

**luna-max, 4.** `bin/test_cpu` dispatches each id by suite, and inside a
suite several ids share one check: all `troe_rate_constant` ids run one fixed
Troe evaluation (k0 2, kinf 4, T 300 K) that never uses MICM's parameters
(`tests/test_cpu.c:21-48,113-115`), and `DAESolveEnforcesConservation` falls
through to a default branch that solves an unrelated three-species constraint
(`tests/test_cpu.c:3332-3426`); the Python driver then overwrites eight rows
with other programs' exit codes (`tests/run_cases.py:129-144`). The sample
grades 2 F, 3 P and 7 S. The norm's `fmax` floor turns a NaN into 1e-10, so
a step with NaN values is accepted (`src/solver.c:80,767-776`). Backward Euler
has no step reductions and restores the state on non-convergence
(`src/backward_euler.c:35-117`), not MICM's method. The reader converts units
with fixed exponents for Troe and ternary only (`tools/generate_c.py:88-94`).
Objects are shared between compilers (README: `make clean` between them);
gcc-15 reports three warnings in the library and the tests add 38 with each
compiler.

**luna-xhigh, 5.** Most suites run per-id checks (`bin/test_state <test>`,
`bin/test_lu <test>`, `bin/test_constraints <test>`...), but every
rate-constant id, the reaction-rate store, `solver_configuration` and the
example programs take one binary's exit status for all their ids
(`tools/test_cases.py:74,145-150,155-157,166-169`). The sample grades
3 F, 4 P and 5 S. The floored norm accepts a NaN step (`src/chem.c:1045,1206`).
The CAMP reader sets every initial concentration to zero and converts units
with fixed exponents for Troe and ternary only (`tools/generate.py:71,131-136`).
The library prints to stderr. The symbolic build and all allocation happen
inside the integrator call, so inside the benchmark's timer. The Makefile
repeats a compile and a link rule for each of 24 test programs. No compiler
warnings with either compiler. The per-id suites, the TS1 gate run at one and
four threads inside the tests, and the absence of a check that reports an
unrelated result put it a step above luna-max.

**luna-high, 2.** `tools/test_cases.py` runs nothing: it prints `PASS` with
a canned reason for each of 510 literal ids once `bin/test_core`, a single
program that stops at its first failing check, and `tests/test_codegen.py`
have exited 0 (`tools/test_cases.py:7-639,640-778`; `Makefile` test-cpu). No
id selects a check; a failure anywhere suppresses all 680 lines. The sample
grades 1 F, 2 P, 8 S and 1 M, every one printed from the list. Two integrators and two rate evaluators duplicate the
numerics (`src/solver.c:489-690,767-945`; `src/mechanism.c:140-225,243-295`).
The generator writes each reaction term's input index into a species table it
has permuted (`tools/generate_c.py:102,126`; `tools/read_config.py:83-86`), so
a CAMP mechanism gets the wrong species wherever the Markowitz order differs
from the input order; TS1, translated by `tools/import_ts1.py` without
indices, is unaffected. The floored norm accepts a NaN step
(`src/solver.c:907,916`). The library and tests compile from one line each,
without objects. No compiler warnings.
## Extensibility

Each exercise was worked through on paper for all nine entries, by the judge
for the integrator, the controller and the linear algebra while reading the
code, and by helpers for the full edit-site accounting (`exercises/E1-E2.md`,
`E3.md`, `E4-E5.md` in the results folder); the judge re-read every site that
a rating depends on. The ratings follow the scale of
[[Tier 1 Judging|Tier-1-Judging]]; the standard applied to each exercise is
stated under its table so that the nine ratings can be compared.

One remark on E1 before the tables: RODAS3 (Sandu et al., 1997) is,
coefficient for coefficient, MICM's four-stage differential-algebraic set
(`include/micm/solver/rosenbrock_solver_parameters.hpp:309-356`), which every
entry already carries. The exercise was therefore judged as what adding RODAS3
would have required: a new four-stage set with its own tables, flags and order,
selectable beside the five.

### E1. Add a Rosenbrock parameter set

| Entry | Ease | Where the set goes | The integrator TS1 uses | Hazard |
|---|---:|---|---|---|
| sol-high | 5 | a `case` of `chem_rosenbrock_parameters` (`src/parameters.c`), or a `ChemRosenbrock` filled by the caller (`include/chem.h:104-114`) | generic to 6 stages (`src/rosenbrock.c:30,61`) | capacity six; benchmark picks 3 at `src/bench.c:187` |
| sol-xhigh | 5 | a `case` of `ck_rosenbrock_parameters_init` (`src/parameters.c`), or a filled `ck_rosenbrock_parameters` (`include/chem.h:275-288`) | generic to 6 (`src/solver.c:46,90`) | capacity six; `src/parameters.c` is generated by `tools/import_parameters.py` from MICM's five |
| sol-max | 5 | a `case` of `chem_rosenbrock_parameters`, assigned to `ChemSolver.parameters` (`include/chem.h:175-183,193`) | generic to 6 (`src/rosenbrock.c:115-116,136`) | capacity six; `tools/translate_parameters.py` regenerates the file |
| astra-high | 5 | a value of `ChemRosenbrockMethod` and a `case` (`include/chem_solver.h:38`, `src/rosenbrock_parameters.c`) | generic to 6 (`src/rosenbrock.c:148,173`) | capacity six; `tests/test_dae.c:119` sweeps the enum by range |
| astra-xhigh | 5 | a value of `chem_ros_method` and a `case` (`include/chem.h:213`, `src/parameters.c`) | generic to 6 (`src/rosenbrock.c:38,71`) | capacity six; the file is the committed output of `tools/import_parameters.py` |
| astra-max | 5 | a `case` of `ch_rosenbrock_default`, or a filled `ch_rosenbrock_parameters` (`include/chem.h:170-181`) | generic to 6 (`src/rosenbrock.c:35,68`) | capacity six; a hand-built table must also set `controller` (`src/rosenbrock.c:36`) |
| luna-xhigh | 5 | a constructor beside the five (`src/chem.c:229-321`), passed to `chem_integrate_group_steps_with_parameters` | generic to 6 (`src/chem.c:1166-1191`) | the benchmark calls the defaults-only wrappers (`src/ts1_bench.c:343,364`); the DAE integrator is a second copy of the loop |
| luna-max | 4 | a constructor beside the five (`src/solver.c:452-604`), passed to `chem_solver_create` | generic to 6 (`src/solver.c:666,709`) | the copy that feeds a stage reusing the previous forcing sits inside the new-evaluation branch (`src/solver.c:713-728`), so a set with two consecutive reuse stages gets a stale stage vector |
| luna-high | 2 | a value of `chem_rosenbrock_parameter_set` and a branch (`include/chemistry.h:619-625`, `src/solver.c:170-252`), reaching the generic single-cell path only | the TS1 group integrator hard-codes three stages and its tables (`src/solver.c:796-804,859-887`) | a second integrator to keep in step; the order `1.0/3.0` is literal (`:908`) |

Standard: a public parameter struct read per call by one generic stage loop
with a documented capacity of six, and a named selector that is one `case` or
one constructor, is a 5. luna-max's loop has a latent flaw that only a new
set's reuse flags can reach (4). luna-high's TS1 path hard-codes the method (2).

### E2. Replace the step-size controller

MICM's controller is the clamped factor `0.9·err^(-1/order)` and two rules
after rejections (`rosenbrock.inl:230-283`). A PI controller needs the previous
accepted step's error and size, and it replaces those rules. In every entry the
accept/reject branch and the rejection rules stay in the integrator; the
difference is whether a callback exists, what it receives, and where its state
can live.

| Entry | Ease | Hook | What it receives | Where a PI controller's state lives |
|---|---:|---|---|---|
| astra-xhigh | 5 | `chem_step_controller`, an argument of `chem_ros_solve` (`include/chem.h:224-226,233-234`; `src/rosenbrock.c:49,97`) | `h`, the error, the parameters, a per-call context | the caller's per-group context (`src/bench.c:193,204` pass NULL today) |
| sol-max | 4 | `ChemSolver.controller` (`include/chem.h:189-195`; `src/rosenbrock.c:160`) | the error and the parameters | the per-group `ChemSolver` (`src/bench.c:143`) |
| astra-high | 4 | `chem_solver_set_controller` (`include/chem_solver.h:46-47,56`; `src/rosenbrock.c:13-14,199`) | the error and the parameters | the per-group `ChemSolver` (`src/ts1_bench.c:133-134`) |
| sol-xhigh | 4 | `step_factor` field (`include/chem.h:281-282`; `src/solver.c:131-132`), tested (`tests/test_integrators.c:182-191`) | the error and the parameters | the parameter struct, one per group if the caller makes it so (`src/bench.c:293-314` share one) |
| astra-max | 4 | `parameters.controller` (`include/chem.h:169,177-178`; `src/rosenbrock.c:94`), tested (`tests/native_cases.c:26-53`) | the error and the parameters | the parameter struct (`tools/ts1_bench.c:153-170` share one; DESIGN.md says per-state contexts are needed) |
| sol-high | 4 | `ChemStepController` field (`include/chem.h:93-94,110-111`; `src/rosenbrock.c:102-104`) | the error, the order, the safety factor and the bounds | the parameter struct (`src/bench.c:187-202` share one) |
| luna-max | 3 | `parameters.controller` (`include/kinetics.h:172-173,207`; `src/solver.c:761-770`) | the error, the order, the safety factor | one `callback_context` shared by the controller, the factor and the solve callbacks and by every group that the one solver runs in parallel (`src/solver.c:981-997`), with no group identity passed |
| luna-xhigh | 2 | none: `chem_step_size_factor` is a static function (`src/chem.c:1048-1051`) | | the accept/reject logic is in two integrators (`src/chem.c:1201-1231,1476-1502`); no object persists between calls |
| luna-high | 2 | none: literal constants in both integrators (`src/solver.c:656,908`) | | two integrators; the TS1 one hard-codes the order too |

Standard: a hook that receives the step size and the error with a per-call
context, so that an exact PI controller is one function and one line in the
caller, is a 5; a factor-only hook with a context is a 4, since a signature
change in one module and two lines in the accept/reject branch complete it; a
hook whose context cannot be made per group without changing the API is a 3;
no hook, with the rules in two integrators, is a 2.

### E3. Swap the linear solver

| Entry | Ease | Interface | What the factor callback receives | Hazard |
|---|---:|---|---|---|
| astra-xhigh | 5 | `chem_state_set_linear_backend(state, factor, solve, context)` (`include/chem.h:199-205`; `src/state.c:75-94`), used by Rosenbrock, backward Euler and constraint initialization | the state (its public LU plan: map, row starts, columns, diagonal) and the packed matrix | a state copy carries the binding and its context (`src/state.c:63`) |
| astra-max | 5 | `ch_state_set_linear_solver(state, factor, solve, context)` (`include/chem.h:204-207`; `src/state.c:72-83`), tested with Doolittle through the hook (`tests/native_cases.c:33-65`) | the `ch_sparse` plan (public CSR getters), cells, the packed matrix, a scratch row | `ch_state_set_lu` afterwards silently drops the callbacks (`src/state.c:69`) |
| astra-high | 5 | `chem_solver_set_linear_solver` with `ChemLinearSolver {factor, solve, context}` (`include/chem_solver.h:48-57`; `src/rosenbrock.c:15-16,170,189`) | the `ChemSparse` plan, cells, the packed matrix | a pair with one NULL member is silently ignored (`src/rosenbrock.c:16`) |
| sol-xhigh | 5 | `ck_linear_operations {factor, solve, context}` through `parameters.linear`, dispatched by `ck_dispatch_factor/solve` (`include/chem.h:268-273,280`; `src/process.c:285-300`) | the model (its map and row starts), cells, the packed matrix, the LU enum | one parameter struct, hence one context, shared by every group in the benchmark |
| luna-max | 5 | `parameters.factor/solve/callback_context` (`include/kinetics.h:174-181,208-210`; `src/solver.c:657-658,705,708,737`) | the public `chem_sparse_pattern`, alpha, the negative Jacobian; writes a separate LU | per cell, scalar; one context for all threads |
| sol-high | 4 | `ChemLinearSolver {factor, solve}` set per state by `chem_set_linear_solver` (`include/chem.h:95-98,209`; `src/linear.c:181-186`) | the mechanism (`chem_jacobian_index` maps storage), cells, the packed matrix | no context pointer: an external solver's handle must be static or thread-local |
| sol-max | 4 | `ChemSolver.factor/.solve` (`include/chem.h:190-198`; `src/rosenbrock.c:37,80,98,152`) | the model, cells, stride, the packed matrix | no context pointer and `factor` returns void: a singular pivot can surface only as NaN in the norm |
| luna-xhigh | 3 | the enum `chem_lu_method` with two values (`include/chem.h:29-32,308`), tested in the step loop (`src/chem.c:1152-1154`) | | a third method edits the integrator's dispatch, the validator and the symbolic plan, whose species order depends on the method (`src/chem.c:502-508`); the plan struct is private |
| luna-high | 2 | none: static `factor`/`solve_factored` called in both integrators (`src/solver.c:305-344,591,624,854,877`) | | DESIGN.md:113 claims an interface the code lacks |

Standard: a run-time function-pointer pair with a context, on the TS1 path,
documented, is a 5; the same without a context is a 4; an enum over built-in
methods whose extension edits the integrator is a 3; no selector, with the LU
hard-coded in two paths, is a 2.
### E4. Add a rate-constant type

The type for definiteness, k = A·exp(C/T)·(T/300)^B·(1 + D·P/101325), needs
a configuration name, a parameter list with defaults, a unit rule (A scaled
by the molecule-per-cm3 factor to the power order − 1), a C tag, a C
parameter layout, an evaluation branch and a test. Every entry keeps the rate
types in both Python (a type table and a unit rule) and C (an enum and an
evaluator), so no entry can score 5 here; 4 is the ceiling for a
configurable type, and a 3 records a hidden coupling beyond the irreducible
change.

| Entry | Ease | Python | C | Hidden coupling | Callback type for C-side rates |
|---|---:|---|---|---|---|
| astra-high | 4 | one `TYPES` entry `(ctype, member, defaults)` and a unit branch (`tools/configure.py:16-27,137-153`); the emitter is generic (`:224`) | an enum value before the `CHEM_RATE_TYPE_COUNT` sentinel, a union member, defaults, a case (`include/chem.h:22-58`; `src/rates.c:22-44,73-119`) | none: the counting sort and the range checks size themselves from the sentinel (`src/mechanism.c:55-67`, `src/rates.c:47`) | `CHEM_CALLBACK` with context, set by name (`include/chem.h:25-27,95`) |
| astra-xhigh | 4 | a branch of `rate_initializer` and the unit lists (`tools/configure.py:95-108,120-158`) | an enum value, a union member, defaults, a case (`include/chem.h:21-55`; `src/rates.c:8-24,50-100`) | none that bites an appended type: `chem_rate_type_order` returns the enum value for unknown types (`src/process.c:64-69`), `chem_rate_prepare` has no range check | `CHEM_LAMBDA` + `chem_model_set_lambda` (`include/chem.h:24,26,245`) |
| sol-high | 3 | a `RATE_PARAMETERS` entry and the unit branch (`tools/configure.py:48-60,143-151`) | an enum value, a default, a case in `p[24]` slots (`include/chem.h:6-10`; `src/rates.c:32-123`) | `valid_data` rejects `type > CHEM_CALLBACK` and `rate_family[]` is a 12-entry table indexed by type (`src/mechanism.c:30,213-221`) | `CHEM_CALLBACK`, settable by name (`include/chem.h:11-12,121`) |
| sol-xhigh | 3 | `TYPES`, a `DEFAULTS` field whose name must match the C struct, the unit set (`tools/configure.py:16-24,140-150,209`) | an enum value before `CK_LAMBDA`, a field in the one flat `ck_rate` struct, a case (`include/chem.h:39-67`; `src/rates.c:66-122`) | `ck_rate_prepare` rejects `type > CK_LAMBDA` (`src/rates.c:47`); `rate_priority[]` is positional (`src/model.c:239-243`) | `CK_LAMBDA` + `ck_model_set_lambda` (`include/chem.h:52,234`) |
| sol-max | 3 | `KINDS`/`C_KINDS` parallel lists, `fields`, `rename`, the unit branch (`tools/configure.py:10-12,89-101,203-204`) | an enum value, fields, `chem_rate`, `chem_rate_value`, the range check in `chem_model_reaction` (`include/chem.h:21-42`; `src/rates.c:11-81`; `src/model.c:92`) | the position in `KINDS` is the reaction sort key in Python and the enum value is the sort key in C (`tools/configure.py:101`; `src/model.c:186-191`), so the two lists must stay parallel | `CHEM_LAMBDA` + `chem_model_lambda` (`include/chem.h:23,136`) |
| astra-max | 3 | `RATE_NAMES` position and a `FIELDS` list, the unit branch (`tools/configure.py:17-33,141-155,168`) | an enum value, `p[32]` slots, defaults, a case, the range check, the `type_order` table (`include/chem.h:28-41`; `src/rates.c:10-84`; `src/model.c:106`; `src/inspect.c:18-22`) | `type = RATE_NAMES.index(kind)` must equal the C enum ordinal; a second order table by position | `CH_LAMBDA` + `ch_model_rate` (`include/chem.h:30,32,94`) |
| luna-max | 3 | a `KINDS` entry and a unit branch that has no notion of reaction order (`tools/generate_c.py:16-30,88-94`) | an enum value before the sentinel, a case in `p[16]` slots (`include/kinetics.h:8-21`; `src/kinetics.c:30-112`) | `valid_mechanism` rejects `kind > CHEM_CALLBACK` (`src/solver.c:615`); the slots have no key in C | `CHEM_CALLBACK` with context (`include/kinetics.h:19,23-26,84-85`) |
| luna-xhigh | 3 | `RATE_ENUM`, `PARAM_INDEX`, `type_map`, the shared `fields` dict, the unit branch (`tools/generate.py:14-34,74-88,115-136`) | an enum value, a case with new slots (`include/chem.h:12-24`; `src/chem.c:65-151`) | all 24 slots are assigned (`tools/generate.py:26-34`); the shared defaults apply Arrhenius' `D = 300` to every type | `CHEM_RATE_CALLBACK` with context (`include/chem.h:23,188,200-201`) |
| luna-high | 3 | `KINDS`, `TYPE_ORDER`, per-kind slot code in the generator and the units in the reader (`tools/generate_c.py:8-37,148-198`; `tools/read_config.py:75-79`) | an enum value, a case in `chem_reaction_rate` and, for a pressure-dependent type, in `chem_rate_constants` (`src/mechanism.c:140-225,252-283`) | two Python files and two C evaluators; pressure is rebuilt from density in one of them | `CHEM_LAMBDA`, a function pointer without context (`include/chemistry.h:25,28-30,52`) |

Standard: the irreducible change (a Python type table, a unit rule, a C
enumerator, a C evaluator case) in one obvious place per side, with nothing
else that depends on enum position, is a 4; one or more hidden couplings
(positional tables indexed by type, parallel lists across the language
boundary, range checks on the last enumerator, slot layouts without a key, or
the logic split over two files on one side) is a 3. The helper rated luna-max 4
on its edit count; the judge kept 3 for the slot layout and the unit rule
without a reaction order, which a new type with the standard rule would have
to introduce.

### E5. Use another mechanism

What a user runs, what must be edited by hand, whether a generated-kernel path
serves the new mechanism, whether the tests run the generator on anything but
TS1, and defects. The helper ran copies of each entry's reader and generator
on MICM's user-defined and Robertson tutorials, on both with the species
listed in reverse order, and on a three-reaction units probe, and compared the
emitted reactions with the configuration by species name
(`exercises/e5-check/`, 45 runs).

| Entry | Ease | What a user runs | Hand edits | Kernels | Tests on other mechanisms | Defect or hazard |
|---|---:|---|---|---|---|---|
| sol-high | 4 | `tools/configure.py CONFIG OUT.c --symbol NAME` (README) | a Makefile rule; the benchmark binds `ts1_data` (`src/bench.c:14`) | none; the generic path serves all | four configurations against hand-built mechanisms (`tests/test_configs.py:71-80`) | format extensions a new file must adopt (`third body`, `rate units`) |
| sol-xhigh | 4 | `tools/configure.py CONFIG OUT.c --symbol NAME` | a Makefile rule; the benchmark binds `ts1` (`src/bench.c:14`) | none | three, with rates, forcing, Jacobians and trajectories (`tests/test_configuration.c:54-134`) | `concentration units` and `parameterization: air-density` are this entry's keys |
| sol-max | 4 | `tools/configure.py CONFIG OUT.c --name NAME`, then `chem_load_NAME(1)` | a Makefile rule; the benchmark calls `chem_load_ts1(1)` (`src/bench.c:103`) | yes, for any mechanism; indices from a Python copy of the ordering, checked at load by the nonzero count only (`tools/configure.py:230`) | three examples, a validation configuration and TS1, generic and generated (`tests/tutorials.c:175-191,261-292`) | the duplicated ordering (consistent on all five probe inputs) |
| astra-high | 4 | `tools/configure.py CONFIG OUT.c --symbol NAME`, then `NAME_create`/`NAME_initialize` | a Makefile rule; the benchmark calls `ts1_create` (`src/ts1_bench.c:12-13`) | none | three tutorials at one and seven cells plus `config_extra.json` (`tests/config_equivalence.c:38-66`) | `rate units` and `parameterization` keys are this entry's |
| astra-xhigh | 4 | `tools/generate.py CONFIG OUT.c --symbol NAME [--no-kernels]` | a Makefile rule; the benchmark binds `ts1_mechanism` (`src/bench.c:12`) | yes, for any kinetic mechanism; indices from a Python copy of the ordering and fill with no run-time check (`src/model.c:201-202`) | three tutorials with kernels on and off and both LU methods, plus Chapman (`tests/config_equivalence.h:12-76`) | the duplicated ordering (consistent on all five probe inputs) |
| astra-max | 4 | `tools/configure.py CONFIG OUT.c --name NAME`, then `ch_build_NAME()`/`ch_init_NAME()` | a Makefile rule; the benchmark calls `ch_build_ts1` (`tools/ts1_bench.c:16-17`) | none | three (`tools/tutorials.c:234-262`) | `rate units` per mechanism |
| luna-max | 3 | `make mechanism CONFIG=path` | none for a library; the symbol is always `chem_model` (`tools/generate_c.py:326`), so one mechanism per link | none | three parsed, two compared with hand-built tables (`tests/test_example_runtime.c`) | wrong rate constants for any reaction whose order differs from the tutorial's: fixed unit exponents (`tools/generate_c.py:88-94`); exactly one `MECHANISM` entry allowed (`:71-72`) |
| luna-xhigh | 3 | `make cpu CONFIG=path` or `tools/generate.py CONFIG OUT.c` | none; the symbol is fixed, `chem_generated_mechanism()` (`tools/generate.py:298-301`) | none | two tutorials (`tests/test_example_descriptor.c:44-147`) | the same fixed unit exponents (`tools/generate.py:131-136`); every initial concentration emitted as 0 (`:71`) |
| luna-high | 2 | `tools/read_config.py CONFIG OUT.json`, then `tools/generate_c.py OUT.json OUT.c` | no Makefile rule but TS1's; fixed symbol `chem_generated_mechanism` | none | two tutorials, checked only for token counts and the two conversions (`tests/test_codegen.py:94-129`) | **wrong reactions**: the reader's input positions are emitted as indices into a Markowitz-ordered table (`tools/read_config.py:83-86`; `tools/generate_c.py:102,126`); on the species-reversed user tutorial all ten reactions name wrong species; the same unit defect (`tools/read_config.py:75-79`) |

Standard: a documented command, one Makefile rule, a symbol chosen by option,
unit conversion by reaction order, species indices that survive any input
order, and tests that run the generator on the tutorials against hand-built
mechanisms, is a 4 (a new mechanism is a generate, compile and link cycle,
not a run-time hook, so not a 5). A reader that converts units with the
tutorial's fixed exponents gives wrong coefficients for other reaction orders,
3. A generator that emits wrong reactions whenever its ordering differs from
the input's, with a two-step pipeline and no rule, 2. The probe's values are
in the helper's report: the six Sol and Astra readers scale a bimolecular
Arrhenius A by the factor to the first power; the three Luna readers leave it.
## The test sample

The twelve cases of [[Tier 1 Judging|Tier-1-Judging]], traced in every
entry from `make test-cpu` to the check that decides the id's verdict. Grades:
**F** faithful (MICM's scenario and assertions, through the entry's C API);
**P** partial (fewer or weaker assertions, different inputs, looser
tolerance); **S** surrogate (a different or shared check stands in); **M**
mapped only (a verdict with no check behind it). One rule was applied
throughout: when an id's verdict is the exit status of a binary shared by many
ids, the grade is S even if a matching check exists inside that binary, and
the check is cited in the helpers' reports. For luna-high every grade carries
the further note that the id's line is printed from a list once one program
has exited 0, so no grade there describes a check selected by the id.

| | Case | sol-high | sol-xhigh | sol-max | astra-high | astra-xhigh | astra-max | luna-max | luna-xhigh | luna-high |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | `analytical_rosenbrock_integration:AnalyticalExamples.Robertson` | F | F | F | F | F | F | P | P | F† |
| 2 | `rosenbrock:RosenbrockSolver.RejectedStepAlphaMatchesInPlaceSolver` | F | F | F | F | F | F | P | F | P† |
| 3 | `lu_decomposition_mozart_in_place:LuDecompositionMozartInPlace.RandomMatrixVectorOrdering` | F | P | F | F | F | F | S | P | S† |
| 4 | `linear_solver_in_place:LinearSolverInPlace.DiagonalMarkowitzReorder` | F | F | F | F | F | P | S | S | P† |
| 5 | `chapman_ode_solver:ChapmanMechanismHardCodedAndGeneral.factored_alpha_minus_jac` | F | F | F | F | F | F | S | F | M† |
| 6 | `regression_test_rosenbrock_dforce_dy:RegressionRosenbrock.VectorJacobian` | F | P | P | F | F | P | S | S | S† |
| 7 | `troe_rate_constant:TroeRateConstant.CalculateWithAllArguments` | F | F | F | F | F | F | S | S | S† |
| 8 | `reaction_rate_store:ReactionRateConstantStore.SurfaceMissingDiffusionCoefficientThrows` | P | F | F | F | F | F | F | S | S† |
| 9 | `state:State.SettingConcentrationsWithInvalidArguementsThrowsException` | F | P | P | F | F | P | P | P | S† |
| 10 | `analytical_backward_euler:AnalyticalExamples.Oregonator` | F | F | F | F | F | F | F | F | S† |
| 11 | `external_model_constraints:ExternalModelConstraints.DAESolveEnforcesConservation` | F | F | F | F | F | F | S | P | S† |
| 12 | `rate_constants_user_defined_example_by_hand` | F | F | F | F | F | F | S | S | S† |
| | **F / P / S / M of 12** | 11 / 1 / 0 / 0 | 9 / 3 / 0 / 0 | 10 / 2 / 0 / 0 | 12 / 0 / 0 / 0 | 12 / 0 / 0 / 0 | 9 / 3 / 0 / 0 | 2 / 3 / 7 / 0 | 3 / 4 / 5 / 0 | 1 / 2 / 8 / 1 (all †) |

† luna-high: printed from a list; the grade is for the nearest check in `tests/test_core.c`.

What the grades mean, in the judge's reading of the code behind them:

- **Six entries port the sample case by case.** The P grades are small:
  sol-xhigh's random LU case uses four lane widths that only change the seed and
  a dominant diagonal instead of lognormal values; sol-max and astra-max test
  the Chapman Jacobian against finite differences rather than MICM's table for
  the `VectorJacobian` id; sol-high asserts the surface error at rate
  construction rather than store construction; the invalid-arguments state case
  loses MICM's second fault (one value per species on three cells) in
  sol-xhigh, sol-max and astra-max.
- **luna-max shares checks across suites.** All 52 LU and linear-solver ids
  run one dense 3 x 3 and 12 x 12 check (`tests/test_cpu.c:4313-4316,2558-2609`);
  the four Troe ids evaluate one fixed Troe case that never uses MICM's
  parameters (`:21-48,113-115`); `DAESolveEnforcesConservation` has no branch
  and falls to an unrelated three-species constraint step (`:3332-3426`); the
  Robertson case runs one method where MICM runs five, with a relative bound
  of 1e-8 for MICM's 1e-6 (`:1002-1059`). The Oregonator and the surface
  store case are faithful.
- **luna-xhigh is split.** The state, LU, constraint, Chapman, analytical and
  regression suites run a program with the test's name; the rate-constant
  suites (30 ids), the reaction-rate store, `solver_configuration` and the
  example programs take one binary's exit status each
  (`tools/test_cases.py:74,145-150,155-157,166-169`). Its Troe check inside
  `bin/test_rates` is faithful; its store check evaluates a rate on a species
  without a diffusion coefficient rather than building a store.
- **luna-high prints its results.** `tools/test_cases.py` holds 510 literal
  ids and prints `PASS` with a reason chosen by prefix for each, provided
  `bin/test_core` and `tests/test_codegen.py` exited 0 earlier in the recipe
  (`tools/test_cases.py:7-639,640-778`). The checks behind the sample range
  from a faithful five-method Robertson (`tests/test_core.c:4484-4547`) to
  none at all for `factored_alpha_minus_jac` (the reason printed is "C
  Rosenbrock alpha-minus-Jacobian transforms checked"; the nearest check is
  MICM's `lin_solve` case, `:4231-4256`). A failure anywhere in the 90-odd
  checks suppresses all 680 lines instead of failing one id.

The full traces, with MICM's assertions for each case and every entry's
check cited, are `tests/T1.md` and `tests/T2.md` in the results folder.
## Findings for the SSE

These bear on the primary gate or on the record, not on the scores above.
Each was confirmed by the judge in the code.

- **luna-high's 662 PASS lines are printed from a list.** `tools/test_cases.py`
  runs no check; it prints `PASS` with a prefix-chosen reason for each of 510
  literal ids once `bin/test_core`, a single program that returns at its first
  failing check, has exited 0 (`tools/test_cases.py:7-639,640-778`;
  `Makefile`, target `test-cpu`). Of the twelve sampled ids, one
  (`chapman_ode_solver:…factored_alpha_minus_jac`) has no corresponding check in
  `tests/test_core.c`; the others have checks of varying fidelity that gate the
  program as a whole.
- **luna-max reports some ids from checks that are not theirs.** All 52 LU and
  linear-solver ids run one dense check (`tests/test_cpu.c:4313-4316`); the four
  Troe ids evaluate one fixed Troe case with none of MICM's parameters
  (`:21-48,113-115`); `external_model_constraints:…DAESolveEnforcesConservation`
  has no branch and falls to an unrelated constraint step (`:3332-3426`).
- **luna-xhigh takes one binary's exit status for 30 rate-constant ids and the
  reaction-rate store ids** (`tools/test_cases.py:74,145-150`).
- **All three Luna entries accept a step whose error norm is NaN.** The norm
  floors with `fmax(…, 1e-10)`, which returns 1e-10 for a NaN argument, so the
  finiteness test after it never fires and the candidate is accepted
  (luna-max `src/solver.c:80,767-776`; luna-xhigh `src/chem.c:1045,1206`;
  luna-high `src/solver.c:907,916`). sol-high floors the same way but inspects
  every candidate value (`src/rosenbrock.c:98-101`); the other five propagate
  the NaN through a conditional or a custom maximum. TS1 does not reach this.
- **luna-high's generator names the wrong species for CAMP inputs.**
  `tools/read_config.py:83-86` records each term's input position and
  `tools/generate_c.py:102,126` uses it as the index into a species table it has
  permuted. The E5 helper's run on MICM's user-defined tutorial with the species
  listed in reverse order emitted all ten reactions with wrong species
  (`exercises/e5-check/run.log`). TS1 is translated by `tools/import_ts1.py`
  without positions and is unaffected.
- **The three Luna readers convert units with fixed exponents.** Troe and
  ternary coefficients are scaled by fixed powers of the molecule-per-cm3
  factor, and Arrhenius, tunneling and branched coefficients are not scaled at
  all (luna-max `tools/generate_c.py:88-94`; luna-xhigh
  `tools/generate.py:131-136`; luna-high `tools/read_config.py:75-79`). A
  bimolecular Arrhenius rate given in molecule cm-3 units stays as written: the
  E5 helper's probe shows luna-high leaving A = 1e-12 for an A + B reaction.
  MICM's tutorials happen to need only the Troe and ternary conversions.
- **luna-max's backward Euler is not MICM's.** One Newton loop over the whole
  interval per cell, with a dense LU; on non-convergence it restores the state
  and returns -1; no step reduction sequence (`src/backward_euler.c:35-117`).
- **Generated kernels and the run-time plan.** astra-xhigh installs generated
  kernels whose flat indices come from a Python copy of the ordering and fill,
  with no check against the model's own (`src/model.c:201-202`;
  `tools/generate.py:13-34,61-75`); sol-max checks only that the nonzero
  counts agree (`tools/configure.py:230`). Both agree with MICM on TS1 and on
  the tutorials; the hazard is for a new mechanism.
- **astra-high recomputes every rate constant inside the timed loop**
  (`src/ts1_bench.c:154`), work that MICM's driver and the other entries do
  before the warm-up; its TS1 time therefore includes it.
- **luna-xhigh builds its symbolic plan and allocates all work arrays inside
  the integrator call**, once per group inside the benchmark's timer
  (`src/chem.c:1080-1098`, `src/ts1_bench.c:364-365`).
## Records

- **Results folder:** [`results/2026-10-24-tier1-judging-fable`](https://github.com/NCAR/TurboMICM/tree/main/results/2026-10-24-tier1-judging-fable):
  `README.md` (procedure), `scores.json`, `notes/<run>.md` (the judge's reading
  notes with citations), `brief/` (the three briefs and the independence
  note), `exercises/E1-E2.md`, `E3.md`, `E4-E5.md` (the helpers' exercise
  reports) with `exercises/e5-check/` (the checker, its inputs and its log),
  and `tests/T1.md`, `T2.md` (the twelve traces).
- **Judge 1:** [Tier-1-Judging-Opus](Tier-1-Judging-Opus) and
  [`results/2026-10-20-tier1-judging-opus`](https://github.com/NCAR/TurboMICM/tree/main/results/2026-10-20-tier1-judging-opus).
- **The entries:** branches `runs/<model>/<effort>` of
  [NCAR/TurboMICM](https://github.com/NCAR/TurboMICM); the final commits are
  sol-high `1d944f4`, sol-xhigh `171f54a`, sol-max `001fb64`, astra-high
  `4055d7c`, astra-xhigh `fbd365c`, astra-max `0dc7f67`, luna-max `da2d722`,
  luna-xhigh `ada55d5`, luna-high `d3c2480`. Citations on this page are
  `path:line` at those commits.
