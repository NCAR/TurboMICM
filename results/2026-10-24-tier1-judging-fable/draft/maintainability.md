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
