# The five exercises, on paper, for the three Terra entries

Worked by Judge 1 itself, from each entry's final commit, read-only, at
`~/NSF-NCAR/tier1-work/verify/<run>/clang/<run>/`, to the brief and rating rules of
`results/2026-10-20-tier1-judging-opus/brief/EXERCISES.md` and the standards stated in that
folder's exercise reports (`exercises/E1-E2.md`, "Conventions used here"; `E3.md`; `E4-E5.md`,
"Rating rules applied to all nine"), so that the three are rated as the nine were. Citations are
`path:line` from the entry's root. Edit-site counts read "total (library and Python + benchmark +
tests)", one site per function, table or declaration. For E5, copies of each entry's reader and
generator were run in scratch (E5.0). Nothing was built or run inside an entry.

## Summary

| Entry | E1 | E2 | E3 | E4 | E5 | Extensibility |
|---|---:|---:|---:|---:|---:|---:|
| terra-max | 5 | 5 | 5 | 3 | 4 | 8.8 |
| terra-xhigh | 5 | 3 | 2 | 3 | 4 | 6.8 |
| terra-high | 5 | 4 | 2 | 3 | 4 | 7.2 |

---

## E1. Add a Rosenbrock parameter set

RODAS3 is already in all three, as MICM's four-stage DA set (terra-max
`src/rosenbrock.c:65-75`; terra-xhigh `src/rosenbrock.c:86-97`; terra-high
`src/rosenbrock.c:74-84`), so, as for the nine, E1 adds one more set of its shape. All three check
the stage count against 6 on every solve and copy a reused stage's raw forcing forward before its
c-terms (terra-max `src/rosenbrock.c:142,192-193`; terra-xhigh `:566,602`; terra-high
`:148,230-231`), so a set with two reused stages in a row is safe in all three.

| Entry | Ease | Edit sites | Extension point | Main hazard |
|---|---:|---|---|---|
| terra-max | 5 | 5 (2 + 1 + 2) | a public constructor returning `ChemRosenbrock` (`src/rosenbrock.c:19-95`; `include/chem.h:187-191`), passed to `chem_solve` | the analytical loops are bounded by a literal 5 beside a method switch (`tests/reference_nonlinear_integration.c:13-22,208`), so an appended set is silently skipped |
| terra-xhigh | 5 | 5 (2 + 1 + 2) | a public constructor returning `ChemRosenbrockParameters` (`src/rosenbrock.c:38-114`; `include/chem.h:111-115`), consumed by both integrators | the tests choose sets by ternary chains bounded by 5 (`tests/test_analytical.c:146,183`) |
| terra-high | 5 | 5 (2 + 1 + 2) | a public constructor returning `ChemRosenbrockMethod` (`src/rosenbrock.c:27-107`; `include/rosenbrock.h:43-48`), passed to `chem_rosenbrock_solve` | in the Oregonator test, slot 5 of the method loop means backward Euler (`tests/test_nonlinear.c:161,168-170`), so a sixth set appended to the array would be run as backward Euler |

The benchmark edit is the same in all three, an option to pick the set where the three-stage set
is hard-coded (terra-max `src/ts1_bench.c:196`; terra-xhigh `src/ts1_bench.c:170`; terra-high
`bench/ts1_bench.c:164`). DESIGN.md says how in each, and the code agrees (terra-max
`DESIGN.md:101`; terra-xhigh `DESIGN.md:77-80`; terra-high `DESIGN.md:25-27`). **Ease 5** for all
three.

---

## E2. Replace the step-size controller

The rating rules are Judge 1's: 5 when the hook already receives what a PI controller needs (the
step size, and the previous accepted step's size and error, kept by the caller) and TS1 goes
through it; 4 when TS1 goes through one place that must grow; 3 when there is no hook and the
accept/reject logic is written out in more than one integrator; 2 when TS1 has its own copy.

| Entry | Ease | Edit sites | Extension point | Main hazard |
|---|---:|---|---|---|
| terra-max | 5 | 3 (1 + 1 + 1) | `ChemStepController(error, step, parameters, context)` returning the step factor, in `ChemRosenbrock` (`include/chem.h:98,109-110`), called by `chem_solve` (`src/rosenbrock.c:217-218`); NULL selects MICM's formula | the context sits in the parameter struct, which the benchmark shares across its groups (`src/ts1_bench.c:196,201,211`), so a stateful controller needs one struct per group |
| terra-high | 4 | 5 (3 + 1 + 1) | none: the controller is the block in `chem_rosenbrock_solve` (`src/rosenbrock.c:250-262`), which the grouped solve also runs (`src/rosenbrock_group.c:125-126`) | a hook (a typedef and an options field, `include/rosenbrock.h:15-29`) and the history must be added in that one integrator |
| terra-xhigh | 3 | 6 (4 + 1 + 1) | none: the controller is written out in `chem_rosenbrock_solve` (`src/rosenbrock.c:617-630`) and again in `chem_rosenbrock_solve_group` (`:706-718`) | two copies of the accept/reject rules to change together |

**terra-max: ease 5.** The hook receives the error and the step size and returns the factor
(`src/rosenbrock.c:217-220`); the solver keeps MICM's rules around it: acceptance when
`error < 1` or `h < h_min` (`:222`), the clamp to `[h_min, h_max]` and no growth after a
rejection (`:226-227`), and the tenfold cut from the third rejection in a row, which overrides
the factor (`:232`). A controller with state keeps the previous accepted error and step in its
context. Edits: an H211b function, `main` giving each group its own `ChemRosenbrock` copy with
its own context (`src/ts1_bench.c:196-213`), and a test. Hazards: acceptance must be inferred
from `error < 1`, which misses the forced acceptance at `h < h_min`, as in astra-xhigh; no test
sets a controller. DESIGN.md says the callback "supplies the next step-size factor; a null
callback uses the default formula" (`DESIGN.md:59-61`), as the code does. Backward Euler has its
own MICM rules (`src/backward_euler.c:24-79`).

```
+typedef struct { double h_prev, err_prev; } pi_state;          /* one per group */
+static double h211b(double err, double h, const ChemRosenbrock *p, void *c) {
+    pi_state *s = c; double k = p->order, f = p->safety / pow(err, 1/k);
+    if (err < 1 && s->h_prev > 0)
+        f = p->safety * pow(err * s->err_prev, -0.25/k) * pow(h / s->h_prev, -0.25);
+    if (err < 1) { s->h_prev = h; s->err_prev = err; }
+    return fmin(p->factor_max, fmax(p->factor_min, f)); }
 /* src/ts1_bench.c:196-213: params[g] = p; params[g].step_controller = h211b;
    params[g].step_controller_context = &pi[g]; pi[g] reset before each chem_solve */
```

**terra-high: ease 4.** There is no hook. DESIGN.md names the block to change and says it is
shared by single-cell and grouped solves (`DESIGN.md:186-191`), which is true: the grouped
integrator builds a block system and calls `chem_rosenbrock_solve` (`src/rosenbrock_group.c:
107-126`). Edits: a typedef and an options field with its default (`include/rosenbrock.h:15-29`;
`src/rosenbrock.c:127-137`), the block and the history it needs (`src/rosenbrock.c:180-181,
250-262`), the benchmark (`bench/ts1_bench.c:164-170`) and a test: a few edits in one module,
with no second copy, the general scale's 4. Backward Euler has its own rules
(`src/backward_euler.c`).

**terra-xhigh: ease 3.** There is no hook, and the error-to-step formula and the accept/reject
rules are written out twice, in the single-cell integrator the benchmark runs
(`src/rosenbrock.c:617-630`) and in the grouped one (`:706-718`); backward Euler's single and
grouped paths have their own (`:726-860`). DESIGN.md says so: "change the error-to-step
calculation and the accept/reject logic in both `chem_rosenbrock_solve` and
`chem_rosenbrock_solve_group`" (`DESIGN.md:80-83`). Edits: a typedef and a parameter field
(`include/chem.h:68-77`), both integrators, the benchmark and a test.

Order: terra-max (5) > terra-high (4) > terra-xhigh (3).

---

## E3. Swap the linear solver

The rating rules are Judge 1's: what separates the entries is the hook and its context: a hook
per state or per solver (5), one shared context (4), a hook with no context (3), and no hook (2).
Every benchmark must be edited to choose a new solver.

| Entry | Ease | Edit sites | Extension point | Main hazard |
|---|---:|---|---|---|
| terra-max | 5 | 3 (1 + 1 + 1) | `ChemLinearSolver {factor, solve}` installed per state with its own context by `chem_set_linear_solver` (`include/chem.h:122-135,251-254`; `src/state.c:295`), dispatched by `chem_linear_factor` and `chem_linear_solve` (`src/linear.c:239-251`) | the system passed is the negative Jacobian in the plan's filled pattern and the generator's ordering; a solver that wants another pattern keeps its own |
| terra-xhigh | 2 | 6 (4 + 1 + 1) | none: `ChemLUAlgorithm` in the parameters (`include/chem.h:67,70`) picks a branch inside the static `factor` (`src/rosenbrock.c:287-331`), and `substitute` is static (`:333-347`) | dense n × n storage with a symbolic pattern; `factor` and `substitute` serve five integration paths (`:517-545,588,605,673,690,767,770,830,832`) |
| terra-high | 2 | 7 (5 + 1 + 1) | none: the integrator chooses among sparse or dense, Mozart or Doolittle, in place or copied, inline (`src/rosenbrock.c:208-216,234-235`); backward Euler repeats the choice (`src/backward_euler.c:64-70`) | two integrators to connect, and an option to validate (`src/rosenbrock.c:157`) |

**terra-max: ease 5.** The callbacks receive a `ChemLinearSystem` view: CSR row starts and
columns, an n × n lookup, the algebraic-row mask, and values stored `[entry·cells + cell]` for the
group's cells (`include/chem.h:122-131`), with alpha to add only to differential diagonals
(`:124`); the default replays the in-place Mozart plan (`src/linear.c:239-251`), and Doolittle and
a searching Mozart are supplied (`:253-345`). The integrator has one factor and one solve call
(`src/rosenbrock.c:177,200`). Edits: the new pair (in the library or the caller), one
`chem_set_linear_solver` per group in the benchmark (`src/ts1_bench.c:186-191`), and a test.
Hazards: the pattern carries the Mozart fill; a rejected step recomputes the negative Jacobian
because the default factors in place (`src/rosenbrock.c:237`), harmless for an out-of-place
solver. DESIGN.md describes the pair, the view and the alpha rule (`DESIGN.md:66-73`), as the code
has them.

```
+static int my_factor(ChemLinearSystem *s, double alpha, void *c) { /* add alpha, factor into c */ }
+static int my_solve(const ChemLinearSystem *s, double *rhs, void *c) { /* solve with c */ }
 /* src/ts1_bench.c:186-191, after each group's state is created: */
+    ChemLinearSolver mine = {my_factor, my_solve};
+    chem_set_linear_solver(groups[g].state, &mine, &contexts[g]);
```

**terra-xhigh: ease 2.** No hook: the LU is an enum choosing between two branches of a static
function inside the integrators' file, on a dense n × n matrix with precomputed row and column
lists (`src/rosenbrock.c:23-36,287-347`). A replacement adds an enum value, branches in `factor`
and `substitute`, its state in the workspace, the benchmark option and a test. DESIGN.md says the
routines "are the boundary for replacing the linear solver" and that "the interface should be
made pluggable before swapping it" (`DESIGN.md:46-48,84-88`), which matches the code.

**terra-high: ease 2.** No hook: `ChemLuAlgorithm` (`include/linear.h:4`) and the optional plan
select one of eight factor calls and two solve calls written out in `chem_rosenbrock_solve`
(`src/rosenbrock.c:208-216,234-235`), and again in backward Euler (`src/backward_euler.c:64-70`).
DESIGN.md: "implement the factor and solve operations for that plan, add an algorithm choice to
`ChemLuAlgorithm`, and connect it at the factor/solve calls in both integrators"
(`DESIGN.md:78-80`), which matches the code.

Order: terra-max (5) > terra-xhigh = terra-high (2).

---

## E4. Add a rate-constant type

The rating rules are Judge 1's (`E4-E5.md`, "Rating rules"): 4 for only the type-dispatch edits
with no hidden coupling (a table indexed by the enum, a bound tied to one enumerator, a second
evaluation path, or a Python list mirroring enum integers or a struct layout); 3 with at least
one; no configuration in any entry can declare a callback rate, so none gets 5.

| Entry | Ease | Edit sites (Python + C + tests) | Extension point | Main hazard |
|---|---:|---|---|---|
| terra-max | 3 | 5 (1 + 3 + 1) | `ChemRateKind` and its `switch` (`include/chem.h:15-19`; `src/rates.c:61-104`), named fields of `ChemRate` (`include/chem.h:43-62`), a branch in the generator's `fields` (`tools/generate.py:69-133`) | the new kind must come before `CHEM_CALLBACK`, which bounds every range check and sizes the rate index (`include/chem.h:83`; `src/rates.c:16-19,140`) |
| terra-xhigh | 3 | 6 (2 + 3 + 1) | `ChemRateType` and the `switch` in `chem_rate_constants_with_callbacks` (`include/chem.h:15-19`; `src/chem.c:94-157`), a reader entry and a generator branch (`tools/mechanism.py:54-55`; `tools/generate.py:84-140`) | parameters are positional in `p[22]` (`include/chem.h:24`): the generator's list order must match the C case's reads; the rate store maps every type to a group (`src/rate_store.c:7-23`; `include/rate_store.h:6-10`) |
| terra-high | 3 | 7 (2 + 4 + 1) | `ChemRateType` and `chem_rate` (`include/rates.h:5-15`; `src/rates.c:15-67`), a generator branch (`tools/generate.py:149-194`) | positional `p[21]` in Python and C (`DESIGN.md:69-70`); the rate store's group table and count (`src/rate_store.c:8-24`; `include/rate_store.h:17`); the reference-order table sorts an unlisted type last (`tools/generate.py:71-76`) |

The new type, k = A·exp(C/T)·(T/300)^B·(1 + D·P/101325), is MICM's Arrhenius with D = 300 and
E = D/101325, so a reader-only alias would serve in all three, as in the nine; it is not counted.

**terra-max: ease 3.** Python: one branch in `fields`, a `mapped` entry and the fields with the
order-dependent unit rule the generator already applies (`tools/generate.py:71-92`); the reader
needs nothing (`tools/config.py:139-170`). C: the enum value, a case in `base_rate_value`
(`src/rates.c:61-104`) and its defaults in `chem_rate_defaults` (`:117-136`); the fields a to e
hold A to D (`include/chem.h:46`). An unhandled kind evaluates to NaN, which `chem_update_rates`
reports as non-finite (`src/rates.c:103,190`). The generator writes kinds and fields by name
(`.kind = CHEM_ARRHENIUS`, `.a = …`; `tools/generate.py:78,86-88`), so nothing mirrors an enum
integer or a layout. The coupling: `ChemRateIndex` holds `offsets[CHEM_CALLBACK+2]` and every
range check is `kind > CHEM_CALLBACK` (`include/chem.h:83`; `src/rates.c:16-19,140`), so a kind
appended after `CHEM_CALLBACK` is rejected as invalid; DESIGN.md says "add a `ChemRateKind`,
implement its calculation and validation in `src/rates.c`, and teach the normalizer and generator
to emit its data" (`DESIGN.md:103-105`), without the position. `CHEM_CALLBACK` and
`chem_override_rate_function` serve C callers (`include/chem.h:21,250`), not configurations.

**terra-xhigh: ease 3.** Python: the reader's list of allowed types (`tools/mechanism.py:54-55`)
and a generator branch with its unit rule (`tools/generate.py:84-140`). C: the enum value, a case
in the switch (`src/chem.c:94-157`; `default: return 1` makes an unhandled type an error), and a
group in the rate store (`src/rate_store.c:7-23`; `include/rate_store.h:6-10`). The couplings:
the parameters are a positional `double p[22]; /* interpreted by the rate type */`
(`include/chem.h:24`), written as a list by the generator (`tools/generate.py:85`) and read by
position in C (`src/chem.c:95-96`), and the rate store's second enum of groups. DESIGN.md: "an
enum member, parameter interpretation in the rate calculation, generator mapping, and tests"
(`DESIGN.md:48-50`), which omits the rate store.

**terra-high: ease 3.** Python: a generator branch with its unit rule, from each reaction's
`rate coefficient units` (`tools/generate.py:143-194`), and the reference-order table
(`:71-76`). C: the enum value, a case in `chem_rate` (`src/rates.c:15-67`), a group in the rate
store (`src/rate_store.c:8-24`) with its count (`include/rate_store.h:17`), and the custom-index
check if the type reads custom values (`src/cpu_solver.c:19-20`). The couplings: `p` "holds
parameters in the order used by `src/rates.c`" (`DESIGN.md:69-70`; `include/rates.h:13`), the
store's table and count, and the order table, which sorts a type it does not list last
(`tools/generate.py:76`). DESIGN.md's recipe (`DESIGN.md:69-74`) omits the store and the table.

Order: terra-max = terra-xhigh = terra-high (3); terra-max has the fewest sites and names its
fields.

---

## E5. Use another mechanism

### E5.0 What was run

- **Copies**, in Judge 1's scratch folder for this judgment (`<scratch>/judging-terra/e5-scratch/`):
  terra-max's `tools/config.py` and `tools/generate.py`; terra-xhigh's `tools/mechanism.py` and
  `tools/generate.py`; terra-high's `tools/config.py` and `tools/generate.py`; `cmp` showed each
  identical to the entry's file.
- **Interpreter**: `/usr/bin/python3`, Xcode's Python 3.9.6, with no PyYAML; Judge 1's stand-in
  `yaml.py` (`exercises/e5-check/yaml.py`), copied next to each script; only JSON inputs. Each
  run was `python3 -I run_copy.py TOOLS SCRIPT ARGS`, which puts only the copied folder on
  `sys.path` (isolated mode drops the script's folder in 3.9).
- **Inputs**: Judge 1's (`exercises/e5-inputs/`): MICM's `rate_constants_user_defined` and
  `robertson`, byte-identical to MICM's files, `permute` (species A-D, Arrhenius, a one-reactant
  Troe and a two-reactant ternary, in the default molecule cm-3 units, whose Markowitz order is
  not the input order), and `permute_m_<entry>` (the same with M first, marked by the entry's
  third-body key, plus B + C + M → D + M, Arrhenius A = 1e-30): terra-max's key is
  `parameterization: air density`, terra-xhigh's and terra-high's `third body: true`.
- **Commands**: terra-max `generate.py IN OUT --symbol S` with and without `--reorder`;
  terra-xhigh `generate.py IN OUT --symbol S` (it always reorders); terra-high
  `generate.py [--reorder | --reference-order] IN OUT`. All 24 runs exited 0.
- **Checker**: `terra_terms_check.py` maps every emitted species index to the name in the emitted
  species table and compares the reactions, as multisets of named reactants and products with
  yields, with the configuration's.
- **Results**: all 24 outputs match their configuration, in input order and reordered. The
  emitted coefficients for `permute` and `permute_m` equal the six Sol and Astra entries' in
  Judge 1's table: A + B → C, A = 1e-12, gives 602214.076; B + C + M → D + M, A = 1e-30, gives
  362661.79; the Troe gives k0_A 6.02214e-3 and kinf_A 5e-4; the ternary k0_A 602.214 and kinf_A
  1e-12. All three leave user-defined scaling factors at 1 for both tutorials.

### E5 summary

The rating rules are Judge 1's: 4 for regenerating under the benchmark's symbol (or repointing one
Makefile rule) through the generic path, with the generator tested on other configurations by
species name and no wrong-result defect; 3 with a hidden coupling that decides correctness, or a
unit convention that differs from MICM's tutorials; 2 if the reader must be fixed first.

| Entry | Ease | Hand edits | Extension point | Main hazard |
|---|---:|---|---|---|
| terra-max | 4 | 1 | `tools/generate.py CONFIG build/ts1.c --reorder` under the default symbol `chem_generated_mechanism`, or the rule at `Makefile:50-51` | none found; M must carry `parameterization: air density` (`tools/config.py:74-76,128-130`) |
| terra-xhigh | 4 | 1 | `tools/generate.py CONFIG build/ts1.c` under the default symbol `generated_mechanism`, or `Makefile:36-37` | none found; the ordering is always applied (`tools/generate.py:141`) |
| terra-high | 4 | 1 | the `build/ts1.c` rule with `TS1_ORDER_FLAGS` (`Makefile:4,171-172`); the symbols are fixed (`tools/generate.py:211-234`) | none found; the reader accepts exactly one MECHANISM entry (`tools/config.py:34-35`) |

**Common to the three.** Each benchmark maps CSV columns by name and rejects unknown ones
(terra-max `src/ts1_bench.c:79-94`; terra-xhigh `src/ts1_bench.c:79-93`; terra-high
`bench/ts1_bench.c:58-83`), so the driver is mechanism-generic. None generates kernels: the
generators write data tables, and the C library builds the pattern and plan at run time, so no
Python copy of C logic decides correctness. The species ordering runs in Python, and the
generated reactions keep their species through it (E5.0). A cells file for a new mechanism must
be written by hand.

**terra-max: ease 4.** The tests generate the two tutorial configurations, Robertson, a
parameterized and a constrained one, ordered and unordered (`Makefile:149-153` and the rules
beside it), and compare generated with hand-built mechanisms and ordered with original forms by
species name (`tests/test_configurations.c:202`; `DECISIONS.md:36-40`). The unit rule counts
reactants, third bodies and parameterized reactants (`tools/generate.py:83-118`); TS1's
configuration sets `rate units: mol m-3` (`DECISIONS.md:23-27`).

**terra-xhigh: ease 4.** The tests generate the two tutorials, Robertson and a parameterized
configuration (`Makefile:187-206`), compare the tutorial mechanism with a hand-built one
(`tests/test_tutorial.c:68-86`) and check the unit conversion (`tests/test_mechanism.py:26`).

**terra-high: ease 4.** The tests generate Robertson and both tutorials from MICM's own tree
(`Makefile:73-74,147-150`) and several fixtures (`tests/callback`, `tests/parameterized`,
`tests/process_fixture`, `tests/rate_extensions`), with both orderings "covered by
generated-mechanism tests" (`DESIGN.md:58-60`). Hazards: one generated mechanism per program,
since the symbols are fixed; a configuration with more than one MECHANISM entry is rejected
(`tools/config.py:34-35`), loudly.

Order: terra-max = terra-xhigh = terra-high (4).

## Defects found

None that gives wrong results. Hazards that fail loudly: terra-max rejects a rate kind appended
after `CHEM_CALLBACK`; terra-high's reader rejects a configuration with several MECHANISM entries.
A hazard that fails silently, in a test: terra-high's Oregonator loop would run a sixth parameter
set appended to its array as backward Euler (`tests/test_nonlinear.c:158-170`).
