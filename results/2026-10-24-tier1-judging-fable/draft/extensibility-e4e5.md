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
