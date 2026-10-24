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
