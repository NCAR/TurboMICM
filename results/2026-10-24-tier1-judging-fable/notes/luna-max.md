# luna-max: Judge 2's reading notes

Final commit da2d722. Worktree ~/NSF-NCAR/tier1-work/verify/luna-max/clang/luna-max.

## Shape
- Table-driven `chem_mechanism` (flattened reactant/product lists, `chem_reaction.p[16]` positional
  parameters, an int `kind`, `include/kinetics.h:71-109`); the generated TS1 is 950 lines of tables
  (`build/ts1_model.c`) in the species order chosen by a Python Markowitz pass (`tools/generate_c.py:165-192`)
  whose cost rule, `max(0, r-1) * max(0, c-1)` with a degree-sum tie-break, is not MICM's; C does no ordering
  (`chem_make_pattern` takes the table's order, `src/solver.c:93-190`).
- One `chem_solver` for all cells, created by the benchmark with the rate constants of every cell cached at
  creation (`src/solver.c:799-855`, `src/benchmark.c:325-326`); `chem_solver_step` runs the groups with an
  `omp parallel for` inside the library (`src/solver.c:979-996`), one group per thread by default.
- Integrator `solve_group` (`src/solver.c:644-797`, 154 lines): cell-major `y[cell*n + v]`; inside each
  attempted step a loop over the group's cells does, per cell, the forcing and sparse negative Jacobian
  (`evaluate_rhs`, `:247-324`, scalar), a copying sparse LU (`factor`, `:326-356`, scalar with `fma`),
  the stage solves, and the candidate; the error norm is then taken over the group (`:752-755`). No SIMD
  across cells anywhere on the TS1 path.
- Controller `p->controller(error, order, safety, ctx)` returns the factor, clamped by the solver
  (`:761-770`); factor/solve callbacks `(pattern, alpha, negative_jacobian, lu, ctx)` per cell in the
  parameters struct (`include/kinetics.h:172-181,207-210`).
- Error norm floors with `fmax(sqrt(sum/count), 1e-10)` (`src/solver.c:80`): a NaN sum becomes 1e-10, the
  finiteness test at `:767` passes, and the step is accepted with NaN values (`:772-776`).
- Backward Euler (`src/backward_euler.c:35-117`): one Newton iteration loop over the whole interval with a
  dense LU, no step reduction sequence; on non-convergence it restores the state and returns -1.
- Reader: `tools/generate_c.py:43-151`; units converted for TROE (`k0_A` by the square of the factor,
  `kinf_A` by the factor) and TERNARY (`kinf_A` divided) only, with fixed exponents (`:88-94`); Arrhenius
  and the others are not converted.

## Readability
- One statement per line, 2-space indent, C89-style declarations at block tops; common count 0.3 % of
  lines with a comment (the lowest of the nine), median function 13 lines, 10 over 60 and 3 over 100.
- Positional `p[16]`: `r->p[1] * exp(-r->p[0] / temperature)` for branched (`src/kinetics.c:70`), slot
  meanings only in the Python `KINDS` table (`tools/generate_c.py:16-30`).
- Status: most functions return -1 for every failure (`src/solver.c:945-1001`); `chem_solver_state` exists
  but the solver returns ints.
- Documents: README 91, DESIGN 101 with a four-item extension list that matches the code, DECISIONS 45,
  RESULTS 33.

## Maintainability
- Tests: `bin/test_cpu` reads MICM's case list and dispatches each id by suite name to a handler
  (`tests/test_cpu.c:4330-4357`, `run_case`); many ids inside a suite share one function (for example the
  LU suites). `tests/run_cases.py` then overwrites some rows from other programs: the FlowTube regression
  row from the 100-cell TS1 gate, two example rows, six tutorial rows (`tests/run_cases.py:129-144`).
- The generator runs on the three examples in the build (`Makefile` tutorial_*_model.c) and
  `verify_configurations` checks counts and the presence of the C table (`tests/run_cases.py:23-41`).
- Makefile: `CC = clang` unconditionally (a `CC=gcc-15` on the command line overrides it), objects in one
  `build/` for both compilers (README: run `make clean` between compilers), each test program linked from
  objects with `-isysroot`.
- The generated symbol is always `chem_model` (`tools/generate_c.py:326`), so two mechanisms cannot be
  linked into one program; `make mechanism CONFIG=...` builds a separate archive.
- Error handling: `chem_solver_step` allocates per call (`src/solver.c:972-973`); `evaluate_rhs` mallocs a
  gradient per cell per evaluation when constraints are present (`:292`).

## Exercises
- E1: fill `chem_rosenbrock_parameters` (public, `include/kinetics.h:183-211`) and pass it to
  `chem_solver_create`; constructors modify a default; generic loop to 6 stages (`src/solver.c:666-670,684-750`);
  DESIGN says how. Ease 5.
- E2: `chem_step_controller(error, order, safety, ctx)` returns the factor (`include/kinetics.h:172-173`,
  `src/solver.c:761-765`); no h; the context lives in the one parameters struct of the one solver shared by
  every group and thread (`src/solver.c:27,644-649`). Ease 3.
- E3: `factor`/`solve` callbacks per cell with the public `chem_sparse_pattern` (row and column offsets,
  positions) (`include/kinetics.h:160-181`, `src/solver.c:657-658,705-708,737`); documented. Ease 5.
- E4: C: enum before the sentinel (`include/kinetics.h:8-21`), a `chem_rate_constant` case
  (`src/kinetics.c:30-112`), the range check `r->kind > CHEM_CALLBACK` (`src/solver.c:615`), slots;
  Python: a `KINDS` entry and the unit branch (`tools/generate_c.py:16-30,88-94`); a `tests/test_cpu.c`
  case. `CHEM_CALLBACK` exists. Ease 3.
- E5: `make mechanism CONFIG=path` or `tools/generate_c.py config out.c`; one `MECHANISM` entry required
  (`tools/generate_c.py:71-72`); the symbol is fixed, so the benchmark runs another mechanism only by
  replacing `mechanisms/ts1.json`; ordering in Python only (consistent with the table); fixed-exponent
  units (wrong for a bimolecular Arrhenius in molecule cm-3 units); tests generate the three examples.
  Ease 3.
