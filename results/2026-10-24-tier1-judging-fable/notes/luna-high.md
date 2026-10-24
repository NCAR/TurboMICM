# luna-high: Judge 2's reading notes

Final commit d3c2480. Worktree ~/NSF-NCAR/tier1-work/verify/luna-high/clang/luna-high.

## Shape
- `include/chemistry.h` (727 lines): matrix and CSR views, species/phase/system records, constraints,
  dense LU, a `chem_reaction` with `p[24]` positional parameters (`:43-53`), and two solver APIs.
- Two Rosenbrock integrators: the generic single-cell `rosenbrock_solve_impl` (`src/solver.c:489-690`,
  any set to 6 stages, dense `n x n` Jacobian and matrix with sparsity lists, natural-order LU
  `factor`/`solve_factored`, `:305-344`) and the TS1 path `chem_rosenbrock3_solve_group_with_workspace`
  (`:767-945`, 179 lines) with the three-stage coefficients inlined (`:796-804`), the controller inlined
  (`:908`) and a per-cell scalar loop that recomputes rates if conditions changed, builds the dense
  Jacobian, factors and solves each cell separately (`:812-888`); the group shares only the step size and
  the error norm (`:890-907`).
- No sparse storage on the solve path: `w->jac` and `w->a` are `n*n` dense (`:90-91`), 210 x 210 doubles per
  workspace, with `row_cols` lists to skip structural zeros (`:254-303`).
- Species order: `tools/generate_c.py` reorders the species table by its own Markowitz pass (`:40-86`), but
  emits each term's index as `term.get("index", names.get(term["name"]))` (`:126`); `tools/read_config.py`
  writes `"index"` as the input position (`:83-86`). With a CAMP input the reactions therefore name the
  input positions in a permuted table. `mechanisms/ts1.json` comes from `tools/import_ts1.py`, and TS1
  passes.
- Rates: `chem_rate_constants` evaluates SURFACE, ARRHENIUS, TAYLOR and LAMBDA itself and defers the rest to
  `chem_reaction_rate`, which also has ARRHENIUS and TAYLOR cases and derives pressure from the air
  density (`src/mechanism.c:140-225,243-295`); PHOTOLYSIS, USER_DEFINED and SURFACE share one case (`:149-153`).
- Error norm floors with `fmax(..., 1e-10)` (`src/solver.c:655,907`): a NaN step is accepted.
- Benchmark `src/main.c`: `--threads` only; one workspace per thread with the cells split evenly (`:291-307`);
  cell-major state; the warm-up and the timed loop call the group solver per step (`:308-348`).

## Readability
- Dense: several statements per line (`src/solver.c:51-58,296,332`), 2-space indent; common count 0.4 %
  comments, median 10 lines, 10 over 60 and 5 over 100.
- Two integrators with the same rules written twice (`:558-680` and `:805-938`); literal 0.2, 6.0, 0.9 in
  both controllers (`:656,908`).
- `chemistry.h` mixes five concerns; names are long and descriptive.
- Documents: README 103, DESIGN 124 (prose; the extension paragraph says the factorization "can be replaced
  behind a linear-solver interface", which does not exist), DECISIONS 120, RESULTS 126.

## Maintainability
- Tests: `tools/test_cases.py` holds a literal set of 510 ids (`implemented`, `:7-639`), runs nothing itself,
  and prints `PASS` with a per-prefix reason for every id in the set (`:640-769`); the checks run earlier in
  `make test-cpu` as one program, `bin/test_core` (`Makefile` test-cpu), whose `main` runs the checks in
  sequence and returns a code at the first failure (`tests/test_core.c:4801-...`), then `tests/test_codegen.py`.
  So one program's exit status stands behind 662 PASS lines, and the mapping from an id to a check exists
  only in the agent's mind. 18 N/A.
- `tests/test_codegen.py` checks the generator's ordering, defaults and labels, and reads MICM's two
  rate-constant tutorials through `read_config.py` (`:94-129`); it compares the Troe and ternary unit
  conversion with fixed exponents (`:120-126`).
- Makefile: `CC ?= cc`; the benchmark and the test program each compile all eleven sources in one line
  (no objects); `clean` removes `generated/`.
- Logic in two places: two Rosenbrock integrators, two rate evaluators, pressure derived from density in
  one and passed in the other.

## Exercises
- E1: the generic API takes a public `chem_rosenbrock_parameters` (`include/chemistry.h:608-628`); the TS1
  path hard-codes three stages and its coefficients (`src/solver.c:796-804,859-887`). Ease 2.
- E2: no hook; the controller formula is inlined in both integrators (`:656,908`). Ease 2.
- E3: static `factor`/`solve_factored` on dense arrays, called directly in both integrators (`:591,624,854,877`);
  no interface. Ease 2.
- E4: C: enum (`include/chemistry.h:12-26`), cases in `chem_reaction_rate` and, for pressure-dependent
  types, `chem_rate_constants` (`src/mechanism.c:140-225,252-283`), slots; Python: `KINDS`, `TYPE_ORDER`,
  the per-kind slot code (`tools/generate_c.py:8-37,148-198`), `read_config.py` units; tests in `test_core.c`
  and `test_codegen.py`. Ease 3.
- E5: `read_config.py` then `generate_c.py`; one fixed symbol `chem_generated_mechanism`; the index defect
  above gives wrong reactions for any CAMP input whose Markowitz order differs from the input order;
  fixed-exponent units. Ease 1 to 2, pending the E5 helper's run.
