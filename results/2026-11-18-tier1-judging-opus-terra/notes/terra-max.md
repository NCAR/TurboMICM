# terra-max (9071941): Judge 1's reading notes, 2026-11-18

Root: ~/NSF-NCAR/tier1-work/verify/terra-max/clang/terra-max/ (read-only).

## Counts (tools/metrics_terra.py, Judge 1's rule)
- C library and benchmark 3,166 (lib 2,932 in 19 files + bench 234); 169 functions, median 8;
  over 60: 5, over 100: 2 (chem_make_symbolic 110, src/linear.c:56; chem_solve 104,
  src/rosenbrock.c:140); comments 2.6 %; tests 5,755 C (33 files) + 128 Python; tooling 522;
  generated build/ts1.c 772 lines (122,561 bytes; data tables only).

## Documents
- README.md 71: build, benchmark, mechanisms. DESIGN.md 116: data flow, numerical components,
  constraints, external models, a 5-line extension paragraph (DESIGN.md:101-105), tuning controls.
  DECISIONS.md 78: format extensions, units, ordering, densities, RNG, API choices.
  RESULTS.md 40: coverage and one local measurement (9,632 vs 23,835 ms).
- DESIGN.md:101-105: "To add a Rosenbrock method, populate a new ChemRosenbrock parameter set. To
  change the controller or linear solver, install the respective callback. To add a rate type, add
  a ChemRateKind, implement its calculation and validation in src/rates.c, and teach the normalizer
  and generator to emit its data. A condition-dependent C function may instead be supplied as a
  callback rate." Matches the code; silent on keeping the new kind before CHEM_CALLBACK.
- DESIGN.md:22 "There is no mechanism-specific runtime code" — true (generator emits data only).
- DESIGN.md:115-116 SIMD width chosen by compiler — true.

## Public API (include/chem.h, 289 lines)
- Comments on most declarations (ownership, borrowing, layout: chem.h:13-14,59,81-83,95-96,122-124,
  137-145,180-182,198-203,219-226,233,238-246,251,255-256,260-261,265-270,273).
- ChemRateKind enum (chem.h:15-19), CHEM_CALLBACK last. ChemRate: one flat struct with named
  fields per family (a..e; k0_*, kinf_*, kint_*; fc, n; branch_*; diffusion, molecular_weight,
  probability; coefficients[16]; function/context; parameterizations) (chem.h:43-62).
- ChemRateIndex offsets[CHEM_CALLBACK+2] (chem.h:83) — table sized by the last enumerator.
- ChemRosenbrock: stages <= 6, a[15], c[15], m/e/gamma[6], new_function[6], tolerances, factor
  bounds, safety, rejection_factor, max_steps, constraint projection settings, and
  step_controller + step_controller_context (chem.h:97-111).
- ChemStepController(double error, double step, const ChemRosenbrock *, void *context) returns the
  step factor (chem.h:98).
- ChemLinearSystem view (CSR row_start, columns, lookup n x n, algebraic mask, values
  [entry*cells+cell]) and ChemLinearSolver {factor(system, alpha, ctx), solve(system, rhs, ctx)}
  (chem.h:122-135); chem_set_linear_solver per state, NULL restores default (chem.h:251-252);
  Doolittle and generic Mozart supplied (chem.h:253-254).
- 14 status codes with chem_error_string (chem.h:281-287).

## Integrator (src/rosenbrock.c, 243 lines)
- Five parameter-set constructors, public, returning ChemRosenbrock by value (:19-95); defaults()
  (:8-17). Dense style: several assignments per line (:12-15,23-26,54-60,69-73).
- chem_solve (:140-243, 104 lines): one general stage loop over p->stages (:180-203); stage count
  checked 1..6 and parameters validated (:142-148); reuse of the previous stage's raw forcing copied
  forward before the c-terms (:192-193), so chains of reused stages work.
- Error norm public chem_rosenbrock_normalized_error (:103-120), floored to 1e-10 with a
  comparison (`error < 1e-10 ? 1e-10 : error`, :118) — NaN passes through (comparison false), and
  chem_solve returns CHEM_NONFINITE on a non-finite norm (:216). MICM's NaN behaviour kept.
- Controller: callback or MICM's formula inline (:217-218); factor checked finite and > 0 (:219);
  accept if error < 1 or h < h_min (:222); clamp [h_min, h_max], no growth after a rejection
  (:226-227); tenfold cut on the third consecutive rejection overrides the factor (:232); rejected
  counted after the first acceptance (:236); Jacobian recomputed after rejection (:237).
- finish_solve clamps negatives except algebraic rows and records final_time for resumption
  (:130-138).
- One Rosenbrock path; backward Euler separate (src/backward_euler.c, 79 lines).

## Benchmark (src/ts1_bench.c, 234 lines)
- Column names mapped to temperature/pressure/species/parameter by name, unknown and duplicate
  columns rejected (:79-94); careful parsing (:24-60,95-123).
- One `ChemRosenbrock p = chem_rosenbrock_three_stage()` shared by all groups (:196), passed to
  chem_solve in both parallel loops (:201,211). Groups: chem_state_create, then create_like
  (:186-191).
- main 77 lines (:158-234), single cleanup label (:227-230).

## Linear solver (src/linear.c, 363 lines)
- chem_linear_factor / chem_linear_solve dispatch to the installed ChemLinearSolver or the default
  plan (:239-251); integrator calls each once per attempt/stage (rosenbrock.c:177,200).
- chem_factor (:167-199) replays the precomputed plan; checks every pivot per cell for zero or
  non-finite (:177-178) and returns CHEM_SINGULAR.
- chem_make_symbolic 110 lines (:56-165).

## Rates (src/rates.c, 368 lines)
- base_rate_value switch (:61-104) returns NAN for an unhandled kind; chem_update_rates returns
  CHEM_NONFINITE on a non-finite rate (:190).
- chem_rate_defaults switch (:117-136); validate_rate bounds kind by CHEM_CALLBACK (:140);
  chem_rate_index_init counts[CHEM_CALLBACK+1] and bound (:16,19).
- Constants: gas_constant from k_B * N_A (:8).

## Generator and reader (tools/generate.py 202, tools/config.py 217)
- Kinds emitted by name `.kind = CHEM_<NAME>` (generate.py:78) — no integer mirror.
- Unit rule by reaction order (reactants + air-density power + parameterized reactants;
  generate.py:83-84): A^(order-1) for Arrhenius/Taylor/tunneling/branched X; Troe k0^order,
  kinf^(order-1); ternary k0^(order-1), kinf^(order-2) (:86-118). User rates unscaled (molar,
  DECISIONS.md:23-27). `rate units: mol m-3` disables conversion (config.py:40-49,157).
- Unknown reaction type raises (generate.py:76-77); unknown camp-data type raises (config.py:94-95);
  unknown species in a reaction raises (config.py:116-117); third body = species with
  `parameterization: air density` (config.py:74-76,129-130).
- reorder_mechanism in Python only (generate.py:12-57); the C library has no reordering of its own
  on the solver path; indices remapped for reactions and constraints together (:46-56).

## Build
- CFLAGS -std=c11 -Wall -Wextra -Wpedantic (Makefile:25); no warnings with either compiler
  (verification). build/ts1.c rule (Makefile:50-51). test-cpu builds 33 programs then
  tests/run_cases.py (Makefile:194-195).
