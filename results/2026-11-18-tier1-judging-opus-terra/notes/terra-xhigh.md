# terra-xhigh (288f274): Judge 1's reading notes, 2026-11-18

Root: ~/NSF-NCAR/tier1-work/verify/terra-xhigh/clang/terra-xhigh/ (read-only).

## Counts (tools/metrics_terra.py)
- C library and benchmark 2,865 (lib 2,658 in 20 files + bench 207); 172 functions, median 7;
  over 60: 6, over 100: 0 (largest chem_rosenbrock_solve_group 88, src/rosenbrock.c:637;
  chem_rate_constants_with_callbacks 83, src/chem.c:84; chem_rosenbrock_solve 75,
  src/rosenbrock.c:561); comments 2.3 %; tests 4,209 C (32 files) + 151 Python; tooling 366;
  generated build/ts1.c 3,703 lines (157,604 bytes, data only, hex floats).
- 101 library and benchmark lines over 120 characters, longest 208 (src/rosenbrock.c 20,
  src/solver.c 18, src/matrix.c 10, src/ts1_bench.c 9, src/sparse.c 9).

## Documents
- README.md 23, DESIGN.md 95, DECISIONS.md 50, RESULTS.md 21.
- DESIGN.md:77-88 per-change guidance, candid and matching the code: a new set is a constructor
  in src/rosenbrock.c consumed by both loops; step control: "change the error-to-step calculation
  and the accept/reject logic in both chem_rosenbrock_solve and chem_rosenbrock_solve_group"; the
  factorization and substitution routines "are the boundary for replacing the linear solver"; "The
  current implementation selects the two built-in LU variants by enum, rather than through a solver
  callback"; DESIGN.md:46-48 "the interface should be made pluggable before swapping it"; new rate
  type: "an enum member, parameter interpretation in the rate calculation, generator mapping, and
  tests" (:48-50).
- DECISIONS.md:16-20 unit rule by reaction order including the falloff units.

## Public API (include/chem.h, 238 lines)
- ChemReaction: type, term ranges, parameter index, `double p[22]; /* interpreted by the rate
  type */` (chem.h:20-25). ChemRateType enum, CHEM_CALLBACK last (:15-19).
- ChemRosenbrockParameters: stages <= 6, new_function[6], a/c[15], m/e[6], scalar gamma, order,
  factor bounds, safety, rejection factor, h limits, tolerances, constraint-init settings,
  `ChemLUAlgorithm linear_solver` enum (:65-77). No controller or linear-solver hook.
- ChemSolverStatus enum 0..6 with strings (:83-88; src/rosenbrock.c:9-21), but the integrators
  return bare integers: `return 1/2/3/4/5` (src/rosenbrock.c:565-588,616,629; 63 `return N;`).
- Two Rosenbrock integrators: chem_rosenbrock_solve (single cell, :561-635) and
  chem_rosenbrock_solve_group (shared step over cells, :637-724), each with the stage loop, the
  controller (:617, :706) and the accept/reject rules (:619-630, :708-718). Backward Euler single
  and grouped (:726-...).
- chem_linear_factor / _in_place / _solve expose the integrators' static factor and substitute
  (src/rosenbrock.c:287-347,365-402), so the LU tests reach the solver's own kernels.
- chem_diagonal_markowitz_order in C (src/rosenbrock.c:412-434): same rule as the generator's
  Python reorder_species (tools/generate.py:17-45), two copies; TS1 uses the Python one (:141).

## Integrator details
- General stage loop over p->stages (:590-606), stage reuse copied forward before the c-terms
  (:602). Error norm public, cell-major (:116-131); returns NAN for a non-finite norm, and the
  solver returns 4 (:614-616): MICM's NaN behaviour kept.
- Controller inline (:617) with MICM's rules (:619-630); a rejection counts a Jacobian update
  without recomputing it, since factor() works on a copy (:292-300,628).
- Grouped integrator mallocs a buffer each call (:649).
- factor(): copy −J with the shift into w->matrix (:292-300), then Doolittle by pattern tests
  (:301-315) or Mozart by precomputed lists (:316-330); singular or non-finite pivot returns 1.

## Rates, forcing, Jacobian (src/chem.c, 264 lines)
- chem_rate_constants_with_callbacks switch on type with `default: return 1` (:94-157); p[]
  positional per type (:95-153). Third bodies and parameterized species folded into k (:158-161).
- One forcing and one Jacobian implementation (chem_forcing :168-188; chem_jacobian_add
  :196-221), used by the integrators and wrapped by the process API (:223-...).

## Benchmark (src/ts1_bench.c, 207 lines)
- Columns mapped by name, unknown column rejected (:79-93); defaults from the generated
  configuration (:99-101); rate constants at read time (:117); one workspace per thread
  (:176-179); `p = chem_rosenbrock_three_stage()` (:170), read-only, shared; argument errors exit
  1 without a message (:149-160).

## Generator and reader (tools/generate.py 189, tools/mechanism.py 83)
- Reader validates camp-data types, species references, qty, units, environment (mechanism.py:
  31-83); allowed reaction types listed (:54-55).
- Generator: positional parameter lists per kind (generate.py:84-118) emitted into p[] (:169);
  kinds by name `CHEM_%s` (:169); unit rule by reaction order counting third bodies and
  parameterized reactants (:121-139); ordering always applied (:141); numbers as hex floats (:10-14).

## Tests
- tools/test_cpu.py runner (66 lines): per-id `<id> 1|0` lines, duplicates and unknown ids raise.
- tests/test_mechanism.py: tutorials, unit conversion, TS1 counts and defaults, parameterized
  species, bad references (:17-74).
- Sample: 9 of 12 faithful.
