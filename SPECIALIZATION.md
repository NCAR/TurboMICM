# TS1 specialization

The generator fixes the reaction order, species indexes, stoichiometric values,
Jacobian dependencies, symbolic fill, and natural-order factorization schedule.
Generated forcing and negative-Jacobian code has no runtime reaction or
stoichiometry traversal. Generated LU and triangular solves have constant row
and column indexes. The specialized binary uses these functions through
`chem_specialized_ops`.

Temperature, pressure, air density, concentrations, rate parameters, duration,
and tolerances remain runtime data. Rate constants are computed by the shared
generic evaluator and cached per solver workspace when their inputs are
unchanged. The three-stage Rosenbrock method, error estimate, controller, and
step-size rules are unchanged. The generic Tier 1 path has no specialization
callbacks and remains the comparison implementation.

## Measured effects

The measurements below use Clang, `-O3 -DNDEBUG`, the 100-cell seed-1 input,
100 steps of 30 seconds, and one thread. Times are elapsed milliseconds from the
solver's timed region.

| Generated operations enabled | Time | Relative to generic |
| --- | ---: | ---: |
| Generic Tier 1 path | 2,861 | 1.00× |
| Straight-line forcing and Jacobian updates | 1,965 | 1.46× |
| Plus generated sparse LU and triangular solves | 843 | 3.39× |

The forcing/Jacobian result was measured before adding the LU schedule; the
final row measures all retained specializations together. No attempted
specialization was dropped. Solver counts were unchanged across these runs.

For the full 10,000-cell, 100-step, 30-second run with seed 2026, the final
specialized executable took 7.76 s on 18 threads and 90.42 s on one thread.
The generic implementation took 21.34 s and 305.34 s, respectively. MICM took
12.74 s and 134.18 s. The full-run counters matched for all three programs.

## Trade-offs

The generated TS1 C file is about 1.12 MB. Generation is mechanical and follows
the normalized mechanism records, so adding a reaction or changing a fixed
coefficient requires editing the mechanism input and regenerating. Runtime
inputs such as user-defined rates and photolysis values are never embedded in
the generated code. See `DESIGN.md` for the emitted schedules and their
relationship to the generic solver.
