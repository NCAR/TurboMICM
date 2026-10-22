# Mechanism specialization

The Python generator reads the same JSON/YAML configuration as the Tier 1
library. There is no TS1-specific logic in the generator. The retained TS1 input
contains the chemistry; adding or changing a reaction requires only editing that
input and rebuilding.

The reader resolves included files, expands repeated reactants and branched
reactions, and converts rate units. The symbolic planner determines variable
species, reference Markowitz order, Jacobian entries, sparse fill, LU
elimination operations and triangular dependencies. The emitter fixes these
indices, stoichiometric coefficients and rate parameters in portable C11. All
reaction, Jacobian and sparse operation lists become straight-line arithmetic
inside independent-cell OpenMP SIMD loops. Triangular solves accumulate a row
in a local variable. The three-stage method fixes its stages and coefficients;
stage vector operations are fused while retaining each scalar's original order.

Temperature, pressure and air density are live state fields. Every named user,
photolysis, emission, loss and surface parameter is read from mutable state.
Concentrations, requested intervals and relative/per-species absolute tolerances
also stay live. Prepared constants in rate formulas use C double expressions;
Python does not evaluate transcendental constants for the target. Rates are
refreshed on demand and at each scenario step. Ordinary benchmark rates are
computed before its warm-up, matching Tier 1.

The library still constructs a generic model for metadata, initialization and
CSV name lookup outside the timed region. Its numerical operation tables are
not interpreted by the specialized solve. Public concentrations and rates keep
the generic species-major layout. The private tiled workspace stores complete
small cell tiles contiguously; it is reusable and owned by the state. Cell tiles
share the parent group's step controller. They do not become independent solver
groups. Dummy lanes have zero rates and are excluded from the error reduction.
The reduction visits real cells in exactly the original species-major order.

## Experiments

All runs use the default compiler flags; no flags were tuned. Measurements here
are local observations, not the evaluator's ranking. `results/tier2/` retains
compact evidence. The 10,000-cell input uses seed 2026, 100 steps of 30 seconds,
and a separate untimed warm-up. Group size is 32 unless stated otherwise.

The initial direct unrolling retained the existing layout and emitted one large
function per kernel. At 18 threads its median was 3,156.6 ms (three runs), versus
3,001.6 ms for Tier 1. One-thread times were 41,750.3 and 38,305.0 ms respectively.
Unrolling alone was therefore not retained as the default strategy.

Splitting generated loops into functions of approximately 256 source lines cut
a rebuild to 45.36 seconds, but its 18-thread median was still 3,189.3 ms. A
10,000-repeat, 32-cell kernel benchmark measured LU factorization at 555.8 ms
versus 490.7 ms for the generic plan, and triangular solves at 142.9 versus
122.7 ms. This motivated changing layout and compiling constant-width tiles,
instead of merely expanding more instructions. Function splitting is retained
to bound compiler work and remains configurable.

A width-8 private tiled layout with direct calls took 307.57 seconds to rebuild
and 3,717.8 ms in one 18-thread run. Allowing the compiler to merge fixed-width
parts was counterproductive. Calling fixed kernel parts through private volatile
function pointers bounds their size using ISO C, with no noinline attributes or
compiler flags. This rebuilt in 48.72 seconds and measured 3,053.9 ms median
(three runs). All five disclosed scenarios passed at 14.05--14.37 digits overall.
Both direct calls and flat layout remain available as experimental options.

A symbolic-only ordering comparison found 6,986 elimination updates and 2,381
filled entries in reference order. Greedy product, minimum-fill and degree
heuristics gave respectively 6,619/2,358, 6,342/2,300 and 6,881/2,375. These modest
reductions did not motivate changing the reference order; no runtime speedup is
claimed for these unimplemented orderings.

The subsequent 18-thread results are medians of three complete 100-step runs:

| Change from the preceding retained variant | Elapsed ms | Decision |
| --- | ---: | --- |
| Four-cell tiles, divisions, one chemistry job per loop | 2,522.9 | Retain four-cell tiles |
| Store reciprocal pivots in packed factors | 2,489.8 | Retain; day-night gates pass |
| Batch eight chemistry/Jacobian jobs per loop | 2,389.4 | Retain |
| Dynamic scheduling, otherwise unchanged | 2,471.3 | Keep static scheduling |
| Doolittle accumulation, static scheduling | 1,993.5 | Retain |
| Two-cell tiles with Doolittle | 2,417.2 | Keep four-cell tiles |
| Sixteen-cell tiles with Doolittle | 2,673.6 | Keep four-cell tiles |

Doolittle uses the same filled structure and evaluates each entry's elimination
contributions in increasing pivot order. Accumulating an entry in a local
variable reduces repeated target loads and stores. Packed factors store inverse
diagonals; all three stage solves multiply by these inverses. This remains a
direct double-precision solve of the original linear system, without extra
factorization reuse. Every attempted step still builds and factors its matrix;
rejected steps follow MICM's original Jacobian-rebuild policy.

The final defaults are chemistry batch 8, algebra batch 8, approximately 256
lines per function, tile width 4, Doolittle with reciprocal pivots, groups of 32
and static group scheduling. Runtime/compiler-dependent alternatives remain
parameters. The public inspection kernels expose the ordinary SoA Mozart factor
layout; the timed packed solve uses the chosen private factor layout. Full
trajectory comparisons exercise that private path with both full and padded
cell tiles, changing rates, conditions and tolerances between calls.
