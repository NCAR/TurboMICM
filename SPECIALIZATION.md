# Mechanism specialization and measurements

The build-time compiler is `tools/specialize.py`. It shares CAMP parsing and
unit conversion with `tools/configure.py`, then independently binds the
mechanism and emits ISO C11. Generated sources are ignored build products.
There is no mechanism-name test, TS1-specific reaction kernel, numerical table
in Python, third-party C code, or runtime Python dependency.

## What becomes code

- Stable rate-type ordering and the reference diagonal Markowitz species
  permutation; species, reaction, custom-parameter and packed matrix indices.
- Reactant occurrences and individual product yields, including repeated
  species, and the removal of each parameterized reactant from the integrated
  state. Its environmental multiplier remains in rate updates.
- Analytic negative-Jacobian contributions, symbolic sparse fill, all 1,244
  below-diagonal divisions and 6,986 LU updates for TS1, and forward/backward
  substitution. TS1 has 209 integrated species and 2,381 packed matrix entries.
- Fixed-parameter rate formulas, including all the configurable rate families.
  C/libm calculates transcendental functions. The C model constructor prepares
  the two branch constants once, retaining the reference's rounding.
- The three Rosenbrock stages and their coefficient values. Independent array
  passes are fused, while every element retains the original arithmetic order.
  The shared-forcing copy, ordered error reduction, rejected-step Jacobian
  reconstruction, factorization frequency and step controller remain MICM's.
- Bounded arithmetic batches and helper functions. Full groups can use constant
  row strides; partial groups retain their actual stride and lane count.

Temperature, pressure, independently supplied air density, concentrations, all
custom rate parameters (including aerosol sizes/numbers and every photolysis
label), external step length and per-state relative/absolute tolerances remain
runtime data. Initial values are defaults only. No small-concentration pruning,
rate freezing across scenario steps, approximate Jacobian, factorization reuse,
precision reduction, or chemistry change is used.

The generic constructor still supplies names, default initialization and
workspace ownership outside timing. Generated arithmetic requires the default
ordered model of its own emitted descriptor. The original generic library and
all its solvers remain available.

## Experiments

Development measurements use clang, `-O3 -DNDEBUG`, the specified OpenMP setup
and this laptop. These single-run exploratory measurements are diagnostic,
not ranking medians. The kernel benchmark uses one 32-cell group and 2,000
calls, with matrix/RHS copies included in the two linear algebra measurements.
Temperature and chemistry are identical between the paired kernel calls.

| Kernel, microseconds/call | Generic, final measurement | Generated, selected defaults |
|---|---:|---:|
| Rate update | 97.628 | 40.492 |
| Forcing | 15.042 | 9.006 |
| Negative Jacobian | 28.706 | 18.342 |
| Factorization plus matrix copy | 45.863 | 48.521 |
| Triangular solve plus RHS copy | 12.332 | 10.985 |

The fixed rate formulas, chemistry indices and yields reduce dispatch and
indirect addressing. Unrolled LU alone does not improve this measurement;
its memory traffic and instruction footprint offset the saved schedule lookups.
The whole-solver results in RESULTS.md include the stage-pass fusion and all
kernel effects; no isolated speedup is claimed for the stage fusion alone.

| Trial | Observed result | Choice |
|---|---|---|
| One cell loop containing a whole LU pivot | factor/copy 100.417 us versus paired generic 46.703 us; full 18-thread run 4,591 ms | Replaced; excessive expression/register pressure |
| Dynamic strides, eight algebra operations per loop | factor/copy 61.586 us | Replaced with shorter groups while investigating strides |
| Dynamic strides, one algebra operation per loop | factor/copy 48.406 us; full run 3,838 ms | Correct, useful comparison; fixed strides improve chemistry |
| Constant stride and constant 32-lane trip counts | forcing 9.502 us and Jacobian 19.189 us, but factor/copy 79.965 us; 6.63 MB executable; full run 5,375 ms | Dropped constant trip counts |
| Constant stride with runtime loop bounds, one algebra operation per loop | factor/copy 49.003 us; forcing 9.467 us; Jacobian 19.373 us | Retain bounded helper dispatch and runtime loop bounds |
| Same layout, eight algebra operations per loop | factor/copy 50.654 us, solve/copy 11.517 us; full run 2,601 ms | Selected; the three solves outweigh the small factor cost |
| Eight-cell tiles within unchanged 32-cell solver groups | factor/copy 58.941 us, solve/copy 13.724 us; full run 3,004 ms | Disabled by default; extra passes/calls outweighed cache locality |
| Row-oriented independent LU updates plus combined initial forcing/Jacobian | factor/copy 51.324 us; full run 2,737 ms, one-thread 37,722 ms | Disabled by default; no full-run benefit |
| Greedy minimum-fill and minimum-work symbolic orderings | 2,293/2,350 nonzeros and 6,314/6,601 LU updates, versus 2,381 and 6,986 | Symbolic exploration only; retain the reference permutation and its error-reduction order |

Exploratory run timings have different thermal/scheduling conditions and should
not be subtracted to assign independent causal speedups. The paired final
benchmark comparison and exact commands are in RESULTS.md. All compiler
optimization flags stayed unchanged throughout.

The generator retains optional tiling, joint forcing/Jacobian emission and
row-oriented LU traversal for other machines. Joint emission interleaves
independent Jacobian columns while preserving each entry's reaction order.
Every LU update within a pivot targets a distinct entry and reads only the
completed pivot row/column, so row and column traversals preserve each entry's
arithmetic. A separate numerical fixture tests these transformations together,
using an odd fixed width, uneven tiles and different batching/splitting values.

Defaults and controls are listed in DESIGN.md. Regenerate through Make, or run
`.conda/bin/python tools/specialize.py CONFIG OUTPUT --symbol NAME`. Changing a
reaction or rate parameter in the configuration changes the schedules/formulas
without generator edits. Tests cover TS1, all three Tier 1 configuration-file
mechanisms, every supported configurable rate, zero-order chemistry, repeated
terms, parameterized species, runtime input changes and both tolerance regimes.
