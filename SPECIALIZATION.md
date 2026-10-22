# Mechanism specialization

The compiler is `tools/specialize.py`. It consumes the same normalized MICM
configuration as the generic constructor generator. Its symbolic `Plan`, rate
formula emitter, and SIMD-loop emitter are separate. TS1 is an input, not a
special case: adding a reaction or changing a fixed coefficient changes the
generated operations on the next build. The generic library remains intact.

## What becomes code

- Species order, all concentration/rate/parameter indices, reaction orders,
  stoichiometric multiplicities, product yields, and parameterized-reactant
  counts become literals. Ordered duplicate contributions are retained.
- Rate formulas have fixed coefficients and no runtime switch. Arrhenius,
  Troe, both activation formulas, tunneling, both branch channels, surface,
  user-defined, Taylor, and reversible rates all have emitters. The compiler
  computes constant branch preparation and surface species coefficients.
- The negative analytical Jacobian has fixed destinations and fill. TS1 has
  209 evolving species, 547 reactions, and 2,381 filled entries. Mozart's plan
  contains 8,439 operations. The same algorithm generates other sizes/patterns.
- LU pivot inverses, scaling, outer products, and triangular rows become
  fixed-index expressions. A second emitter implements exact Doolittle
  factorization, used in the experiments below. Neither changes the system.
- The three Rosenbrock stages call those kernels directly. Stage copies and
  vector updates are fused while retaining the reference's arithmetic order.
  The existing coefficient table, ordered error norm, acceptance/rejection
  rules, and step-size controller remain authoritative.
- Default-sized groups have constant array strides. Cells remain contiguous
  within each species/entry. The compiler selects SIMD lanes; no target-specific
  instructions, types, attributes, or optimization flags occur in the source.

Generated constructors still use the generic API for setup and metadata. This
cost is outside the specified timer. Generated kernels do not traverse reaction,
derivative, sparse-operation, or species-name metadata during integration.

## What remains data

Temperature, pressure, air density, every concentration and named custom rate
parameter, step length, relative tolerance, and species absolute tolerances
remain mutable. A custom rate is always loaded from the state, including every
photolysis, heterogeneous rate, surface radius/number, and other user-defined
parameter. Conditions and parameters are re-evaluated at every scenario step.
There is no assumption about seed, latitude, season, initial NOx, or day length.

Each state retains one adaptive history for its entire group. A tile cannot
change the error norm or split a group. Short groups use emitted runtime-stride
loops. Factors are recomputed whenever MICM recomputes them, and every solve and
intermediate value uses double precision. There is no approximate Jacobian,
approximate solve, clipping inside stages, or precomputed trajectory.

## Experiments and selection

All experiments used clang and `-O3 -DNDEBUG`; compiler flags were never tuned.
The short tuning workload was 100 supplied seed-1 cells, one thread, 100 timed
30-second steps after one warm-up. It includes three full groups and a short
group, so its speedup should not be extrapolated to the primary workload.

Initial single observations (milliseconds):

| Implementation | elapsed_ms | Outcome |
| --- | ---: | --- |
| Retained generic solver | 465.191 | Baseline |
| Fixed chemistry, LU, solves, and fused stages; 32-unit loops, runtime strides | 766.111 | Slower: excessive address and register pressure |
| Add constant strides for full groups | 479.074 | Retained the constant-stride option |
| Doolittle, eight-unit loops, constant strides | 438.323 | Explored further below |
| Four-cell tiles, 32-unit Mozart loops | 681.944 | Dropped as default |
| Four-cell tiles, 256-unit Doolittle loops | 793.870 | Dropped as default |
| Separate dispatch/kernel translation units, four-cell tiles | 671.495 | Splitting alone did not fix the regression |
| Separate units, runtime loop bound, one operation unit per loop | 400.068 | Retained bounded loops and separate compilation |

The last change retains a constant stride but passes the loop trip count to the
kernel. It prevents the caller's group size from expanding every SIMD loop into
a large fixed-trip instruction sequence. Separate translation units express
this using ordinary C, without compiler-specific inlining directives.

The subsequent controlled sweep took two samples per setting. Every setting
passed the TS1 native parity test and the MICM final-state gate. Times below are
the arithmetic mean of the two samples; raw samples, component profiles, and
incremental build times are in
[`validation/specialization-experiments.json`](validation/specialization-experiments.json).

| LU | Function block | Units per cell loop | Cell tile | Mean elapsed_ms |
| --- | ---: | ---: | ---: | ---: |
| Doolittle | 32 | 1 | 32 | 389.877 |
| Mozart | 8 | 1 | 32 | 403.635 |
| Mozart | 128 | 1 | 32 | 394.471 |
| Doolittle | 8 | 1 | 32 | 396.174 |
| Doolittle | 128 | 1 | 32 | 390.293 |
| **Mozart** | **32** | **4** | **32** | **374.818** |
| Mozart | 32 | 1 | 16 | 449.432 |
| Mozart | 32 | 8 | 32 | 417.518 |
| Mozart | 32 | 16 | 32 | 439.447 |

For the selected setting, the 2,000-call, 32-cell component profile was:

| Kernel | Generic ms | Generated ms |
| --- | ---: | ---: |
| Forcing | 31.155 | 17.280 |
| Negative Jacobian | 59.264 | 38.059 |
| Shift and LU, including input copy | 100.349 | 85.044 |
| Triangular solve, including RHS copy | 24.708 | 30.211 |

Triangular row unrolling did not improve its standalone time; the full generated
path was selected on whole-solve timing. Fixed rate formulas run outside the
steady-condition benchmark's timer and inside the scenario timer. Stage fusion
was measured as part of the complete implementation, not in an isolated ablation;
no separate speedup is claimed for it. The selected setting's full 10,000-cell,
18- and one-thread measurements and scenario scores are in [RESULTS.md](RESULTS.md).

The alternatives remain reproducible through `KERNEL_LU`, `KERNEL_BLOCK`,
`KERNEL_FUSION`, `KERNEL_FIXED`, and `KERNEL_WIDTH`; their valid ranges and defaults
are in [DESIGN.md](DESIGN.md). Different cell counts, group sizes, and compilation
settings retain the same chemistry and numerical method.
