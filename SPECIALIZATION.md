# Mechanism specialization

`tools/specialize.py` reads the Tier 1 YAML/JSON format through `tools/configure.py`. The generator fixes species/parameter indices, the Markowitz permutation, ordered stoichiometry, negative-Jacobian derivatives, symbolic fill, exact elimination/substitution operations, built-in kinetic parameters, and the three-stage execution shape. It emits ordinary C11 and CPU OpenMP loops. Generated C is an ignored build product; Make tracks configuration, included configuration files, both generators, templates, and generation options.

Runtime data are conditions (including independently supplied air density), concentrations, all custom rate parameters, step length, tolerances, cell/group counts, worker count, and scheduling. The six unusual TS1 photolysis defaults are overridden only by input, exactly as in Tier 1. Nothing in the generator identifies or approximates TS1 chemistry. Built-in rate coefficients change through configuration and regeneration; user parameters remain writable at every call. Derived branched/surface constants use the generic loader's once-computed C doubles to retain libm roundoff. There is no runtime configuration parser or Python dependency.

The default linear backend emits the same Mozart reciprocals, scalings, ordered updates, and triangular substitutions as the generic library. It fixes every sparse index and partitions the generated code into helpers. An optional left-looking backend accumulates each entry with its identical ascending-pivot subtractions. No backend drops fill or uses an approximate Jacobian, factorization reuse, or reduced precision. The three-stage coefficients and controller are shared with the generic library.

The useful layout change is a four-cell execution tile inside each 32-cell adaptive group. A tile packs the already computed rates, old concentrations, and retained initial forcing into a small workspace, evaluates its Jacobian, factors, completes the three stages, and copies its new state/error terms back. Independent cells permit this stage interleaving. The full group's species-major serial sum still controls acceptance and the next step. A rejected attempt recomputes the Jacobian/factorization and retains initial forcing, exactly as before. Scratch fits the three unused stage vectors for a full TS1 group; smaller groups or larger tiles use a checked allocation. Group boundaries, padding exclusion, warm-up behavior, and logical statistics remain unchanged.

Chemistry kernels can fuse a chunk of ordered reaction/derivative blocks inside one cell loop. First-write initialization replaces the first accumulation with `0.0 +/- term`, preserving arithmetic order, and explicitly zeros untouched species and matrix fill. The stage pipeline fuses vector passes, computes squared scaled error terms with SIMD, and retains the original serial summation order. These options, rate specialization, linear layout, helper sizes, literal cell-stride backend, and tiling are build parameters in `DESIGN.md`.

## Experiments

All local timings below use clang with `-O3 -DNDEBUG`, without extra optimization/target/floating-point flags. The small benchmark is 100 cells from the supplied seed-1 CSV, one worker, one warm-up and 100 timed steps of 30 seconds. The multicore smoke test is seed 2026, 10,000 cells, 18 workers, one warm-up and 20 timed steps of 30 seconds. These are individual exploratory runs; a few percent can be run variation. Every listed completed variant passed the supplied 100-cell MICM gate, with maximum relative difference `3.48e-14`.

| Variant | Small elapsed ms | Multicore smoke ms | Decision |
| --- | ---: | ---: | --- |
| Tier 1 baseline | 377.591 | 590.111 | Retained independently |
| Unrolled LU, fused 8-pivot/row loops, variable stride | 886.012 | — | Excessive generated working set; dropped as default |
| Separate operation loops, chunk 8, variable stride | — | — | Clang compilation stopped after more than 260 seconds; no runtime claim |
| Separate loops with distinct restricted element pointers, chunk 1 | 388.182 | — | Correct, little benefit; retained as optional backend |
| Left-looking register accumulation, chunk 8, variable stride | 396.439 | 669.730 | No improvement; optional backend |
| Fused 8-pivot/row loops with a literal 32-cell backend | 502.501 | 768.356 | Smaller addressing cost, still slower than baseline |
| Literal 32-cell left-looking backend, chemistry fusion/first writes/vector error | 414.442 | 630.427 | Combined transformation, not an isolated attribution |
| Four-cell tiles, row loops, chunk 8, fused initialized chemistry | 317.886 | 484.252 | Main improvement; retained |
| Two-cell tiles, otherwise same | 457.955 | 723.620 | Too small on this compiler/machine; dropped as default |
| Four-cell tiles with left-looking LU | 329.369 | 500.600 | Mozart row loops faster; dropped as default |
| Four-cell tiles, one pivot/row per helper | 309.089 | 466.880 | Selected default |
| Same, chemistry helper chunk reduced from 32 to 8 | 306.338 | 477.512 | No consistent benefit; retain chunk 32 |

Diagnostic kernel timings used 5,000 repetitions on one 32-cell state. The initial fused 8-pivot backend took 991.942 ms for factorization versus 222.192 ms for Tier 1, and 118.011 versus 63.367 ms for a triangular solve. Separate restricted-pointer loops reduced those to 215.624 and 62.821 ms (matching Tier 1 rather than materially improving it). Left-looking accumulation took 260.501 and 69.775 ms. These measurements motivated blocking the complete stage pipeline rather than only changing sparse operation syntax.

Rate specialization reduced 5,000 complete 32-cell rate updates from 593.684 to 128.645 ms in the restricted-pointer experiment, about 4.6 times faster. These updates occur outside the primary benchmark's timer and inside the scenario timer. The speedup therefore helps day-night runs and setup rather than explaining primary benchmark improvements. Native `tests/kernel_bench.c`, built with `make kernel-bench`, can repeat this diagnostic; no Python executes in its timed kernels.

The unblocked, separate-loop, literal-stride, left-looking, and smaller-tile experiments remain selectable to support other mechanisms/architectures and reproduce decisions. No ordering search, approximate chemistry, coefficient fitting, compiler-flag tuning, or architecture-specific implementation was used. Final full-size measurements, build costs, both-compiler gates, and scenario scores are in `RESULTS.md`.
