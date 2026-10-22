# Mechanism specialization

The generator specializes the chemistry and exact sparse algebra for the input
configuration. TS1 is an input, not a branch in the generator: it has 209 solution
species and 547 reactions, and its ordered LU pattern contains 2,381 entries.
The same generator is tested on the three Tier 1 tutorial configurations, the
additional-rate configuration, and an automatically modified configuration.

## What becomes code

The normalized reaction list fixes species indices, reactant occurrences,
parameterized-reactant multipliers, yields, rate types and fixed coefficients.
Python computes the generic library's Markowitz permutation and all LU fill.
Generated forcing and Jacobian expressions contain literal offsets; LU and both
triangular sweeps contain the exact elimination operations instead of traversing
sparse plans. No dense inverse, approximate Jacobian or reduced precision is used.

The main group width is also a compile-time stride. A group with a different
width uses the generic kernels; the adaptive controller still sees exactly its
real cell count. `--width 0` emits variable-stride specialized kernels instead.
The normal driver default and the scenario gate are both 32 cells per group.

Rate formulas have literal coefficients, enabling ordinary constant folding of
zero exponents and unit factors. Two prepared branched constants and the surface
speed constant are computed once in C and read from immutable metadata. This
keeps the platform's libm evaluations consistent with the generic library.

The three-stage integration loop and its final-state/error loops are specialized
through the same source template as the generic solver. The final-state and
error additions are fused across traversal while retaining their original order
for each element. Error normalization retains the generic, sequential RMS sum.
The numerical coefficients, step controller, rejection behavior and interval-end
clamp remain MICM's.

## What stays data

Temperature, pressure, air density, concentrations, every custom rate parameter,
interval length and tolerances remain runtime data. This includes all photolysis,
heterogeneous and surface parameters, even when the benchmark leaves them fixed.
The exact-input rate cache snapshots conditions and the complete custom matrix.
It skips rate evaluation only for byte-identical inputs, and can be disabled with
`--no-cache`. It never caches stage concentration products, Jacobians or LU
factors. The generic updater invalidates the cache and every state owns its cache.

Metadata construction, name-based initialization, allocation and CSV handling
reuse the Tier 1 library. They are outside the benchmark timer. Solver metadata
is immutable after attachment to the generated functions.

## Measured experiments

These are development measurements on the supplied 18-core machine with Apple
Clang 21, `OPTFLAGS='-O3 -DNDEBUG'`, 10,000 cells from seed 2026, one warm-up and
100 timed 30-second intervals. Each row is one run with 18 threads. No compilation
or other validation process ran concurrently with these measurements. They are
not idle-machine medians, and small differences should not be treated as robust.
The final independently rerun benchmark results are in RESULTS.md.

| Variant | elapsed_ms | Decision |
| --- | ---: | --- |
| Generic Tier 1 library, groups of 32 | 3,927.40 | Baseline |
| Generated rates, forcing, Jacobian, LU and solves; variable stride | 4,131.42 | Unrolling alone was slower |
| Variable stride, fuse eight operations per cell loop | 5,068.31 | Drop as default |
| Fixed stride 32, one reaction/derivative/solve row per loop | 3,327.11 | Retain |
| Fixed stride 32, fuse eight operations | 3,616.27 | Drop as default |
| Fixed stride 32, four-cell helper tiles | 5,385.95 | Remove; also much slower to compile |
| Fixed stride 32, separate LU-update loops, exact rate cache, dynamic stages | 5,617.32 | Keep fused-pivot LU |
| Fixed stride 32, fused-pivot LU, rate cache and fused three-stage estimate | 3,045.69 | Retain as default |
| Same, groups/stride 16 | 3,029.91 | Near tie; keep 32 to match the scenario gate |
| Default plus combined initial forcing/Jacobian pass | 3,122.68 | Optional, off by default |

The initial generation experiment specialized all four arithmetic components
together; it did not separately attribute timing to rate folding, fixed chemistry
indices, unrolled LU and unrolled substitution. Likewise the cache and stage
specialization were measured together in the retained variant. The table does
not claim independent effects for those components. Fixed stride improved the
initial generated variant by about 19%; the retained cache/stage bundle improved
the fixed-stride variant by about 8%. The four-cell and split-LU experiments
illustrate why larger generated code is not automatically faster.

The initial variable-stride, fixed-stride, fusion, tiling and split-LU variants
passed the supplied 100-cell final-state check. Cache/stage variants also passed
the generated-versus-generic integration tests. The retained default passes the
full benchmark and scenario gates on both compilers, as recorded in RESULTS.md. Tests also ran variable-stride
kernels with 17 operations per helper and three fused operations, covering
nondefault chunking/fusion and all tested cell counts.

## Reproduction and extension

Generated C is a build product under `build/specialized/`. `make bench` regenerates
it from CAMP files and compiles it with `CC` and `OPTFLAGS`. Transitive include
files and generator settings are make dependencies. See DESIGN.md for every
parameter and README.md for the command lines.

The generator is organized into validation/normalization (shared with Tier 1),
symbolic planning, rate expression emission, chemistry emission, factorization,
triangular solves and integrator instantiation. Adding a reaction or changing a
fixed coefficient only changes the configuration. The mutation test does both
and also exercises a new species, repeated reactants, a density multiplier,
a source reaction, fractional yields and a shared custom label.
