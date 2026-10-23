# TS1 specialization

The generator reads the same YAML/JSON mechanism format used by the generic
library and emits both the mechanism descriptor and a solver entry point.
`make bench` regenerates the source from `config/ts1.json`; generated C is a
build product and is not edited by hand. The two small configuration
mechanisms in `make test-specialized` use the same generator and solver path.

## What becomes code

- Reaction topology, reactant/product indices, stoichiometry, yields, and
  integrated species mapping are fixed in the generated source.
- Each supported rate kind becomes its own formula. Fixed Arrhenius, falloff,
  branching, tunneling, Taylor-series, reversible, surface, and user-rate
  constants are embedded in expressions. Cell temperature, pressure, air
  density, and user-defined parameters remain inputs.
- The Jacobian sparsity pattern and a deterministic low-fill species order are
  calculated at generation time. The generator simulates symbolic elimination,
  then emits the exact Mozart factor updates, row structure, and triangular
  solves. No numerical Jacobian or LU factorization is approximated.
- Reaction forcing and Jacobian contributions are emitted per reaction rather
  than dispatched through runtime reaction descriptors. Species state and
  solver work are species-major across cells in each group; the benchmark's
  cell-major input and output are converted at the group boundary.
- The three-stage Rosenbrock coefficients, error weights, and step controller
  match the generic implementation. State values, parameters, tolerances,
  interval length, and scenario photolysis values remain runtime data.

## Choices and measured effects

The final measured results are recorded in [RESULTS.md](RESULTS.md). The
baseline is the Tier 1 generic solver on the same 10,000-cell input, compiler,
thread count, and optimization flags. The table includes specialized versus
generic and MICM timings at 18 and one thread, plus generated source size and
clean build duration.

The main experiments were:

| Specialization | Decision and observed effect |
| --- | --- |
| Generated rate formulas, reaction forcing, and Jacobian | Kept. This removes runtime rate-kind and reaction descriptor dispatch while preserving per-cell conditions and parameters. The combined generated solver's measured gain over Tier 1 is 7.0% at 18 threads and 18.1% at one thread. The timing does not isolate this part from the other specializations. |
| Fixed sparse pattern, symbolic ordering, and triangular solve rows | Kept. The minimum-actual-fill ordering is deterministic and computed once by the generator. Its measured timing contribution is included in the combined solver result; it was not benchmarked as a separate ablation. |
| Direct generated pivot updates | Kept. Each pivot has a straight-line exact update sequence. In an earlier single-run comparison before the final ordering change, table-driven sparse pivots took 10,174.840 ms and per-pivot generated updates took 10,245.466 ms at 18 threads. The 0.7% difference was too small to establish a reliable speedup; the fixed operations were retained as part of specializing the exact solve. The final combined result is 9,975.272 ms. |
| Species-major state inside a cell group | Kept. It makes each cell-vector loop contiguous while retaining the public cell-major CSV interface. Its measured timing contribution is included in the combined solver result; it was not benchmarked as a separate ablation. |
| Runtime tolerances and group size | Kept as run settings. Scenario groups default to the gate's 32 cells; changing tolerances or grouping does not require regeneration. |
| One fully unrolled monolithic factor routine | Dropped. Its generated translation unit exceeded 2.8 MB and clang had not finished optimizing it after about 70 seconds. There is no runtime measurement for that version. Splitting the fixed pivot operations into per-pivot functions produced the final 2.4 MB source and a 62.03-second clean clang build. |

The generated TS1 source is about 2.4 MB. It is substantially larger than the
generic C library, and its clean build takes longer because the compiler sees
the fixed factorization and reaction operations. This is an accepted tradeoff
for the frozen TS1 production mechanism; other mechanisms can still use the
generic library directly.
