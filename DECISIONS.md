# Decisions

- The implementation uses double precision and C11. The numerical formulas are
  ported from the specified reference revision; no reference code is linked into
  the C library or benchmark.
- Configuration reading and code generation run only at build time, in the local
  prefix environment. JSON and YAML share one normalized representation.
- The TS1 translation preserves reaction order and floating-point constants.
  Extensions to CAMP configuration: species `third body`, `initial concentration`,
  mechanism `initial conditions`, and `initial rate parameters`. These represent
  the reference's parameterized air species and InitState defaults explicitly.
- The user corrected the delivery instruction: commit locally and do not push.
  The supplied repository has no remote.
- PROMPT.md is an experiment input and is ignored rather than committed because
  it contains identifying names forbidden in deliverable files.
- Species reordering is a build-time generic Markowitz pass. Stable original-index
  ties match the reference default; the C library itself also accepts unreordered
  descriptors and computes their sparse fill.
- CAMP built-in numeric coefficients default to molecule cm-3 units, inferred
  from the tutorial's explicit conversions. `rate units: mol m-3` is an extension
  on a mechanism or reaction for already converted inputs, including translated
  TS1. Troe low/high coefficients scale with reaction order/order minus one;
  ternary activation scales with order minus one/order minus two. User-defined
  parameters are supplied directly in library units, as in the Robertson tutorial.
- Reactant and product lists are an extension permitting repeated reactants and
  preservation of source arithmetic order. Original mapping forms remain valid.
- The reader preserves source reaction order. Branch expansion occurs at the
  source position. Independent emissions and losses can have different ordering
  in examples and hand-built mechanisms; equivalence is tested numerically.
- Air density is explicit state data, as in the reference. Call
  `chem_calculate_air_density` when ideal-gas density is desired. The benchmark
  does this after applying each CSV row.
- Post-solve nonnegativity is part of the public solver behavior. A retained TS1
  regression cell detects omission of this wrapper operation; checking convergence
  status and the usual 100-cell fixture alone did not detect that omission.
- Pending case ports are printed as FAIL. N/A will be used only for a documented
  case with no meaningful C counterpart. Current results are explicitly partial.
