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

- Matrix template policies become runtime C element types and group lengths.
  CSR/CSC sparse matrices expose a blocks-by-nonzeros array, sharing callback and
  view primitives with dense matrices and vectors. This preserves shape errors,
  aliasing, padding, and strict/non-strict reductions without C++ type machinery.
- Matrix function creation captures columns and types but does not require equal
  row counts. Invocation performs row-count validation. This distinction is
  explicitly tested in the reference and is retained in C.
- Five cases are N/A because they contain exclusively C++ template/concept tag
  assertions. Their individual IDs and reasons are in tests/not_applicable.json.
  Runtime vector, dense, sparse, and grouped-view cases are ported and executed.
- Sparse/array invalid dimensions and structural-zero accesses return CHEM_INVALID
  (or SIZE_MAX from index lookup) instead of an exception or a debug-only death.
  Checks remain enabled in Release. Read-only metadata and explicit callback
  write masks replace C++ const-reference overload selection.
- Seeded random test inputs reproduce the supplied host build's minstd_rand,
  MT19937, canonical uniform, and cached polar-normal sequences in local Python
  tooling. This fixes the source tests' implementation-defined random choices
  across C compilers and operating systems while retaining their matrix sizes,
  distributions, block counts, and numerical assertion bounds.
- Separate LU factors share the in-place algorithm's operation plan. They have
  an explicit unit lower diagonal and zero unused triangles, preserve the input,
  and overwrite every fill value. The temporary factor buffer belongs to the
  invocation, so independent calls can run concurrently.
- Species and phase metadata use owned C objects and typed properties. Missing
  properties and unsupported property types have distinct error codes. An empty
  state is represented by creation with a NULL mechanism and zero cells; explicit
  clone, assignment, and move operations replace the corresponding C++ mechanics.
- The `parameterization` species extension contains `constant`, `temperature`,
  `pressure`, and `air density` coefficients. Such species do not occupy evolving
  state columns; their reactant contributions become rate multipliers. The
  existing `third body` extension keeps its direct air-density multiplication.
  `configs/parameterized` exercises the full generated-versus-code path.
- Version getters describe the supplied reference API and use the values in
  its actual include/micm/version.hpp (3.14.0), rather than inferring generated
  version metadata from the experiment description.
