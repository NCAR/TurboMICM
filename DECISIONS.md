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
- Missing case ports are reported as FAIL. All applicable cases are now ported;
  N/A is reserved for five documented cases with no meaningful C counterpart.

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
- Analytical test fixtures retain the original source invocations and bounds.
  C++ matrix-policy variants that have the same C representation share an
  execution of the corresponding method/LU/group setting at the tightest bound.
  One-cell helpers share a runtime representation for every group capacity.
  CSR/CSC array storage and separate/combined LU behavior are also independently
  exercised by their full matrix and linear-algebra case ports.
- Analytical expected concentration tables are retained as numeric test data.
  The tests run all original output times, including 999 one-second chain
  solves, the long Robertson/E5 intervals, and 18,000 backward Euler substeps per
  Oregonator output. Test execution never loads the reference implementation.

- External model configuration is a two-phase C interface: declare state names,
  then resolve indices against the merged descriptor. A finalization callback
  runs after sparse construction. Builder-owned merged tables make repeated
  builds independent, while borrowed model contexts and equations retain explicit
  caller-managed lifetimes. The C descriptor preserves the supplied state order;
  generic Markowitz ordering is separately available at build time and at runtime.
- Algebraic filtering applies after external processes and before constraints.
  External Jacobian coordinates remain in the structure even when their numeric
  process contributions are cleared. Finite-difference, composition, and ordering
  cases exercise this distinction.
- The aerosol invocation fixture records every original method, LU, group, and
  bound. Initial vectors and temperature-dependent parameters are copied from
  the source helpers. The flow-tube CSV fixtures come from the source regression
  expected_results directory, including all 3,601 rows.
- The legacy Chapman regression oracle is ported independently in test-only
  Python, with closed-form derivatives, dense Gaussian elimination, and the
  original three-stage controller, start step, and tolerances. Production C has
  no Chapman-specific numerical kernels. Random inputs retain the source's
  extra concentration fill and padded vector storage consumption.
- The six remaining tutorial programs are C programs using the public library.
  Their output is compared with retained output of the supplied executables.
  Only measured wall-clock durations and unsynchronized thread startup messages
  are excluded from exact text comparison; concentrations and solver statistics
  match the printed reference results. The shared-solver case separately requires
  bit-identical results from eight concurrent states.
- The vectorized tutorial's source Solve loop initially tests a false condition.
  Its C port preserves the resulting parameter updates and unchanged state,
  while substantive vectorized integrations are checked in the other case ports.
  Integer labels encode its sparse string values so both storage orderings use
  the existing typed C sparse API and print the same labels.

- Tier 2 reuses the configuration normalizer and the Rosenbrock algorithm. The
  generator emits numerical kernels rather than maintaining a second adaptive
  controller. Generation-time symbolic fill is independently checked against
  every coordinate of the generic library's runtime plan.
- Fixed rate descriptors are immutable in a specialized solver. Runtime custom
  parameters, including surface, emission, loss, heterogeneous and photolysis
  inputs, remain mutable arrays. Changing a fixed coefficient or adding a
  reaction means editing the configuration and rebuilding.
- Full groups can have constant strides in generated kernels. Other cell counts
  use the generic kernels on the same descriptor and CSR representation. This
  preserves arbitrary group sizes and partial groups without padding the error
  norm or adding artificial cells.
- The scenario file format is interpreted as documented: paths are relative to
  the scenario, interpolation is linear with constant endpoints, comments start
  at `#`, and absent snapshot intervals request final-state output only. An
  explicit interval always includes the last step, even when it is not divisible
  by that interval. Scenario defaults retain 32-cell groups even on small inputs.
- Snapshot copies and output are outside timed segments. Each segment runs its
  steps within a static parallel group loop, preserving independent per-group
  trajectories while timing actual rate updates and solves.
- Empty reaction sets get a C11 dummy descriptor element with a reaction count of
  zero. Tests include inert species, high-order repeated reactants, self-products,
  source/loss rates, condition multipliers, and every supported rate family.
