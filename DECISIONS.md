# Decisions

- C11 is the runtime language. Python is confined to generation and testing in
  the local prefix environment. The C library owns copies of mechanism metadata;
  states borrow an immutable mechanism and own their numerical storage.
- Stable rate-type sorting follows MICM's type order. Species initially retain
  configuration order, excluding parameterized species. A parameterized reactant
  multiplies its rate constant once per occurrence.
- States store each variable contiguously across cells. This supports portable
  OpenMP SIMD and will allow independent groups with separate adaptive steps.
- The configuration reader accepts the tutorial CAMP YAML and JSON structure.
  Default rate-coefficient units are molecule cm-3, converted to mol m-3 using
  1e-6 times Avogadro's constant. Custom rates retain caller units, matching the
  by-hand tutorials. Surface diffusion and molecular weight use their stated SI
  units. Branched reactions expand into alkoxy and nitrate reactions.
- Format extensions: a mechanism or reaction may declare `rate units: mol m-3`;
  species may declare `parameterization: air density`; reactants may be an
  ordered list of names and products an ordered list of `{species, yield}`
  entries, retaining repeated terms and arithmetic order; a mechanism may carry
  `initial state`, a mapping of temperature, pressure, concentrations and custom
  rate labels to values. These extensions represent TS1's existing data without
  converting coefficients twice or losing order. Additional analytic rate types
  use the parameter names of the corresponding C structures.
- `tools/translate_ts1.py` extracts data from the supplied C++ declaration forms,
  asserts the full species/reaction counts, and writes the retained JSON input.
  It neither compiles nor executes source code. The generic configuration reader
  then generates C builders through the same public API used by hand-built tests.
- A missing test implementation is a FAIL, never an N/A. N/A will be used only
  for individually reviewed cases with no meaningful C behavior.
- The supplied task prompt is an input, excluded from the submitted source tree.
- C++ standard random engines and distributions are implementation-defined.
  The C LU/linear test ports use a fixed Park-Miller uniform generator and
  Box-Muller lognormal draws with the original distribution parameters, matrix
  dimensions, sparsity probability and bounds. They test reconstruction of
  every matrix element or right-hand-side element, including structural zeros.
  Vector-policy cases test groups of 1, 2, 3 and 4 cells with partial tails.
- Sparse and dense layouts are selected by the C API instead of C++ template
  policies. State group size replaces compile-time vector length and is a
  runtime parameter. No floating-point compiler options are added.

## Constraint representation

The C constraint API uses tagged linear/equilibrium definitions with named
terms. A mechanism owns deep copies and resolved variable indices. Installing
constraints is transactional and must precede state and solver creation, just
as mechanism reordering does. The mass matrix is represented by an algebraic
row mask. No mechanism-specific code is generated for DAE enforcement.

Constraint tests retain all original cells, coefficients, temperatures and
bounds. Vector policy cases use the C contiguous-cell layout; group-size cases
split the original cells into groups of the stated size, including partial
groups. There are no inactive padded lanes to include in reductions. The
initialization failure status string preserves MICM's exact spelling.

## Separate factors and process accumulation

Out-of-place LU keeps input and output buffers separate, while storing L and U
in one output buffer with an implicit unit lower diagonal. This retains MICM's
Mozart and Doolittle arithmetic and avoids duplicating numerical kernels solely
for storage layout. Triangular accessors support reconstruction tests. The
rejected-step regression installs separate Doolittle factors using the public
linear-solver callbacks and compares its statistics with in-place Mozart.

The public kinetic accumulation functions leave the destination's initial
values intact and skip algebraic rows. Full forcing/Jacobian functions clear
outputs and then apply kinetics and constraints. This preserves ProcessSet's
additive contract without forcing callers to allocate a temporary array.

## State ownership and dimensions

C pointer-to-pointer assignment/move helpers replace C++ ownership operators;
column accessors replace square-bracket proxies. States borrow immutable
mechanism metadata, which must outlive them, and own all numerical buffers.
The default-constructor tests use a zero-variable mechanism and zero-cell
state. Other state tests create ten reactions with the original variable and
custom-label dimensions.

MICM's named setters reject a scalar for a multi-cell state. The initial C API
broadcast scalars implicitly; the state-case port corrected that mismatch.
Explicit fill APIs now provide broadcasting for configuration defaults and
hand-built configuration comparison fixtures. The original error checks and
all state tests pass without relaxed bounds.

## Matrix utility API

`ChemMatrix` selects real or integer storage and row-major (`group=0`) or
grouped column-vector (`group>0`) layout at runtime. Grouped allocations retain
padding so the original storage offsets and non-strict reductions are observable.
The solver's existing contiguous-cell buffers are unaffected by this utility API.
Owning outputs start zero-initialized and must be freed before reuse; borrowed
arrays/views never free data. Typed element pointers provide direct native
integer/real access. Numeric callbacks use double values, with checked integer
conversion on write.

C callbacks plus borrowed context replace lambdas and their captures. Function
creation records argument kinds, scalar types, column counts and mutability;
invocation checks equal row counts while allowing them to differ from creation.
Row callbacks receive a snapshot of arguments and write back mutable arguments.
Their default iteration excludes padding so caller-owned vectors need no extra
allocation. Reductions explicitly select logical or padded rows. Const view
writes and out-of-range columns return errors, including release builds where
the corresponding C++ assertion would be disabled. Tests retain all original
checks, add those release-mode errors, and use groups 1, 2, 3 and 4 as specified.

## Species, phase and system metadata

Owning `ChemSpeciesInfo` property maps retain separate real/string/integer/bool
entries, including same-name entries of different types. An unsupported enum
tag tests MICM's invalid-property-type error without using extended floating
point, which the portability rules forbid. Affine species parameterizations
retain the original constant, temperature, pressure and density coefficients.
Phase entries own species copies and optional diffusion/density values. Named
reaction definitions own their phase, indices and rate label. Their copy and
assignment helpers replace C++ value semantics.

`chem_system_mechanism` resolves reaction names against the system gas phase
and calls the existing generic mechanism constructor. The resulting mechanism
owns a separate system snapshot; parameter callbacks point into that snapshot,
so no builder metadata must outlive it. Both construction paths therefore use
the same runtime rate and solver code. The OpenMP test shares this immutable
mechanism and gives each thread its own state and solver scratch workspace,
matching the original state-isolation behavior with the C ownership design.

The supplied `micm/include/micm/version.hpp` returns version `3.14.0` and numeric
components 3, 14, 0, 0, despite the older version description in the task prompt.
The compatibility version API follows those actual supplied source values.

## Sparse utility layouts and functions

`ChemSparseMatrix` preserves the declared structure, with no forced diagonal or
factorization fill. Its builder deduplicates entries and sorts them by compressed
row or column order. Values use the existing dense matrix buffer with one row per
block and one column per stored element; group zero is standard block-major and
positive groups retain the original padded offsets. Sparse lookup uses a binary
search over the sorted structure. Print output uses logical row/column order.

Sparse functions use the same callback/view engine as dense functions. They
require matching matrix group widths while allowing vectors without padding,
matching MICM's ordering validation. They accept differing sparse structures and
validate block counts at invocation. A raw block-zero vector index identifies a
block view; missing entries and missing block indices return distinct errors.
All debug-only source bounds assertions are also checked in the C release tests.
Empty function invocations validate their shapes and skip callback execution.

## Grouped views and category tests

`chem_view_group` selects one physical group and precomputes its element base
pointer. A category enum distinguishes vector, dense, sparse, grouped dense and
grouped sparse views. It replaces the source's template category traits with
observable C constructor/classification behavior, including static assertions
for distinct enum values. No category case is declared N/A.

Raw-view operations still visit only logical rows. Explicit grouped-view fills,
copies and non-strict row loops visit the allocated padded lanes, matching the
source grouped primitives. Strict row iteration and strict reductions exclude
those lanes. Tests verify group boundaries, pointer offsets, untouched neighbor
entries and partial tails, with additional checks of padded versus strict visits.


## External processes and constraints

C model descriptors replace the source's model template packs. Each build owns
resolved index arrays and a separate context copy. This preserves repeatable
builds without retaining descriptor strings or moving data out of the caller.
A callback gets the current binding instead of caching a particular sparse
layout. Positive derivative accumulation uses a destination view with a sign,
so a single callback serves dense numerical checks and the negative sparse
solver Jacobian. Generic callers can supply arbitrary nonlinear models.

External process entries in algebraic rows are retained structurally, then
cleared numerically before applying constraints. This covers both externally
controlled rows and rows controlled by built-in constraints. Built-in constraint
replacement retains model parameters and algebraic rows and rejects duplicate
controlled variables. Constraint count/row enumeration includes both sources;
only built-in rows have a ChemConstraint definition.

The 18 aerosol cases retain the source initial values, quantities and analytical
tolerances. The C tests exercise both LU algorithms with separate/in-place
factors, all five Rosenbrock methods, backward Euler, and batch widths 1–4.
As in the source, the vector Rosenbrock policies use ROS3; other Rosenbrock
methods use the full original one- or three-cell state. An exploratory ROS6
split into independent one/two-cell states changed its adaptive error norm and
missed the original three-cell analytical tolerance. That combination is not
one of the source cases and is not used to represent it. No method coefficient,
solver tolerance or original acceptance bound was changed. Solver matrices have
a canonical layout; CSR/CSC data-view behavior is covered by the separate sparse
utility ports rather than duplicate integrator implementations.


The Terminator integration port uses the same portable Park–Miller/Box–Muller
random convention as the linear-algebra ports, with the original lognormal mean
-2 and standard deviation 2 and concentration scaling. It retains all cell
counts, twenty longitude samples, latitude-dependent photolysis, 30-second
intervals, 1e-8 solver relative tolerance, 1e-20 absolute tolerance and the
original closed-form concentration checks (`1e-14 + 1e-6*abs(reference)`). The
scalar and vector cases both exercise the column-major cell loop; the scalar
case retains the source's repeated fixture runs.
