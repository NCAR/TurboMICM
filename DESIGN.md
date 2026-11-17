# Design (in progress)

For mechanisms constructed in C, `ChemProcess` owns copied reactants,
products, phase metadata, and a rate descriptor. `chem_build_mechanism`
resolves names against a gas system, constructs indexed reaction tables,
carries species tolerances into state metadata, and turns parameterized
reactants into condition-dependent multipliers.

`tools/config.py` reads CAMP configuration manifests and their JSON or YAML data.
`tools/generate.py` emits static C tables for tutorial reaction types,
including rate parameters and labels for caller-supplied values. `src/rates.c`
calculates rate constants from conditions and custom parameters. `src/chem.c` evaluates mass-action
forcing and its dense Jacobian. The derivative is formed by products excluding
one factor, so zero concentrations are handled without division. The interface
separates the mechanism from the kernel and leaves the choice of solver open.

`src/rosenbrock.c` provides a single-cell autonomous ODE/DAE solver using callback
functions for forcing and Jacobian, five reference Rosenbrock parameter sets,
and reference step-size defaults. To add a parameter set, construct a
`ChemRosenbrockMethod` with stage coefficients, evaluation flags, gamma, and
error order; the existing integrator consumes it. The parameter sets can use
algebraic constraints. It supports dense or sparse no-pivot Mozart and Doolittle
elimination. The solver clamps negative ODE values after each public solve, as
the reference CPU solver does. It allocates scratch space per solve; multi-cell scheduling,
compact sparse matrix storage and additional solver features are not yet implemented.

This is a starting point, not the final architecture.

`tools/import_ts1.py` translates the reference generated C++ into the checked-in
CAMP-style JSON files. The generator excludes parameterized species from the state and emits their
condition-dependent functions; `chem_rates_parameterized` applies their product
to reaction constants. A third body is the special case equal to air density.
Species and custom defaults are also emitted as C arrays. See DECISIONS.md for
the configuration extensions.

The preliminary benchmark schedules individual cells with OpenMP static
scheduling; `--threads` sets the thread count (default: OpenMP maximum). Each
cell currently owns independent solver state and error control. Timings exclude
input, rate setup, warm-up, and output. A mechanism-generic symbolic plan
calculates fill once, then skips structural
zeros during factorization and triangular solves. Matrix values still occupy
dense row-major storage. Sparse-plan clearing and sign changes visit only
structural and fill positions; the benchmark accumulates its Jacobian into
those cleared entries. The plan must cover every possible derivative. The
`--reorder` generator option applies a diagonal
Markowitz permutation before emitting the mechanism; the tests exercise that
path with a separate small mechanism.

`src/backward_euler.c` implements the implicit Newton iteration and reference
step-reduction sequence. Its grouped interface shares convergence and step
acceptance across independent cell blocks. Both solver families
use the same dense or sparse factorization functions; set `lu_algorithm` to
select Mozart or Doolittle. The default is Mozart.

Rate evaluation dispatches through `ChemRateType`; `p` holds parameters in the
order used by `src/rates.c`, and the Taylor type uses `custom_index` for its
coefficient count. To add a rate family, add its enum and evaluation formula,
then add parsing and unit conversion in `tools/generate.py` and exercise both
the direct and generated paths. Callback rates use a caller-owned binding table
indexed by generated callback labels; rebinding does not require regeneration.

The sparse plan can be built from a mechanism or an arbitrary structural
pattern. The plan holds fill and row/column index lists; it is immutable and
may be shared by threads. To replace the linear solver, implement the factor
and solve operations for that plan, add an algorithm choice to `ChemLuAlgorithm`,
and connect it at the factor/solve calls in both integrators.

Constraints are borrowed indexed descriptions validated by `chem_constraint_validate`.
`chem_constraint_residuals` replaces algebraic forcing rows, while
`chem_constraint_subtract_jacobian` adds the negative constraint derivatives
to the supplied matrix. Equilibrium constants can be supplied per cell or
evaluated using the Van’t Hoff parameters. The Rosenbrock solver projects
inconsistent algebraic state with a damped
Newton iteration, then replaces algebraic forcing and Jacobian rows, suppresses
those rows in the mass matrix, and uses algebraic step changes in error control.
The projection rolls back on failure. A caller-supplied sparse plan must include
the union of kinetic and constraint structure; `chem_constraint_sparsity` can
add constraint positions to an existing pattern.

The additive process APIs preserve existing forcing and negative Jacobian
entries and can suppress algebraic output rows without suppressing their
reactant columns. `include/jacobian_check.h` provides finite-difference
Jacobian generation and comparison against a declared sparsity pattern.

The optional `ChemState` owns cell-major concentration, rate, and custom
parameter arrays, conditions, names, and tolerances. Deep clone and explicit
move functions make ownership visible in C. Bulk named updates validate every
column before changing the state; the solver interfaces can also use the arrays
directly without allocating a `ChemState`.

`chem_jacobian_pattern` and `chem_used_species` OR structural information for
state species into caller-owned masks. Sparse plan construction shares the
same Jacobian pattern operation. Parameterized species are represented through
rate multipliers and are not included in the state-species mask.

The generated `ChemSpeciesCatalog` retains configuration species names and
marks those occurring in reactions, including parameterized reactants and
products that are not stored as state variables. It is independent of the
solver permutation. Callers constructing mechanisms in C can provide the same
borrowed catalog representation.

Custom constraints borrow a callback context and dependency list through
`ChemCustomConstraint`. Their residual and negative Jacobian callbacks receive
the state, temperature, and optional per-cell parameter override. They can be
combined with built-in constraints in one set; declared dependencies are added
to the sparse pattern before constructing a plan. Context and dependency
storage must outlive the solve.

`ChemMatrix` provides optional owned row-major and interleaved group layouts.
A group width of zero selects row-major storage; a positive width pads the
final group and interleaves cells within columns. Row and column operations
validate dimensions, while the low-level index accessor supports direct use
of the layout. This is a separate utility; the current solver keeps its
existing per-cell layout.

Borrowed column views can refer to a matrix column or a caller-owned vector.
`chem_columns_apply` checks row counts before invoking a reusable row callback.
`chem_column_reduce` supports either logical rows or padded group storage,
letting callers choose the identity values used in padded reductions.

`ChemMatrixFunction` records the column counts and number of matrix and
vector arguments at creation. Calls may use a different row count, but all
arguments must agree on it; validation completes before the callback runs.
The caller owns the function and destroys it when it is no longer needed.

`ChemSparseMatrix` owns one shared sparsity pattern and a `ChemMatrix` of
block values. Its CSR or CSC order controls entry numbering, and the matrix
group width controls how block values are interleaved. Duplicate coordinates
are coalesced; absent and out-of-range entries return an invalid index. This
utility is separate from the current LU plan, which retains dense value buffers.

Sparse entries expose borrowed block views that share the dense column-view
interface. Views can be selected by row and column or by the physical index
of an entry in block zero. Column fill, copy, and row callbacks operate on
logical blocks, including across sparse, dense, and vector storage.

Sparse block functions use the shared matrix-function signature over their
value matrices. Signatures validate compatible group widths at creation and
again at invocation. An empty invocation returns without calling the callback.

Group spans expose the contiguous lanes of a dense column or sparse block
entry, with separate logical length and padded capacity. Const spans allow
read-only access; group fill and copy cover the allocated lanes, including
tail padding, while callers can restrict calculations to the logical length.
