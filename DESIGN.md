# Design (in progress)

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
dense row-major storage. The `--reorder` generator option applies a diagonal
Markowitz permutation before emitting the mechanism; the tests exercise that
path with a separate small mechanism.

`src/backward_euler.c` implements the implicit Newton iteration and reference
step-reduction sequence. Its convergence function supports multiple cell rows,
while the integration function currently advances one cell. Both solver families
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
