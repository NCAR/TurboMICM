# Design

## Runtime layout

The public C interface is in `include/kinetics.h`. A generated mechanism is
immutable table data: species, flattened reactant and product lists, yields,
rate parameters, and custom parameter defaults. Concentrations and conditions
belong to the caller and use cell-major concentration rows. Rate constants are
computed once when a solver is created because cell temperature, pressure, and
custom rate parameters remain fixed during chemistry integration.

The rate evaluator handles Arrhenius, Troe, ternary chemical activation,
JPL-19 ternary activation, branched, tunneling, Taylor-series, reversible,
user-defined, surface, and C callback rates. A third-body species has no state
column; its reactant count multiplies the rate by the corresponding power of
ideal-gas molar density. A surface rate uses the named species' molecular
weight and diffusion coefficient and two custom parameters for particle radius
and number concentration.

`include/matrix.h` provides strided dense views and compressed CSR/CSC storage
with grouped block values. These routines expose the storage and shape rules to
C callers. The chemistry solver uses a specialized sparse pattern and factor
workspace in its hot path.

## Configuration to C

`generate_c.py` follows a configuration's `camp-files`, reads its species and
mechanism objects, resolves names to integer indexes, expands integer reactant
quantities, and emits static C tables. Product quantities remain explicit
yields. The reader rejects unknown species, unsupported rate names, noninteger
reactant quantities, and surface reactions without required properties.

The TS1 translator reads the generated species list, reaction builders, and
initial-state setters from the read-only TS1 header. It writes the same
configuration schema consumed by the regular generator. No TS1-specific
calculation is emitted into the C runtime.

## Integration and linear algebra

The solver implements MICM's default three-stage Rosenbrock method, embedded
error estimate, step-size limits, acceptance logic, and tolerances. Its public
coefficient structure accepts up to six stages, a replaceable controller, and
replaceable factor and solve callbacks. For each
mechanism, it builds a sparse Jacobian graph from reaction dependencies, then
adds symbolic fill for no-pivot Gaussian elimination. A diagonal Markowitz
ordering is selected automatically when C tables are generated. At run time,
the solver assembles the sparse negative Jacobian, factors
`alpha * I - J` in place, and reuses the factorization for all three stages.

Cells are split into contiguous groups. By default, the group count equals the
OpenMP maximum thread count. `--threads N` controls that count, while
`--group-size G` sets a fixed number of cells per group. A static OpenMP loop
assigns each group to one thread. Every group owns its adaptive step-size
history and statistics; statistics are summed after all groups finish.

## Runtime tuning parameters

| Parameter | Default | Override |
|---|---:|---|
| Compiler optimization flags | `-O3 -DNDEBUG` | `OPTFLAGS=...` |
| OpenMP thread count | runtime maximum | `--threads N` |
| Cell group size | balanced, one group per thread | `--group-size G` |
| OpenMP group scheduling | static | compiler/runtime OpenMP setting |
| Rosenbrock method | MICM three-stage default | C solver parameters |
| Relative tolerance | `1e-6` | C solver parameters |
| Absolute tolerance | `1e-3` per species | C solver parameters |
| Initial internal step | `1e-6 * DT` | C solver parameters |
| Minimum internal step | `1e-15 * DT` | C solver parameters |

These settings affect step selection or scheduling, not reaction definitions.
The gate checks the numerical result independently of compiler optimization
flags or cell grouping.

## Extending the implementation

- **Add a rate type:** add its enum value and scalar calculation in `kinetics.c`,
  map the configuration type and parameters in `generate_c.py`, and add a
  formula case in `tests/test_cpu.c`.
- **Add a Rosenbrock parameter set:** initialize a `chem_rosenbrock_parameters`
  value with its stage count, coefficients, new-function-evaluation flags,
  and estimator order, then pass it to `chem_solver_create`.
- **Replace the step-size controller:** set the `controller` function pointer
  and shared context in `chem_rosenbrock_parameters`. The callback returns a
  suggested factor; the solver applies the configured factor limits and
  rejection handling.
- **Swap the linear solver:** set the `factor` and `solve` callbacks in
  `chem_rosenbrock_parameters`. They receive the symbolic sparse pattern and
  shared context, and can reuse the existing fill graph or build another
  representation.

The coefficient tables and rate data are separate from the mechanism-independent
runtime; generated files contain data only.
