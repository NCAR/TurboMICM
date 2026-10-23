# Design

## Generic library

The public C interface is in `include/kinetics.h`. A mechanism is immutable
table data: species, flattened reactants and products, yields, rate parameters,
and custom parameter defaults. The generic solver caches environment-dependent
rate constants when it is created. Concentrations and conditions belong to the
caller and use cell-major concentration rows.

The rate evaluator handles Arrhenius, Troe, ternary chemical activation,
JPL-19 ternary activation, branched, tunneling, Taylor-series, reversible,
user-defined, surface, and C callback rates. Third-body reactants contribute
the corresponding power of ideal-gas molar density but have no state column.
Surface rates use the named species' molecular weight and diffusion
coefficient, plus custom particle radius and number concentration.

The generic chemistry solver implements MICM's three-stage Rosenbrock method,
embedded error estimate, step-size limits, acceptance logic, and default
tolerances. It builds the sparse Jacobian graph and symbolic no-pivot LU fill
at solver creation. The public API also supports alternate Rosenbrock
parameter sets, linear constraints, mass-action equilibrium constraints, and
caller-defined residuals with Jacobians.

## Generated specialization

`tools/generate_c.py --specialized` reads the same configuration format as the
generic reader and emits C functions with species indexes, reaction parameters,
and matrix positions fixed at generation time:

- Rate formulas are emitted per reaction with parameters and third-body order
  as literals. Temperature, pressure, air density, and custom rate values stay
  runtime inputs. Cached user-defined rates have a separate update kernel for
  scenarios that change photolysis values between steps.
- Reaction flux and negative-Jacobian updates are emitted in reaction order
  with reactant and product indexes, yields, and fill positions fixed. The
  generated source has no run-time reaction or stoichiometry loops.
- A diagonal Markowitz ordering is selected at generation time. The generator
  computes the symbolic no-pivot fill graph and emits direct positions for
  every elimination update.
- LU factorization and forward/back substitution are emitted as fixed
  operations using the generic solver's row and column traversal order. The
  generated factorization preserves the same `fma` updates and pivot checks.

The generated kernel table is optional on `chem_mechanism`. The generic solver
uses a generated kernel only when the run-time sparse pattern has the expected
dimensions. Constraint installation changes the pattern and therefore uses
the generic constraint-aware evaluator and linear algebra. The numerical
controller, error norm, stage coefficients, per-cell state, and step-size
history remain shared with the generic solver.

`--without-rate-kernel`, `--without-chemistry-kernel`,
`--without-factor-kernel`, and `--without-solve-kernel` are generator ablation
options for benchmarking individual kernels. Production Make targets generate
all four kernels. Generated C is never edited by hand.

## Configuration to C

The generator follows a configuration's `camp-files`, reads its species and
mechanism objects, resolves names to indexes, expands integer reactant
quantities, and emits static C data. It rejects unknown species, unsupported
rate types, noninteger reactant quantities, and surface reactions without
required properties. YAML input uses PyYAML; JSON generation otherwise needs
only the standard library.

The TS1 translator reads species, reaction builders, and initial-state setters
from the read-only TS1 header. It writes the normalized JSON consumed by the
ordinary generator. TS1's rate constants, species ordering, reaction graph,
and symbolic fill are not hard-coded in the solver implementation.

## Cell groups and scenarios

Cells are split into contiguous groups. By default, the group count equals the
OpenMP maximum thread count. `--threads N` controls that count and
`--group-size G` sets a fixed number of cells per group. A static OpenMP loop
assigns groups to threads. Each group owns its adaptive step-size history and
statistics; statistics are summed after the groups finish.

The scenario driver reads relative paths from the scenario file, initializes
the state from the cell CSV, linearly interpolates photolysis factor rows, and
updates all custom parameter names beginning with `j` at the start of each
step. It has no warm-up interval. The measured interval includes photolysis
parameter updates and chemistry integration, not CSV output.

## Runtime and build parameters

| Parameter | Default | Override |
|---|---:|---|
| Compiler optimization flags | `-O3 -DNDEBUG` | `OPTFLAGS=...` |
| OpenMP thread count | runtime maximum | `--threads N` |
| Cell group size | one group per thread | `--group-size G` |
| Rosenbrock method | MICM three-stage default | C solver parameters |
| Relative tolerance | `1e-6` | C API or scenario `--rtol` |
| Absolute tolerance | `1e-3` per species | C API or scenario `--atol` |
| Initial internal step | `1e-6 * DT` | C solver parameters |
| Minimum internal step | `1e-15 * DT` | C solver parameters |

No compiler target, floating-point, LTO, profile-guided, or additional
optimization flags are set by the Makefile. Generated and handwritten C use
the same `CC` and `OPTFLAGS`.

## Regeneration and extension

`make bench` regenerates the TS1 generic and specialized C files when the input
or generator changes. `make test-specialized` generates and compiles separate
models for TS1 and the three Tier 1 configuration examples, then compares each
solver with the generic library.

To add a supported reaction rate, add its configuration mapping and emitted
formula in `tools/generate_c.py`, then test it against the generic runtime.
The rate tables and Rosenbrock controller remain separate: mechanism-specific
chemistry is generated, while tolerances and integration controls remain
runtime settings.
