# Decisions

- **Mechanism input:** support MICM's tutorial `camp-files` YAML/JSON shape and
  a normalized JSON representation for TS1. PyYAML is pinned in the local
  conda environment; JSON and TS1 generation use only the Python standard
  library otherwise.
- **Third bodies:** encode TS1's `M` as `third_body: true`. It is omitted from
  state columns and contributes one power of ideal-gas molar density for each
  occurrence in a reactant list. Its product occurrences do not create state
  variables.
- **Surface input extension:** store `surface_species` explicitly and map the
  reaction label to `<label>.effective radius [m]` and
  `<label>.particle number concentration [# m-3]`. For tutorial CAMP input the
  gas-phase reactant supplies the surface label; the display name is not used
  as the parameter key. The species record carries its molecular weight and
  diffusion coefficient.
- **Rate units:** tutorial CAMP files use molecule/cm3 falloff coefficients.
  The reader converts Troe `k0_A` by the square of `1e-6 * Avogadro` and
  `kinf_A` by that factor; ternary `kinf_A` uses its reciprocal. Normalized JSON
  such as the translated TS1 table already uses runtime units and passes
  through unchanged.
- **Branched input extension:** expand a `BRANCHED` configuration into two
  processes with the same reactants: an alkoxy product list and a nitrate
  product list. The generated rate parameter records which branch is active.
- **Taylor input extension:** represent the coefficient list as `coefficients`
  in the reaction parameters, with up to 16 coefficients, following MICM's
  rate data structure.
- **Species ordering:** calculate a diagonal Markowitz elimination order from
  the generic reaction graph during C generation. This reduces TS1 symbolic
  fill substantially while preserving name-based CSV access and generic
  mechanism behavior.
- **Cell grouping:** default to one contiguous solver group per OpenMP thread,
  matching the benchmark reference behavior. A fixed group size is available
  as a run-time parameter for experiments.
- **Numerical defaults:** use the MICM three-stage Rosenbrock coefficients,
  controller bounds, `rtol=1e-6`, `atol=1e-3`, and relative internal-step
  defaults. Clamp negative ODE concentrations to zero after each top-level
  solve to match MICM; leave algebraic concentrations to their constraints. The
  TS1 reference comparison verifies the combined behavior.
- **External constraints:** provide linear, mass-action equilibrium, and
  callback-defined constraint types. Callback definitions include explicit
  species dependencies so their Jacobian structure can be reported without
  finite-difference probing.
- **Python environment:** create a prefix environment at `.conda` and keep the
  conda package cache in `.conda-pkgs`; both paths are excluded from Git.
- **Generated solver interface:** keep the generic solver as the numerical
  controller and expose generated rate, RHS/Jacobian, factor, and solve
  functions through an optional mechanism kernel table. This shares the MICM
  Rosenbrock implementation and lets configurations without generated kernels
  use the Tier 1 path unchanged.
- **Generated LU order:** generate direct factorization and triangular-solve
  statements in the same traversal order as the generic no-pivot sparse LU.
  The generated kernel records its species and nonzero counts; the solver
  falls back to generic linear algebra if runtime constraints change the
  pattern.
- **Runtime custom-rate updates:** scenarios scale all user-defined rate
  parameters whose names begin with `j` by the interpolated photolysis factor.
  `chem_solver_update_custom_rates` refreshes only `CHEM_USER_DEFINED` rates,
  since the scenario's temperature, pressure, and other custom parameters stay
  fixed.
- **Scenario grouping and tolerances:** retain run-time grouping rather than a
  build-time cell-group constant. The scenario command accepts explicit group
  size, relative tolerance, and absolute tolerance so it can use the MICM gate
  settings (32 cells, `1e-6`, and `1e-15`).
- **Specialization ablations:** keep `--without-*-kernel` generator options and
  `make bench-ablation` to measure the effect of each generated hot-path
  kernel. These options affect build-time code selection only and do not change
  the production default, which emits all four.
