# Decisions

- The public library is called `chem`. Its implementation is ISO C11, with double
  chemistry and native int/size_t utility arrays. Python is restricted to
  build-time configuration and test tooling.
- Numerical formulas, operation ordering where significant, solver parameters,
  error control, and LU algorithms follow the supplied reference source.
- C++ ownership and template tests will exercise equivalent C allocation, copying,
  destruction, layout, and error-return behavior. Only cases with no meaningful C
  equivalent will be marked N/A, individually and with reasons.
- The supplied experiment prompt remains an untracked local input. It is excluded
  because it contains experiment identity information.
- Reference paths default to sibling directories, overridable in Make variables.
  Builds and runs of reference executables will write only inside this repository.

## Configuration conventions and extensions

- Legacy configuration rate coefficients use molecule/cm3 concentration units.
  The reader converts them to mol/m3 according to reaction order. TS1's supplied
  C++ constants are already in mol/m3; an explicit reaction field
  `concentration units: mol m-3` suppresses conversion. The other accepted value
  is `molecule cm-3` (the legacy default).
- Reactants and products can use the original mappings, including `qty` and
  `yield`. Ordered arrays are also accepted: reactants are species-name strings,
  products are `{species: NAME, yield: VALUE}` objects. This extension preserves
  the exact order and duplicated terms of the C++ TS1 declarations.
- Species may specify `parameterization: air density` for a third body, plus
  `initial concentration [mol m-3]`. The document can supply `initial conditions`
  (temperature and pressure) and `initial rate parameters` (label/value mapping).
  These retain all of TS1 InitState, including its original six photolysis rates
  of 1/s; the benchmark input overrides those values as required.
- `USER_DEFINED` names the general custom-rate type. `PHOTOLYSIS`, `EMISSION`,
  and `FIRST_ORDER_LOSS` retain the legacy convenience forms. The original
  `MUSICA name` key is a data-format label only; no such package is imported.
- Tutorial surface labels differ between file examples (`C surface`) and the
  hand-written programs (`C`). Equivalence tests must map the same physical
  radius/number parameters by their respective labels.
- Species ordering uses the reference's diagonal Markowitz selection and stable
  tie-breaking. Automatic generated kernels use the same symbolic ordering as
  the runtime builder. `preserve_order` supports mechanisms that explicitly need
  their declared variable order.
- `make env` disables conda environment registration and keeps its package cache
  local. This prevents environment creation from updating the global environment
  registry on subsequent runs.

## C array interfaces

- C++ view-category and const-signature tests are mapped to explicit C descriptor
  tags, typed native storage, callback signatures, and rejection of writes through
  const spans. They have meaningful C equivalents and are not N/A.
- Invalid matrix indices return status codes in optimized and debug builds. This
  is stronger than reference checks that are assertion-only under `!NDEBUG`.
- Dense and sparse callbacks share the same grouped span operations. Sparse
  rows/columns map to cached physical offsets while each argument retains its
  own sparsity structure; matrices in a callback need not share a structure.

## Reproducible numerical tests

- The random matrix tests implement the supplied reference build's random engines
  and distributions in test-only C: minstd_rand, MT19937, canonical uniforms, and
  polar normal sampling. This reproduces the reference inputs without a C++
  dependency and preserves the original tight reconstruction tolerances.
- The compatibility version API reports 3.14.0, as declared in the supplied
  source's version.hpp, with a zero tweak component.

The legacy hard-coded Chapman solver is retained as a test-only scalar C oracle,
independent of the generic library's forcing, Jacobian, linear solve and adaptive
loop. It is never linked into the library or benchmark. Its unit cases exercise
the corresponding generic C operations and their original expected quantities;
reciprocal U diagonals are converted only for the legacy factor-storage checks.
Chapman regression cases additionally exercise configuration-generated kernels.

External model descriptions compose extra species, custom parameters, kinetic
patterns and algebraic rows before ordering and symbolic factorization. Contexts
can be copied as plain values, cloned/destroyed by user callbacks, or explicitly
borrowed. Finalization binds each instance after ordering. Kinetic callbacks keep
incoming algebraic rows unchanged, while constraint callbacks supply those rows.
Per-state temporary storage and const model contexts support shared solvers.

## Tutorial output

All eight tutorial entry points are available through `bin/tutorial NAME`. The
case runner executes those same functions against captured MICM output. Every
scientific value, solver statistic, heading and formatting byte must match.
`solver_configuration` also prints measured nanoseconds: this is validated as a
finite nonnegative duration rather than compared to a previous measurement.
`openmp_tutorial` prints three scheduling-dependent thread announcements: their
identities and multiplicity are checked independent of order. The remaining
output is exact. These incidental nondeterministic fields are not chemistry
results; no case is N/A.

The reference vectorized tutorial initializes `solver_state` to Converged but
requires it to be non-Converged before entering its inner solve loop. Its printed
concentrations therefore remain initial values. The port preserves that observed
behavior; the analytical and regression cases exercise actual vector solves.
Sparse string examples store native size_t handles into owned C string storage,
using the same sparse layout indices rather than requiring C++ string objects.

Each of the three tutorial configuration directories is compared against the
independent hand-built mechanism: all rate constants, nonzero forcing/Jacobian
entries, and ten multi-cell integration intervals. Both generated and generic
kernels and both LU algorithms run. Surface labels are mapped as documented
above. The generic configuration reader preserves reaction declaration order;
comparisons match reaction terms rather than assuming equal storage positions.

The benchmark applies the same differential-variable post-solve clamp as the
public solver wrapper after warm-up and each timed solve. Omitting that wrapper
step changed tiny residual concentrations, although the earlier outputs already
passed the specified gate. The final benchmark and tutorial paths now agree on
this behavior.

A short 10,000-cell scheduling/group-size pilot found only small differences
between groups of 16 and 32 and between static and dynamic scheduling. The
32-cell/static default is retained; larger groups degraded performance on this
host. Static, dynamic and guided schedules and explicit chunk sizes are exposed
for other machines. Numerical checks cover all three policies, irregular tails,
and groups larger than the domain. Full measurements and sample counts are in
RESULTS.md; reference single samples are not presented as medians of repeated runs.

## Tier 3 research setup

- Keep the original integration functions and coefficients unchanged. Add a
  separate autonomous-ODE API and a scenario runner; new solvers explicitly
  reject algebraic mechanisms rather than silently integrating constraints.
- Use all six disclosed scenarios. Initial structure scans use the first 64
  cells and matching reference rows, with the supplied scorer. These are
  screening data, not final operating-point evidence. Final validation uses
  every disclosed cell, snapshots, and both compilers.
- Projection of negative candidate concentrations is charged to the local error
  norm before acceptance. Count projected values and summed corrections. Never
  suppress NaN/Inf by projection; reject and retry. Preserve the reference
  wrapper's behavior on the separate compatibility path.
- Snapshot CSV output occurs outside each timed interval. Each interval includes
  photolysis interpolation, rate updates, integration, parallel scheduling, and
  finite/nonnegative validation. There is no warm-up in scenario runs.
- Add linearly implicit Euler extrapolation with sequences 1,2,...,p, p=2..4.
  Freeze the Jacobian during a macro step; extrapolate to zero substep size with
  the Neville recurrence. The preceding diagonal supplies an embedded estimate.
  It is an experimental method, not claimed to be uniformly stiffly accurate.
