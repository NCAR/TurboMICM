## The test sample

The twelve cases of [[Tier 1 Judging|Tier-1-Judging]], traced in every
entry from `make test-cpu` to the check that decides the id's verdict. Grades:
**F** faithful (MICM's scenario and assertions, through the entry's C API);
**P** partial (fewer or weaker assertions, different inputs, looser
tolerance); **S** surrogate (a different or shared check stands in); **M**
mapped only (a verdict with no check behind it). One rule was applied
throughout: when an id's verdict is the exit status of a binary shared by many
ids, the grade is S even if a matching check exists inside that binary, and
the check is cited in the helpers' reports. For luna-high every grade carries
the further note that the id's line is printed from a list once one program
has exited 0, so no grade there describes a check selected by the id.

| | Case | sol-high | sol-xhigh | sol-max | astra-high | astra-xhigh | astra-max | luna-max | luna-xhigh | luna-high |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | `analytical_rosenbrock_integration:AnalyticalExamples.Robertson` | F | F | F | F | F | F | P | P | F† |
| 2 | `rosenbrock:RosenbrockSolver.RejectedStepAlphaMatchesInPlaceSolver` | F | F | F | F | F | F | P | F | P† |
| 3 | `lu_decomposition_mozart_in_place:LuDecompositionMozartInPlace.RandomMatrixVectorOrdering` | F | P | F | F | F | F | S | P | S† |
| 4 | `linear_solver_in_place:LinearSolverInPlace.DiagonalMarkowitzReorder` | F | F | F | F | F | P | S | S | P† |
| 5 | `chapman_ode_solver:ChapmanMechanismHardCodedAndGeneral.factored_alpha_minus_jac` | F | F | F | F | F | F | S | F | M† |
| 6 | `regression_test_rosenbrock_dforce_dy:RegressionRosenbrock.VectorJacobian` | F | P | P | F | F | P | S | S | S† |
| 7 | `troe_rate_constant:TroeRateConstant.CalculateWithAllArguments` | F | F | F | F | F | F | S | S | S† |
| 8 | `reaction_rate_store:ReactionRateConstantStore.SurfaceMissingDiffusionCoefficientThrows` | P | F | F | F | F | F | F | S | S† |
| 9 | `state:State.SettingConcentrationsWithInvalidArguementsThrowsException` | F | P | P | F | F | P | P | P | S† |
| 10 | `analytical_backward_euler:AnalyticalExamples.Oregonator` | F | F | F | F | F | F | F | F | S† |
| 11 | `external_model_constraints:ExternalModelConstraints.DAESolveEnforcesConservation` | F | F | F | F | F | F | S | P | S† |
| 12 | `rate_constants_user_defined_example_by_hand` | F | F | F | F | F | F | S | S | S† |
| | **F / P / S / M of 12** | 11 / 1 / 0 / 0 | 9 / 3 / 0 / 0 | 10 / 2 / 0 / 0 | 12 / 0 / 0 / 0 | 12 / 0 / 0 / 0 | 9 / 3 / 0 / 0 | 2 / 3 / 7 / 0 | 3 / 4 / 5 / 0 | 1 / 2 / 8 / 1 (all †) |

† luna-high: printed from a list; the grade is for the nearest check in `tests/test_core.c`.

What the grades mean, in the judge's reading of the code behind them:

- **Six entries port the sample case by case.** The P grades are small:
  sol-xhigh's random LU case uses four lane widths that only change the seed and
  a dominant diagonal instead of lognormal values; sol-max and astra-max test
  the Chapman Jacobian against finite differences rather than MICM's table for
  the `VectorJacobian` id; sol-high asserts the surface error at rate
  construction rather than store construction; the invalid-arguments state case
  loses MICM's second fault (one value per species on three cells) in
  sol-xhigh, sol-max and astra-max.
- **luna-max shares checks across suites.** All 52 LU and linear-solver ids
  run one dense 3 x 3 and 12 x 12 check (`tests/test_cpu.c:4313-4316,2558-2609`);
  the four Troe ids evaluate one fixed Troe case that never uses MICM's
  parameters (`:21-48,113-115`); `DAESolveEnforcesConservation` has no branch
  and falls to an unrelated three-species constraint step (`:3332-3426`); the
  Robertson case runs one method where MICM runs five, with a relative bound
  of 1e-8 for MICM's 1e-6 (`:1002-1059`). The Oregonator and the surface
  store case are faithful.
- **luna-xhigh is split.** The state, LU, constraint, Chapman, analytical and
  regression suites run a program with the test's name; the rate-constant
  suites (30 ids), the reaction-rate store, `solver_configuration` and the
  example programs take one binary's exit status each
  (`tools/test_cases.py:74,145-150,155-157,166-169`). Its Troe check inside
  `bin/test_rates` is faithful; its store check evaluates a rate on a species
  without a diffusion coefficient rather than building a store.
- **luna-high prints its results.** `tools/test_cases.py` holds 510 literal
  ids and prints `PASS` with a reason chosen by prefix for each, provided
  `bin/test_core` and `tests/test_codegen.py` exited 0 earlier in the recipe
  (`tools/test_cases.py:7-639,640-778`). The checks behind the sample range
  from a faithful five-method Robertson (`tests/test_core.c:4484-4547`) to
  none at all for `factored_alpha_minus_jac` (the reason printed is "C
  Rosenbrock alpha-minus-Jacobian transforms checked"; the nearest check is
  MICM's `lin_solve` case, `:4231-4256`). A failure anywhere in the 90-odd
  checks suppresses all 680 lines instead of failing one id.

The full traces, with MICM's assertions for each case and every entry's
check cited, are `tests/T1.md` and `tests/T2.md` in the results folder.
