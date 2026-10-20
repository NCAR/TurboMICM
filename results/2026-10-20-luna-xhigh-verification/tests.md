# Test results

`make test-cpu` for luna-xhigh, built from its branch with each compiler; the
`<id> PASS|FAIL|N/A` lines are counted against `reference/micm_cases.txt` (680 ids).

| Entry | Commit | Compiler | PASS | FAIL | N/A | Ids reported exactly once |
|---|---|---|---:|---:|---:|---:|
| luna-xhigh | `ada55d5` | clang | 662 | 0 | 18 | 680 of 680 |
| luna-xhigh | `ada55d5` | gcc-15 | 662 | 0 | 18 | 680 of 680 |

The same 18 cases are N/A in both builds, with the entry's reasons:
- `reaction_rate_store:ReactionRateConstantStore.LambdaEntriesRcIndex`, "lambda
  callback entries are not part of the C named-rate descriptor"
- `solver_builder:SolverBuilder.ThrowsUnusedSpecies`, "C has no solver builder
  or compile-time matrix policy dispatch"
- `solver_builder:SolverBuilder.FinalizeReceivesDenseMatrixPolicy`, the same
- `solver_builder:SolverBuilder.FinalizeSupportsSparseOnlySignature`, the same
- `view_category_sfinae:ViewCategorySFINAE.VectorTypes`, "C has no template
  categories or SFINAE type traits"
- `view_category_sfinae:ViewCategorySFINAE.MatrixViews`, the same
- `view_category_sfinae:ViewCategorySFINAE.SparseMatrixViews`, the same
- `view_category_sfinae:ViewCategorySFINAE.MatrixTypesExcluded`, the same
- `grouped_views:GroupedView.CategoryTagsAreDistinct`, "C uses explicit
  pointer-and-stride views instead of compile-time C++ view categories"
- `version:Version.FullVersion`, "release metadata is not part of the numerical
  C library interface"
- `version:Version.VersionMajor`, the same
- `version:Version.VersionMinor`, the same
- `version:Version.VersionPatch`, the same
- `version:Version.VersionTweak`, the same
- `external_model_constraints:ExternalModelConstraints.ProcessJacobianElementInAlgebraicRowSurvivesFiltering`,
  "the C constrained integrator uses dense matrices, so sparse algebraic-row
  filtering has no equivalent"
- `external_model_constraints:ExternalModelConstraints.ProcessJacobianElementInBuiltInAlgebraicRowSurvivesFiltering`,
  the same
- `external_model_constraints:ExternalModelConstraints.ConstraintActiveWhenProcessOnlyModelAddedFirst`,
  "C passes one flat mechanism and constraint set; it has no external-model
  registry whose constraint mask depends on model insertion order"
- `external_model_constraints:ExternalModelConstraints.BuilderCanBuildMultipleSolvers`,
  "C passes immutable mechanism and constraint descriptors directly instead of
  exposing a reusable solver builder"

Five of them (the four `ViewCategorySFINAE` cases and
`GroupedView.CategoryTagsAreDistinct`) are the five the Sol entries mark N/A;
the two `SolverBuilder.Finalize*` cases are also N/A in luna-max; the other 11
(the five `Version` cases, the lambda rate entry, `ThrowsUnusedSpecies`, and
four external-model constraint cases) are N/A in no other entry. The Astra
entries pass all 680.
