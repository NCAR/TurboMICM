# Test results

`make test-cpu` for luna-max, built from its branch with each compiler; the
`<id> PASS|FAIL|N/A` lines are counted against `reference/micm_cases.txt` (680 ids).

| Entry | Commit | Compiler | PASS | FAIL | N/A | Ids reported exactly once |
|---|---|---|---:|---:|---:|---:|
| luna-max | `da2d722` | clang | 671 | 0 | 9 | 680 of 680 |
| luna-max | `da2d722` | gcc-15 | 671 | 0 | 9 | 680 of 680 |

The same nine cases are N/A in both builds, with the entry's reasons:
- `grouped_views:GroupedView.CategoryTagsAreDistinct`, "no C equivalent for this
  C++ compile-time or CUDA-only facility"
- `view_category_sfinae:ViewCategorySFINAE.MatrixTypesExcluded`, the same
- `view_category_sfinae:ViewCategorySFINAE.MatrixViews`, the same
- `view_category_sfinae:ViewCategorySFINAE.RuntimeVectorUse`, the same
- `view_category_sfinae:ViewCategorySFINAE.SparseMatrixViews`, the same
- `view_category_sfinae:ViewCategorySFINAE.VectorTypes`, the same
- `process_set:ProcessSetAlgebraicVariables.CudaMatrix`, the same
- `solver_builder:SolverBuilder.FinalizeReceivesDenseMatrixPolicy`, "C has no
  compile-time dense or sparse policy overload to inspect"
- `solver_builder:SolverBuilder.FinalizeSupportsSparseOnlySignature`, the same

The first, and four of the five `ViewCategorySFINAE` cases, are the five the Sol
entries mark N/A; the Astra entries pass all 680.
