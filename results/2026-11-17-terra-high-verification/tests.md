# Test results

`make test-cpu` for terra-high, built from its branch with each compiler; the
`<id> PASS|FAIL|N/A` lines are counted against `reference/micm_cases.txt` (680 ids).

| Entry | Commit | Compiler | PASS | FAIL | N/A | Ids reported exactly once |
|---|---|---|---:|---:|---:|---:|
| terra-high | `a8a2ea0` | clang | 675 | 0 | 5 | 680 of 680 |
| terra-high | `a8a2ea0` | gcc-15 | 675 | 0 | 5 | 680 of 680 |

The same 5 cases are N/A in both builds, with the entry's reasons:
- `view_category_sfinae:ViewCategorySFINAE.VectorTypes`, `.MatrixViews`,
  `.SparseMatrixViews`, and `.MatrixTypesExcluded`, "C has no template
  concepts or SFINAE; row views use explicit descriptors"
- `grouped_views:GroupedView.CategoryTagsAreDistinct`, "C row and group spans
  use explicit functions rather than template tag dispatch"
