# Test results

`make test-cpu` for terra-max, built from its branch with each compiler; the
`<id> PASS|FAIL|N/A` lines are counted against `reference/micm_cases.txt` (680 ids).

| Entry | Commit | Compiler | PASS | FAIL | N/A | Ids reported exactly once |
|---|---|---|---:|---:|---:|---:|
| terra-max | `9071941` | clang | 676 | 0 | 4 | 680 of 680 |
| terra-max | `9071941` | gcc-15 | 676 | 0 | 4 | 680 of 680 |

The same 4 cases are N/A in both builds, with the entry's reason, "C has no
template concepts or SFINAE":
- `view_category_sfinae:ViewCategorySFINAE.VectorTypes`
- `view_category_sfinae:ViewCategorySFINAE.MatrixViews`
- `view_category_sfinae:ViewCategorySFINAE.SparseMatrixViews`
- `view_category_sfinae:ViewCategorySFINAE.MatrixTypesExcluded`
