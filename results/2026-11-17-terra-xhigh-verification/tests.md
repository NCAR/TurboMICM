# Test results

`make test-cpu` for terra-xhigh, built from its branch with each compiler; the
`<id> PASS|FAIL|N/A` lines are counted against `reference/micm_cases.txt` (680 ids).

| Entry | Commit | Compiler | PASS | FAIL | N/A | Ids reported exactly once |
|---|---|---|---:|---:|---:|---:|
| terra-xhigh | `288f274` | clang | 676 | 0 | 4 | 680 of 680 |
| terra-xhigh | `288f274` | gcc-15 | 676 | 0 | 4 | 680 of 680 |

The same 4 cases are N/A in both builds, with the entry's reasons:
- `view_category_sfinae:ViewCategorySFINAE.VectorTypes`, "C has no C++ concepts or SFINAE"
- `view_category_sfinae:ViewCategorySFINAE.MatrixViews`, "C has no template view-category traits"
- `view_category_sfinae:ViewCategorySFINAE.SparseMatrixViews`, "C has no template view-category traits"
- `view_category_sfinae:ViewCategorySFINAE.MatrixTypesExcluded`, "C has no C++ concepts or SFINAE"
