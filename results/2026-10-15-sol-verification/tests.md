# Test results

`make test-cpu` for each entry, built from its branch with each compiler; the
`<id> PASS|FAIL|N/A` lines are counted against `reference/micm_cases.txt` (680 ids).

| Entry | Commit | Compiler | PASS | FAIL | N/A | Ids reported exactly once |
|---|---|---|---:|---:|---:|---:|
| sol-high | `1d944f4` | clang | 675 | 0 | 5 | 680 of 680 |
| sol-high | `1d944f4` | gcc-15 | 675 | 0 | 5 | 680 of 680 |
| sol-xhigh | `171f54a` | clang | 675 | 0 | 5 | 680 of 680 |
| sol-xhigh | `171f54a` | gcc-15 | 675 | 0 | 5 | 680 of 680 |
| sol-max | `001fb64` | clang | 675 | 0 | 5 | 680 of 680 |
| sol-max | `001fb64` | gcc-15 | 675 | 0 | 5 | 680 of 680 |

The same five cases are N/A in all six builds:
- `grouped_views:GroupedView.CategoryTagsAreDistinct`
- `view_category_sfinae:ViewCategorySFINAE.MatrixTypesExcluded`
- `view_category_sfinae:ViewCategorySFINAE.MatrixViews`
- `view_category_sfinae:ViewCategorySFINAE.SparseMatrixViews`
- `view_category_sfinae:ViewCategorySFINAE.VectorTypes`
