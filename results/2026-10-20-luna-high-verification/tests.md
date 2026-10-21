# Test results

`make test-cpu` for luna-high, built from its branch with each compiler; the
`<id> PASS|FAIL|N/A` lines are counted against `reference/micm_cases.txt` (680 ids).

| Entry | Commit | Compiler | PASS | FAIL | N/A | Ids reported exactly once |
|---|---|---|---:|---:|---:|---:|
| luna-high | `d3c2480` | clang | 662 | 0 | 18 | 680 of 680 |
| luna-high | `d3c2480` | gcc-15 | 662 | 0 | 18 | 680 of 680 |

The same 18 cases are N/A in both builds, with the entry's reasons:
- `state:State.CopyConstructor`, "C state views are non-owning, so they have no
  deep-copy constructor"
- `state:State.CopyAssignmentOperator`, "C state views are non-owning, so they
  have no deep-copy assignment"
- `state:State.MoveConstructor`, "C state views are non-owning and have no move
  constructor"
- `state:State.MoveAssignmentOperator`, "C state views are non-owning and have
  no move assignment"
- `solver_builder:SolverBuilder.FinalizeReceivesDenseMatrixPolicy`, "C has no
  templated Finalize overload dispatch or matrix policy types"
- `solver_builder:SolverBuilder.FinalizeSupportsSparseOnlySignature`, the same
- `sparse_matrix_vector_compressed_column_ordering:SparseVectorCompressedColumnMatrix.IncompatibleOrdering`,
  "C uses runtime scalar views and has no compile-time ordering mismatch"
- `sparse_matrix_vector_compressed_column_ordering:SparseVectorCompressedColumnMatrix.IncompatibleVectorOrdering`,
  the same
- `sparse_matrix_vector_compressed_column_ordering:SparseVectorCompressedColumnMatrix.IncompatibleSparseOrdering`,
  the same
- `sparse_matrix_vector_compressed_row_ordering:SparseVectorCompressedRowMatrix.IncompatibleOrdering`,
  the same
- `sparse_matrix_vector_compressed_row_ordering:SparseVectorCompressedRowMatrix.IncompatibleVectorOrdering`,
  the same
- `sparse_matrix_vector_compressed_row_ordering:SparseVectorCompressedRowMatrix.IncompatibleSparseOrdering`,
  the same
- `view_category_sfinae:ViewCategorySFINAE.VectorTypes`, "C has no template
  type-trait or SFINAE behavior"
- `view_category_sfinae:ViewCategorySFINAE.MatrixViews`, the same
- `view_category_sfinae:ViewCategorySFINAE.SparseMatrixViews`, the same
- `view_category_sfinae:ViewCategorySFINAE.MatrixTypesExcluded`, the same
- `view_category_sfinae:ViewCategorySFINAE.RuntimeVectorUse`, the same
- `grouped_views:GroupedView.CategoryTagsAreDistinct`, "C has no distinct
  compile-time view category tags"

Five of them (four `ViewCategorySFINAE` cases and
`GroupedView.CategoryTagsAreDistinct`) are the five the Sol entries mark N/A;
`ViewCategorySFINAE.RuntimeVectorUse` is also N/A in luna-max, and the two
`SolverBuilder.Finalize*` cases in luna-max and luna-xhigh. The other ten, the
four `State` copy and move cases and six incompatible-ordering cases of MICM's
vector sparse matrices, are N/A in no other entry. The Astra entries pass all
680.
