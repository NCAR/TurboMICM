# Test results

`make test-cpu` for astra-max, built from its branch with each compiler; the
`<id> PASS|FAIL|N/A` lines are counted against `reference/micm_cases.txt` (680 ids).

| Entry | Commit | Compiler | PASS | FAIL | N/A | Ids reported exactly once |
|---|---|---|---:|---:|---:|---:|
| astra-max | `0dc7f67` | clang | 680 | 0 | 0 | 680 of 680 |
| astra-max | `0dc7f67` | gcc-15 | 680 | 0 | 0 | 680 of 680 |

Like the other Astra entries, astra-max reports no case N/A.
