# Test results

`make test-cpu` for astra-xhigh, built from its branch with each compiler; the
`<id> PASS|FAIL|N/A` lines are counted against `reference/micm_cases.txt` (680 ids).

| Entry | Commit | Compiler | PASS | FAIL | N/A | Ids reported exactly once |
|---|---|---|---:|---:|---:|---:|
| astra-xhigh | `fbd365c` | clang | 680 | 0 | 0 | 680 of 680 |
| astra-xhigh | `fbd365c` | gcc-15 | 680 | 0 | 0 | 680 of 680 |

astra-xhigh, like astra-high, reports no case N/A: the five C++-only cases that
the Sol entries mark N/A are reported as PASS. astra-high's tests were verified the
same way on 2026-10-15 (680 PASS with both compilers) and were not rerun here;
for this session it was rebuilt with clang only, for timing.
