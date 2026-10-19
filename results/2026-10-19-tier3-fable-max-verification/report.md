# Verification of tier3/fable-max

Branch `tier3/fable/max`, commit ec34c99 Final: build times with both compilers, two gcc warnings fixed, the flags rule clears objects on a compiler change; compilers clang, gcc-15; OPTFLAGS `-O3 -DNDEBUG`; 2026-10-19 07:34.

## Build

| Compiler | make env | make cpu bench |
|---|---|---|
| clang | exit 0 (14 s) | cpu exit 0 (127 s) bench exit 0 (0 s); 0 warnings |
| gcc-15 | exit 0 (23 s) | cpu exit 0 (54 s) bench exit 0 (0 s); 0 warnings |

## MICM's test cases

| Compiler | make test-cpu | PASS | FAIL | N/A | Ids missing | Ids repeated | Unknown ids |
|---|---|---:|---:|---:|---:|---:|---:|
| clang | exit 0 (134 s) | 675 | 0 | 7 | 0 | 0 | 2 |
| gcc-15 | exit 0 (69 s) | 675 | 0 | 7 | 0 | 0 | 2 |

## Tier 1 TS1 gate (10,000 cells, seed 2026, 100 steps of 30 s)

| Compiler | ts1_bench_cpu | elapsed_ms | compare_final.py (1e-10) |
|---|---|---:|---|
| clang | exit 0 | 4596.604 | PASS: PASS: 2110000 values in 10000 cells x 211 columns; largest relative difference 1.04e-13 (BCARYO2VBS, cell 5558) (rtol 1e-10, atol 1e-30) |
| gcc-15 | exit 0 | 7089.028 | PASS: PASS: 2110000 values in 10000 cells x 211 columns; largest relative difference 9.97e-14 (BCARYO2VBS, cell 5558) (rtol 1e-10, atol 1e-30) |

## Operating points on the scenarios (--threads 10, 1,000 cells)

- `L2`: `bin/ts1_scenario --method rodas3 --option rtol=5e-3 --option atol=1e-20 --option quantize=1 --option controller=micm --option safety=0.9 SCENARIO`
- `balanced`: `bin/ts1_scenario --method rodas3 --option rtol=2e-3 --option atol=1e-20 --option quantize=1 --option controller=h211b --option safety=0.9 SCENARIO`
- `L3`: `bin/ts1_scenario --method extrap --option order=3 --option sequence=shifted --option rtol=1e-3 --option atol=1e-20 --option quantize=1 --option controller=micm --option safety=0.9 SCENARIO`
- `recommended`: `bin/ts1_scenario --method extrap --option order=4 --option sequence=shifted --option rtol=5e-4 --option atol=1e-22 --option quantize=1 --option controller=micm --option safety=0.9 SCENARIO`
- `accurate`: `bin/ts1_scenario --method extrap --option order=4 --option sequence=shifted --option rtol=2e-4 --option atol=1e-20 --option quantize=1 --option controller=micm --option safety=0.9 SCENARIO`

| Compiler | Point | Scenario | Exit | elapsed_ms | Digits | Worst species | Gate | Level 2 | Level 3 |
|---|---|---|---|---:|---:|---|---|---|---|
| clang | L2 | clean-10min | 0 | 184.102 | 2.585 | ISOP 1.521 | PASS | PASS | FAIL |
| clang | L2 | diurnal-10min | 0 | 207.903 | 2.559 | ISOP 1.392 | PASS | PASS | FAIL |
| clang | L2 | diurnal-30min | 0 | 165.761 | 2.515 | ISOP 1.366 | PASS | PASS | FAIL |
| clang | L2 | diurnal-5min | 0 | 285.237 | 2.734 | NO 1.869 | PASS | PASS | FAIL |
| clang | L2 | polluted-10min | 0 | 240.693 | 2.423 | ISOP 1.294 | PASS | PASS | FAIL |
| clang | L2 | tier1-30s | 0 | 91.050 | 5.091 | C6H5O2 4.402 | PASS | PASS | PASS |
| clang | balanced | clean-10min | 0 | 364.767 | 3.26 | ISOP 2.256 | PASS | PASS | PASS |
| clang | balanced | diurnal-10min | 0 | 409.581 | 3.248 | ISOP 2.169 | PASS | PASS | PASS |
| clang | balanced | diurnal-30min | 0 | 336.511 | 3.302 | ISOP 2.32 | PASS | PASS | PASS |
| clang | balanced | diurnal-5min | 0 | 492.270 | 3.749 | ISOP 2.887 | PASS | PASS | PASS |
| clang | balanced | polluted-10min | 0 | 488.183 | 3.175 | ISOP 2.175 | PASS | PASS | PASS |
| clang | balanced | tier1-30s | 0 | 129.559 | 5.803 | ISOP 5.019 | PASS | PASS | PASS |
| clang | L3 | clean-10min | 0 | 411.593 | 3.612 | ISOP 2.491 | PASS | PASS | PASS |
| clang | L3 | diurnal-10min | 0 | 458.450 | 3.842 | ISOP 2.583 | PASS | PASS | PASS |
| clang | L3 | diurnal-30min | 0 | 389.492 | 3.645 | ISOP 2.515 | PASS | PASS | PASS |
| clang | L3 | diurnal-5min | 0 | 660.354 | 4.121 | ISOP 3.002 | PASS | PASS | PASS |
| clang | L3 | polluted-10min | 0 | 542.308 | 3.946 | ISOP 2.847 | PASS | PASS | PASS |
| clang | L3 | tier1-30s | 0 | 218.516 | 5.866 | ISOPO2VBS 5.316 | PASS | PASS | PASS |
| clang | recommended | clean-10min | 0 | 468.155 | 3.619 | MALO2 2.865 | PASS | PASS | PASS |
| clang | recommended | diurnal-10min | 0 | 485.433 | 3.649 | OCLO 2.692 | PASS | PASS | PASS |
| clang | recommended | diurnal-30min | 0 | 375.407 | 4.237 | ISOP 3.098 | PASS | PASS | PASS |
| clang | recommended | diurnal-5min | 0 | 840.698 | 3.865 | BRO 2.996 | PASS | PASS | PASS |
| clang | recommended | polluted-10min | 0 | 621.922 | 4.52 | ISOP 3.405 | PASS | PASS | PASS |
| clang | recommended | tier1-30s | 0 | 305.723 | 6.525 | MEKO2 5.858 | PASS | PASS | PASS |
| clang | accurate | clean-10min | 0 | 541.725 | 4.809 | ISOP 3.663 | PASS | PASS | PASS |
| clang | accurate | diurnal-10min | 0 | 570.458 | 4.716 | ISOP 3.494 | PASS | PASS | PASS |
| clang | accurate | diurnal-30min | 0 | 485.839 | 4.571 | ISOP 3.535 | PASS | PASS | PASS |
| clang | accurate | diurnal-5min | 0 | 872.102 | 3.865 | BRO 2.996 | PASS | PASS | PASS |
| clang | accurate | polluted-10min | 0 | 717.350 | 4.545 | ISOP 3.41 | PASS | PASS | PASS |
| clang | accurate | tier1-30s | 0 | 309.744 | 6.845 | MEKO2 6.3 | PASS | PASS | PASS |

- `L2`: gate on every scenario: yes; level 2 on every scenario: yes; level 3 on every scenario: NO
- `balanced`: gate on every scenario: yes; level 2 on every scenario: yes; level 3 on every scenario: yes
- `L3`: gate on every scenario: yes; level 2 on every scenario: yes; level 3 on every scenario: yes
- `recommended`: gate on every scenario: yes; level 2 on every scenario: yes; level 3 on every scenario: yes
- `accurate`: gate on every scenario: yes; level 2 on every scenario: yes; level 3 on every scenario: yes

## Build audit

`make -n -B cpu bench CC=clang`: exit 0; flags used: `--chunk --file --mode --name --ordering --prefix --width --yes -DNDEBUG -DTMC_DEFAULT_GROUP_SIZE=32 -DTMC_WIDTH=2 -Igenerated -Iinclude -L/opt/homebrew/opt/libomp/lib -L/opt/homebrew/opt/libomp/lib -lomp -lm -O3 -Wall -Wextra -Wpedantic -Xpreprocessor -c -f -fopenmp -isystem -lm -lomp -p -s -std=c11 -std=c11 -Wall -Wextra -Wpedantic -O3 -DNDEBUG -Xpreprocessor -fopenmp -isystem /opt/homebrew/opt/libomp/include -Iinclude -DTMC_WIDTH=2 -DTMC_DEFAULT_GROUP_SIZE=32`

Flags outside the rules (review): none

- `bin/test_cpu` links: /opt/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/ts1_bench_cpu` links: /opt/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/ts1_scenario` links: /opt/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib

C sources scanned: 116 (tracked and generated, outside conda); lines matching a banned construct: 18 (all in `logs/banned.txt`)

    tests/sparse_policy.c:1415: static int vector_size_mismatch(tmc_sparse_ordering ordering, int L, int n_blocks, int size)
    tests/sparse_policy.c:1425: return vector_size_mismatch(ordering, L, 3, 2);
    tests/sparse_policy.c:1430: return vector_size_mismatch(ordering, L, 3, 4);
    tests/sparse_policy.c:1435: return vector_size_mismatch(ordering, L, 3, 0);
    tests/sparse_policy.c:1449: return vector_size_mismatch(ordering, L, 0, 3);
    tests/sparse_policy.c:1581: int sparse_policy_vector_size_matches_one_sparse_matrix_only(tmc_sparse_ordering ordering, int L)
    tests/test_sparse_matrix_vector_compressed_row_ordering.c:226: WIDTH_CASE(vector_size_matches_one_sparse_matrix_only, sparse_policy_vector_size_matches_one_sparse_matrix_only, 1)
    tests/test_sparse_matrix_vector_compressed_row_ordering.c:279: TMC_CASE(ID("VectorSizeMatchesOneSparseMatrixOnly"), vector_size_matches_one_sparse_matrix_only),
    tests/test_sparse_matrix_vector_compressed_column_ordering.c:226: WIDTH_CASE(vector_size_matches_one_sparse_matrix_only, sparse_policy_vector_size_matches_one_sparse_matrix_only, 1)
    tests/test_sparse_matrix_vector_compressed_column_ordering.c:279: TMC_CASE(ID("VectorSizeMatchesOneSparseMatrixOnly"), vector_size_matches_one_sparse_matrix_only),
    tests/test_state.c:574: TMC_TEST(setting_custom_rate_parameter_with_invalid_vector_size_throws_exception)
    tests/test_state.c:765: TMC_CASE("state:State.SettingCustomRateParameterWithInvalidVectorSizeThrowsException", setting_custom_rate_parameter_with_invalid_vector_siz
    tests/sparse_policy.h:57: int sparse_policy_vector_size_matches_one_sparse_matrix_only(tmc_sparse_ordering ordering, int L);
    tests/test_sparse_matrix_standard_compressed_column_ordering.c:148: STANDARD_CASE(vector_size_matches_one_sparse_matrix_only, sparse_policy_vector_size_matches_one_sparse_matrix_only)
    tests/test_sparse_matrix_standard_compressed_column_ordering.c:197: TMC_CASE(ID("VectorSizeMatchesOneSparseMatrixOnly"), vector_size_matches_one_sparse_matrix_only),
    tests/test_sparse_matrix_standard_compressed_row_ordering.c:148: STANDARD_CASE(vector_size_matches_one_sparse_matrix_only, sparse_policy_vector_size_matches_one_sparse_matrix_only)
    tests/test_sparse_matrix_standard_compressed_row_ordering.c:197: TMC_CASE(ID("VectorSizeMatchesOneSparseMatrixOnly"), vector_size_matches_one_sparse_matrix_only),
    tests/test_system.c:84: /* MICM requests "name [units]" as long double, a type no property can hold,

