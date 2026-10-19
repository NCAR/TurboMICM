# Verification of tier3/fable-max

Branch `tier3/fable/max`, commit cb0f8a3 L3r: ROS3 under H211b at 8e-3 / 1e-20 replaces the balanced point as a level-3 point (S15, S15b, S16); the judgment in DECISIONS.md; plots; compilers clang, gcc-15; OPTFLAGS `-O3 -DNDEBUG`; 2026-10-19 08:04.

## Build

| Compiler | make env | make cpu bench |
|---|---|---|
| clang | exit 0 (13 s) | cpu exit 0 (124 s) bench exit 0 (0 s); 0 warnings |
| gcc-15 | exit 0 (12 s) | cpu exit 0 (54 s) bench exit 0 (0 s); 0 warnings |

## MICM's test cases

| Compiler | make test-cpu | PASS | FAIL | N/A | Ids missing | Ids repeated | Unknown ids |
|---|---|---:|---:|---:|---:|---:|---:|
| clang | exit 0 (131 s) | 675 | 0 | 7 | 0 | 0 | 2 |
| gcc-15 | exit 0 (69 s) | 675 | 0 | 7 | 0 | 0 | 2 |

## Tier 1 TS1 gate (10,000 cells, seed 2026, 100 steps of 30 s)

| Compiler | ts1_bench_cpu | elapsed_ms | compare_final.py (1e-10) |
|---|---|---:|---|
| clang | exit 0 | 4623.524 | PASS: PASS: 2110000 values in 10000 cells x 211 columns; largest relative difference 1.04e-13 (BCARYO2VBS, cell 5558) (rtol 1e-10, atol 1e-30) |
| gcc-15 | exit 0 | 7103.614 | PASS: PASS: 2110000 values in 10000 cells x 211 columns; largest relative difference 9.97e-14 (BCARYO2VBS, cell 5558) (rtol 1e-10, atol 1e-30) |

## Operating points on the scenarios (--threads 10, 1,000 cells)

- `L2`: `bin/ts1_scenario --method rodas3 --option rtol=5e-3 --option atol=1e-20 --option quantize=1 --option controller=micm --option safety=0.9 SCENARIO`
- `L3r`: `bin/ts1_scenario --method ros3 --option rtol=8e-3 --option atol=1e-20 --option quantize=1 --option controller=h211b --option safety=0.9 SCENARIO`
- `L3`: `bin/ts1_scenario --method extrap --option order=3 --option sequence=shifted --option rtol=1e-3 --option atol=1e-20 --option quantize=1 --option controller=micm --option safety=0.9 SCENARIO`
- `recommended`: `bin/ts1_scenario --method extrap --option order=4 --option sequence=shifted --option rtol=5e-4 --option atol=1e-22 --option quantize=1 --option controller=micm --option safety=0.9 SCENARIO`
- `accurate`: `bin/ts1_scenario --method extrap --option order=4 --option sequence=shifted --option rtol=2e-4 --option atol=1e-20 --option quantize=1 --option controller=micm --option safety=0.9 SCENARIO`

| Compiler | Point | Scenario | Exit | elapsed_ms | Digits | Worst species | Gate | Level 2 | Level 3 |
|---|---|---|---|---:|---:|---|---|---|---|
| clang | L2 | clean-10min | 0 | 185.727 | 2.585 | ISOP 1.521 | PASS | PASS | FAIL |
| clang | L2 | diurnal-10min | 0 | 208.596 | 2.559 | ISOP 1.392 | PASS | PASS | FAIL |
| clang | L2 | diurnal-30min | 0 | 161.347 | 2.515 | ISOP 1.366 | PASS | PASS | FAIL |
| clang | L2 | diurnal-5min | 0 | 290.321 | 2.734 | NO 1.869 | PASS | PASS | FAIL |
| clang | L2 | polluted-10min | 0 | 257.970 | 2.423 | ISOP 1.294 | PASS | PASS | FAIL |
| clang | L2 | tier1-30s | 0 | 91.627 | 5.091 | C6H5O2 4.402 | PASS | PASS | PASS |
| clang | L3r | clean-10min | 0 | 416.713 | 3.561 | ISOP 2.822 | PASS | PASS | PASS |
| clang | L3r | diurnal-10min | 0 | 435.249 | 3.506 | HYDRALD 2.806 | PASS | PASS | PASS |
| clang | L3r | diurnal-30min | 0 | 371.058 | 3.269 | ISOP 2.373 | PASS | PASS | PASS |
| clang | L3r | diurnal-5min | 0 | 512.405 | 3.471 | HYDRALD 2.761 | PASS | PASS | PASS |
| clang | L3r | polluted-10min | 0 | 530.890 | 3.434 | MVK 2.736 | PASS | PASS | PASS |
| clang | L3r | tier1-30s | 0 | 91.412 | 5.774 | ISOP 4.954 | PASS | PASS | PASS |
| clang | L3 | clean-10min | 0 | 408.710 | 3.612 | ISOP 2.491 | PASS | PASS | PASS |
| clang | L3 | diurnal-10min | 0 | 451.042 | 3.842 | ISOP 2.583 | PASS | PASS | PASS |
| clang | L3 | diurnal-30min | 0 | 388.907 | 3.645 | ISOP 2.515 | PASS | PASS | PASS |
| clang | L3 | diurnal-5min | 0 | 670.018 | 4.121 | ISOP 3.002 | PASS | PASS | PASS |
| clang | L3 | polluted-10min | 0 | 609.803 | 3.946 | ISOP 2.847 | PASS | PASS | PASS |
| clang | L3 | tier1-30s | 0 | 221.573 | 5.866 | ISOPO2VBS 5.316 | PASS | PASS | PASS |
| clang | recommended | clean-10min | 0 | 475.627 | 3.619 | MALO2 2.865 | PASS | PASS | PASS |
| clang | recommended | diurnal-10min | 0 | 490.306 | 3.649 | OCLO 2.692 | PASS | PASS | PASS |
| clang | recommended | diurnal-30min | 0 | 391.247 | 4.237 | ISOP 3.098 | PASS | PASS | PASS |
| clang | recommended | diurnal-5min | 0 | 850.528 | 3.865 | BRO 2.996 | PASS | PASS | PASS |
| clang | recommended | polluted-10min | 0 | 674.143 | 4.52 | ISOP 3.405 | PASS | PASS | PASS |
| clang | recommended | tier1-30s | 0 | 306.010 | 6.525 | MEKO2 5.858 | PASS | PASS | PASS |
| clang | accurate | clean-10min | 0 | 539.643 | 4.809 | ISOP 3.663 | PASS | PASS | PASS |
| clang | accurate | diurnal-10min | 0 | 573.285 | 4.716 | ISOP 3.494 | PASS | PASS | PASS |
| clang | accurate | diurnal-30min | 0 | 499.695 | 4.571 | ISOP 3.535 | PASS | PASS | PASS |
| clang | accurate | diurnal-5min | 0 | 877.388 | 3.865 | BRO 2.996 | PASS | PASS | PASS |
| clang | accurate | polluted-10min | 0 | 720.750 | 4.545 | ISOP 3.41 | PASS | PASS | PASS |
| clang | accurate | tier1-30s | 0 | 313.347 | 6.845 | MEKO2 6.3 | PASS | PASS | PASS |

- `L2`: gate on every scenario: yes; level 2 on every scenario: yes; level 3 on every scenario: NO
- `L3r`: gate on every scenario: yes; level 2 on every scenario: yes; level 3 on every scenario: yes
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

