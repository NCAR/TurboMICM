# Verification of tier2/luna/xhigh

Branch `tier2/luna/xhigh`, commit 1ad13ab Add mechanism-specialized Tier 2 solver; compilers clang, gcc-15; OPTFLAGS `-O3 -DNDEBUG`; 18 threads; 2026-10-22 23:06.

## Build (fresh worktree)

| Compiler | make env | make cpu | make bench | Clean build | Warnings | Generated C |
|---|---|---|---|---:|---:|---|
| clang | exit 0 (12 s) | exit 0 (66.8 s) | exit 0 (0.9 s) | 67.7 s | 0 | 1 files, 48,068 lines, 2.40 MB |
| gcc-15 | exit 0 (10 s) | exit 0 (62.8 s) | exit 0 (0.5 s) | 63.3 s | 0 | 1 files, 48,068 lines, 2.40 MB |

clang: generated C: `build/generated/ts1.c`

gcc-15: generated C: `build/generated/ts1.c`

## MICM's test cases (the Tier 1 gate)

| Compiler | make test-cpu | PASS | FAIL | N/A | Ids missing | Ids repeated | Unknown ids |
|---|---|---:|---:|---:|---:|---:|---:|
| clang | exit 0 (13 s) | 662 | 0 | 18 | 0 | 0 | 0 |
| gcc-15 | exit 0 (17 s) | 662 | 0 | 18 | 0 | 0 | 0 |

## Tier 1 TS1 gate: bin/ts1_bench_cpu (10,000 cells, seed 2026, 100 steps of 30 s)

| Compiler | Exit | elapsed_ms | compare_final.py (1e-10) |
|---|---|---:|---|
| clang | 0 | 10674.357 | PASS: PASS: 2110000 values in 10000 cells x 211 columns; largest relative difference 1.77e-13 (BCARY, cell 6848) (rtol 1e-10, atol 1e-30) |
| gcc-15 | 0 | 10718.003 | PASS: PASS: 2110000 values in 10000 cells x 211 columns; largest relative difference 2.33e-13 (BCARY, cell 4179) (rtol 1e-10, atol 1e-30) |

## Specialized TS1 gate: bin/ts1_bench_specialized (10,000 cells, 100 steps of 30 s)

| Compiler | Input | Exit | elapsed_ms | Internal steps | compare_final.py (1e-10) |
|---|---|---|---:|---:|---|
| clang | seed 2026 | 0 | 10096.074 | 16200 | PASS: PASS: 2110000 values in 10000 cells x 211 columns; largest relative difference 2.15e-13 (BCARY, cell 4179) (rtol 1e-10, atol 1e-30) |
| clang | input 1 | 0 | 10056.181 | 16200 | PASS: PASS: 2110000 values in 10000 cells x 211 columns; largest relative difference 1.61e-13 (BCARYO2VBS, cell 9174) (rtol 1e-10, atol 1e-30) |
| gcc-15 | seed 2026 | 0 | 10009.745 | 16200 | PASS: PASS: 2110000 values in 10000 cells x 211 columns; largest relative difference 2.42e-13 (BCARY, cell 6848) (rtol 1e-10, atol 1e-30) |
| gcc-15 | input 1 | 0 | 10031.206 | 16200 | PASS: PASS: 2110000 values in 10000 cells x 211 columns; largest relative difference 1.91e-13 (BCARYO2VBS, cell 2723) (rtol 1e-10, atol 1e-30) |

## The generator's tests: make test-specialized

| Compiler | Exit | Time | PASS lines | FAIL lines | Last lines |
|---|---|---:|---:|---:|---|
| clang | 0 | 2 s | 3 | 0 | PASS specialized equivalence: rate_constants_no_user_defined (7 species, 7 reactions, 9 internal steps) / bin/test_specialized_user / PASS specialized equivalence: rate_constants_user_defined (7 speci |
| gcc-15 | 0 | 2 s | 3 | 0 | PASS specialized equivalence: rate_constants_no_user_defined (7 species, 7 reactions, 9 internal steps) / bin/test_specialized_user / PASS specialized equivalence: rate_constants_user_defined (7 speci |

## Day-night scenarios (--threads 18, --rtol 1e-6 --atol 1e-15, 32-cell groups)

Agreement is score.py against MICM's result at the same settings, with `--level 10`; accuracy is score.py against the reference solution, for scale.

| Compiler | Scenario | Exit | elapsed_ms | Internal steps (MICM's) | Agreement with MICM | Worst species | Gate | Level 10 | Accuracy (MICM's) |
|---|---|---|---:|---|---:|---|---|---|---|
| clang | arctic-summer-10min | 0 | 1392.069 | 103211 (103211) | 14.218 | ISOP 13.714 | PASS | PASS | 5.045 (5.04) |
| clang | night-start-5min | 0 | 1838.413 | 132855 (132855) | 14.306 | ISOP 13.83 | PASS | PASS | 4.729 (4.73) |
| clang | south-48h-30min | 0 | 1257.994 | 95163 (95163) | 14.229 | ISOP 13.754 | PASS | PASS | 4.966 (4.97) |
| clang | spring45-20min | 0 | 800.613 | 60435 (60435) | 14.253 | ISOP 13.655 | PASS | PASS | 4.842 (4.84) |
| clang | tropics-polluted-15min | 0 | 995.756 | 75293 (75293) | 14.249 | ISOP 13.491 | PASS | PASS | 5.1 (5.10) |
| clang | clean-10min | 0 | 1135.816 | 84273 (84273) | 14.19 | ISOP 13.6 | PASS | PASS | 4.907 (4.91) |
| clang | diurnal-10min | 0 | 1219.025 | 90980 (90980) | 14.184 | ISOP 13.637 | PASS | PASS | 5.052 (5.05) |
| clang | diurnal-30min | 0 | 651.712 | 48619 (48619) | 14.258 | ISOP 13.689 | PASS | PASS | 4.748 (4.75) |
| clang | diurnal-5min | 0 | 1993.976 | 145891 (145891) | 14.151 | ISOP 13.587 | PASS | PASS | 5.214 (5.21) |
| clang | polluted-10min | 0 | 1408.352 | 106393 (106393) | 14.012 | ISOP 13.415 | PASS | PASS | 5.398 (5.40) |
| gcc-15 | arctic-summer-10min | 0 | 1640.001 | 103211 (103211) | 14.23 | ISOP 13.736 | PASS | PASS | 5.045 (5.04) |
| gcc-15 | night-start-5min | 0 | 2175.533 | 132855 (132855) | 14.307 | ISOP 13.836 | PASS | PASS | 4.729 (4.73) |
| gcc-15 | south-48h-30min | 0 | 1479.011 | 95163 (95163) | 14.237 | ISOPAO2 13.623 | PASS | PASS | 4.966 (4.97) |
| gcc-15 | spring45-20min | 0 | 951.141 | 60435 (60435) | 14.251 | ISOP 13.668 | PASS | PASS | 4.842 (4.84) |
| gcc-15 | tropics-polluted-15min | 0 | 1186.985 | 75293 (75293) | 14.248 | ISOP 13.491 | PASS | PASS | 5.1 (5.10) |
| gcc-15 | clean-10min | 0 | 1360.551 | 84273 (84273) | 14.202 | ISOP 13.616 | PASS | PASS | 4.907 (4.91) |
| gcc-15 | diurnal-10min | 0 | 1451.324 | 90980 (90980) | 14.181 | ISOP 13.623 | PASS | PASS | 5.052 (5.05) |
| gcc-15 | diurnal-30min | 0 | 757.606 | 48619 (48619) | 14.255 | ISOP 13.674 | PASS | PASS | 4.748 (4.75) |
| gcc-15 | diurnal-5min | 0 | 2355.316 | 145891 (145891) | 14.146 | ISOP 13.598 | PASS | PASS | 5.214 (5.21) |
| gcc-15 | polluted-10min | 0 | 1679.170 | 106393 (106393) | 13.991 | ISOPNO3 13.306 | PASS | PASS | 5.398 (5.40) |

- clang: level 10 and the gate on every scenario: yes (10 scenarios)
- gcc-15: level 10 and the gate on every scenario: yes (10 scenarios)

## Build audit

`make -n -B cpu bench CC=clang`: exit 0; flags used: `-DCHEM_SCENARIO -DCHEM_USE_SPECIALIZED -DNDEBUG -Iinclude -L~/homebrew/opt/libomp/lib -O3 -Wall -Wextra -Wpedantic -Xpreprocessor -c -fopenmp -isystem -lm -lomp -p -std=c11`

Flags outside the rules (review): none

- `bin/test_aerosol_equivalents` links: ~/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/test_analytical_additional` links: ~/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/test_analytical_simple` links: ~/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/test_backward_euler` links: ~/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/test_chapman` links: ~/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/test_constraints` links: ~/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/test_example_no_user` links: ~/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/test_example_user` links: ~/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/test_external_constraints` links: ~/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/test_grouped_views` links: ~/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/test_integration_equivalents` links: ~/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/test_jacobian_verification` links: ~/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/test_lu` links: ~/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/test_matrix` links: ~/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/test_objects` links: ~/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/test_processes` links: ~/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/test_rates` links: ~/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/test_regression` links: ~/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/test_rosenbrock_utilities` links: ~/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/test_solver_configuration` links: ~/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/test_sparse_matrix` links: ~/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/test_specialized_no_user` links: ~/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/test_specialized_ts1` links: ~/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/test_specialized_user` links: ~/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/test_state` links: ~/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/test_tutorial_programs` links: ~/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/ts1_bench_cpu` links: ~/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/ts1_bench_specialized` links: ~/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/ts1_scenario_specialized` links: ~/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib

C sources scanned: 43 (tracked and generated, outside conda); lines matching a banned construct or starting a process: 3 (all in `logs/banned.txt`)

    tests/test_grouped_views.c:164: const size_t vector_size = blocks * columns;
    tests/test_grouped_views.c:169: if (chem_group_view_init(dense, vector_size, group_block * columns, lanes, columns, &output) != 0 ||
    tests/test_grouped_views.c:170: chem_const_group_view_init(dense, vector_size, group_block * columns + 1u, lanes, columns, &dense_input) != 0 ||

