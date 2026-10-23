# Verification of tier2/luna/high

Branch `tier2/luna/high`, commit 2142282 Add generated TS1 solver specialization; compilers clang, gcc-15; OPTFLAGS `-O3 -DNDEBUG`; 18 threads; 2026-10-23 00:19.

## Build (fresh worktree)

| Compiler | make env | make cpu | make bench | Clean build | Warnings | Generated C |
|---|---|---|---|---:|---:|---|
| clang | exit 0 (12 s) | exit 0 (1.8 s) | exit 0 (10.8 s) | 12.6 s | 0 | 1 files, 31,335 lines, 1.12 MB |
| gcc-15 | exit 0 (10 s) | exit 0 (4.2 s) | exit 0 (26.8 s) | 31.0 s | 0 | 1 files, 31,335 lines, 1.12 MB |

clang: generated C: `generated/ts1.c`

gcc-15: generated C: `generated/ts1.c`

## MICM's test cases (the Tier 1 gate)

| Compiler | make test-cpu | PASS | FAIL | N/A | Ids missing | Ids repeated | Unknown ids |
|---|---|---:|---:|---:|---:|---:|---:|
| clang | exit 0 (4 s) | 662 | 0 | 18 | 0 | 0 | 0 |
| gcc-15 | exit 0 (6 s) | 662 | 0 | 18 | 0 | 0 | 0 |

## Tier 1 TS1 gate: bin/ts1_bench_cpu (10,000 cells, seed 2026, 100 steps of 30 s)

| Compiler | Exit | elapsed_ms | compare_final.py (1e-10) |
|---|---|---:|---|
| clang | 0 | 20580.240 | PASS: PASS: 2110000 values in 10000 cells x 211 columns; largest relative difference 5.34e-14 (BCARYO2VBS, cell 9034) (rtol 1e-10, atol 1e-30) |
| gcc-15 | 0 | 25785.932 | PASS: PASS: 2110000 values in 10000 cells x 211 columns; largest relative difference 1.43e-13 (BCARYO2VBS, cell 4054) (rtol 1e-10, atol 1e-30) |

## Specialized TS1 gate: bin/ts1_bench_specialized (10,000 cells, 100 steps of 30 s)

| Compiler | Input | Exit | elapsed_ms | Internal steps | compare_final.py (1e-10) |
|---|---|---|---:|---:|---|
| clang | seed 2026 | 0 | 7562.616 | 16200 | PASS: PASS: 2110000 values in 10000 cells x 211 columns; largest relative difference 5.34e-14 (BCARYO2VBS, cell 9034) (rtol 1e-10, atol 1e-30) |
| clang | input 1 | 0 | 7412.077 | 16200 | PASS: PASS: 2110000 values in 10000 cells x 211 columns; largest relative difference 4.1e-14 (BCARYO2VBS, cell 7528) (rtol 1e-10, atol 1e-30) |
| gcc-15 | seed 2026 | 0 | 9129.558 | 16200 | PASS: PASS: 2110000 values in 10000 cells x 211 columns; largest relative difference 1.43e-13 (BCARYO2VBS, cell 4054) (rtol 1e-10, atol 1e-30) |
| gcc-15 | input 1 | 0 | 9305.486 | 16200 | PASS: PASS: 2110000 values in 10000 cells x 211 columns; largest relative difference 1.52e-13 (BCARYO2VBS, cell 2723) (rtol 1e-10, atol 1e-30) |

## The generator's tests: make test-specialized

| Compiler | Exit | Time | PASS lines | FAIL lines | Last lines |
|---|---|---:|---:|---:|---|
| clang | 0 | 7 s | 1 | 0 | CC="clang" OPTFLAGS="-O3 -DNDEBUG" OMPFLAGS="-Xpreprocessor -fopenmp -isystem ~/homebrew/opt/libomp/include" OMPLIBS="-L~/homebrew/opt/libomp/lib -lomp" ~/NSF |
| gcc-15 | 0 | 10 s | 1 | 0 | CC="gcc-15" OPTFLAGS="-O3 -DNDEBUG" OMPFLAGS="-fopenmp" OMPLIBS="-fopenmp" ~/NSF-NCAR/tier2-work/verify/luna-high/gcc-15/luna-high/.conda/bin/python tests/test_specialized.py / PASS: 211 |

## Day-night scenarios (--threads 18, --rtol 1e-6 --atol 1e-15, 32-cell groups)

Agreement is score.py against MICM's result at the same settings, with `--level 10`; accuracy is score.py against the reference solution, for scale.

| Compiler | Scenario | Exit | elapsed_ms | Internal steps (MICM's) | Agreement with MICM | Worst species | Gate | Level 10 | Accuracy (MICM's) |
|---|---|---|---:|---|---:|---|---|---|---|
| clang | arctic-summer-10min | 0 | 3150.466 | 766 (103211) | 14.305 | HYDRALD 13.824 | PASS | PASS | 5.045 (5.04) |
| clang | night-start-5min | 0 | 4086.731 | 291 (132855) | 14.463 | ISOP 13.998 | PASS | PASS | 4.729 (4.73) |
| clang | south-48h-30min | 0 | 2777.192 | 1743 (95163) | 14.343 | ISOPO2VBS 13.766 | PASS | PASS | 4.966 (4.97) |
| clang | spring45-20min | 0 | 1785.747 | 1455 (60435) | 14.341 | ISOP 13.843 | PASS | PASS | 4.842 (4.84) |
| clang | tropics-polluted-15min | 0 | 2234.616 | 292 (75293) | 14.34 | MTERP 13.637 | PASS | PASS | 5.1 (5.10) |
| clang | clean-10min | 0 | 2550.007 | 288 (84273) | 14.34 | ISOP 13.801 | PASS | PASS | 4.907 (4.91) |
| clang | diurnal-10min | 0 | 2732.512 | 288 (90980) | 14.288 | ISOP 13.691 | PASS | PASS | 5.052 (5.05) |
| clang | diurnal-30min | 0 | 1422.977 | 320 (48619) | 14.374 | ISOP 13.817 | PASS | PASS | 4.748 (4.75) |
| clang | diurnal-5min | 0 | 4503.248 | 288 (145891) | 14.243 | ISOP 13.69 | PASS | PASS | 5.214 (5.21) |
| clang | polluted-10min | 0 | 3207.176 | 294 (106393) | 14.081 | ISOP 13.461 | PASS | PASS | 5.398 (5.40) |
| gcc-15 | arctic-summer-10min | 0 | 3882.367 | 766 (103211) | 14.305 | ISOP 13.831 | PASS | PASS | 5.045 (5.04) |
| gcc-15 | night-start-5min | 0 | 5117.135 | 291 (132855) | 14.456 | ISOP 13.996 | PASS | PASS | 4.729 (4.73) |
| gcc-15 | south-48h-30min | 0 | 3519.563 | 1743 (95163) | 14.342 | ISOPO2VBS 13.726 | PASS | PASS | 4.966 (4.97) |
| gcc-15 | spring45-20min | 0 | 2234.809 | 1455 (60435) | 14.347 | ISOP 13.781 | PASS | PASS | 4.842 (4.84) |
| gcc-15 | tropics-polluted-15min | 0 | 2793.338 | 292 (75293) | 14.337 | MTERP 13.626 | PASS | PASS | 5.1 (5.10) |
| gcc-15 | clean-10min | 0 | 3188.275 | 288 (84273) | 14.344 | ISOP 13.829 | PASS | PASS | 4.907 (4.91) |
| gcc-15 | diurnal-10min | 0 | 3429.270 | 288 (90980) | 14.295 | ISOP 13.689 | PASS | PASS | 5.052 (5.05) |
| gcc-15 | diurnal-30min | 0 | 1776.837 | 320 (48619) | 14.365 | ISOPNO3 13.61 | PASS | PASS | 4.748 (4.75) |
| gcc-15 | diurnal-5min | 0 | 5576.335 | 288 (145891) | 14.245 | NC4CH2OH 13.532 | PASS | PASS | 5.214 (5.21) |
| gcc-15 | polluted-10min | 0 | 4040.732 | 294 (106393) | 14.073 | ISOPNO3 13.36 | PASS | PASS | 5.398 (5.40) |

- clang: level 10 and the gate on every scenario: yes (10 scenarios)
- gcc-15: level 10 and the gate on every scenario: yes (10 scenarios)

## Build audit

`make -n -B cpu bench CC=clang`: exit 0; flags used: `-DCHEM_SPECIALIZED -DNDEBUG -Iinclude -L~/homebrew/opt/libomp/lib -O3 -Wall -Wextra -Wpedantic -Xpreprocessor -fopenmp -isystem -lm -lomp -p -std=c11`

Flags outside the rules (review): none

- `bin/test_core` links: /usr/lib/libSystem.B.dylib, ~/homebrew/opt/libomp/lib/libomp.dylib
- `bin/ts1_bench_cpu` links: /usr/lib/libSystem.B.dylib, ~/homebrew/opt/libomp/lib/libomp.dylib
- `bin/ts1_bench_specialized` links: /usr/lib/libSystem.B.dylib, ~/homebrew/opt/libomp/lib/libomp.dylib
- `bin/ts1_scenario_specialized` links: /usr/lib/libSystem.B.dylib, ~/homebrew/opt/libomp/lib/libomp.dylib

C sources scanned: 16 (tracked and generated, outside conda); lines matching a banned construct or starting a process: 75 (all in `logs/banned.txt`)

    include/chemistry.h:320: size_t vector_size);
    include/chemistry.h:370: size_t vector_size);
    include/chemistry.h:373: size_t vector_size);
    include/chemistry.h:376: size_t vector_size);
    include/chemistry.h:385: size_t vector_size);
    include/chemistry.h:388: size_t vector_size);
    include/chemistry.h:391: size_t vector_size);
    include/chemistry.h:394: size_t vector_size);
    include/chemistry.h:418: chem_vector_matrix_view *view, const int *vector, size_t vector_size,
    include/chemistry.h:449: double *vector, size_t vector_size);
    include/chemistry.h:452: size_t vector_size);
    include/chemistry.h:455: size_t vector_size);
    include/chemistry.h:458: size_t vector_size);
    include/chemistry.h:461: size_t vector_size);
    include/chemistry.h:465: size_t vector_size);
    include/chemistry.h:474: size_t vector_size, double scale);
    include/chemistry.h:558: size_t vector_size);
    include/chemistry.h:561: size_t vector_size);
    include/chemistry.h:564: size_t vector_size,
    src/matrix.c:106: size_t vector_size) {
    src/matrix.c:108: vector_size != view->columns || (vector_size && vector == NULL)) return 1;
    src/matrix.c:346: static int vector_matrix_vector_size(size_t rows, size_t columns,
    src/matrix.c:347: size_t vector_size) {
    src/matrix.c:348: return columns == 0 ? rows == 0 && vector_size == 0
    src/matrix.c:349: : rows <= (size_t)-1/columns && rows*columns == vector_size;
    src/matrix.c:354: size_t vector_size) {
    src/matrix.c:356: !vector_matrix_vector_size(view->rows,view->columns,vector_size) ||
    src/matrix.c:357: (vector_size && vector == NULL)) return 1;
    src/matrix.c:368: size_t vector_size) {
    src/matrix.c:370: !vector_matrix_vector_size(view->rows,view->columns,vector_size) ||

