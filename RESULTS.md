# Validation and measurements

Validated on the provided 10-core M5 machine, macOS 27.0.1, with Apple clang 21.0.0 and GCC 15.2.0. Both C builds and the independent clang++ reference use `OPTFLAGS='-O3 -DNDEBUG'`. The C build adds C11, warnings, include paths, shared-library support, and the prescribed OpenMP flags. There are no target, floating-point, LTO, or profile-guided optimization flags.

## Required cases

| Build | PASS | N/A | FAIL | Reported case lines |
| --- | ---: | ---: | ---: | ---: |
| `make test-cpu CC=clang` | 675 | 5 | 0 | 680 |
| `make test-cpu CC=gcc-15` | 675 | 5 | 0 | 680 |
| `make test-cpu CC=clang KERNEL_CHUNK=7` | 675 | 5 | 0 | 680 |

The identifiers and ordering match `../reference/micm_cases.txt` exactly. Of the passing cases, 667 are GoogleTest equivalents and eight are native C tutorial programs. The five N/A cases are the four compile-only `ViewCategorySFINAE` cases and `GroupedView.CategoryTagsAreDistinct`: C has no corresponding concepts, template category tags, or overload resolution. Runtime grouped-view and vector behavior remains tested.

The suite checks all reference rate families, dense/grouped/sparse arrays and views, copying and error returns, both LU algorithms, finite-difference Jacobians, all five Rosenbrock tables, backward Euler, analytical and stiff reference problems, algebraic initialization and rollback, external model composition, conservation, overshoot, and shared-model OpenMP execution. The eight tutorial outputs match the recorded reference stdout exactly after excluding elapsed time and interleaved thread-startup diagnostics, as detailed in `DECISIONS.md`.

All three supplied YAML/JSON examples are compared with independently code-built mechanisms. Tests compare rate constants, forcing, every Jacobian entry, and integrated species values with a bound of `1e-30 + 1e-10 * abs(reference)`. Both generic and generated paths run. An additional configuration exercises every serializable rate family, a third body, repeated reactants, a 16-coefficient Taylor polynomial, and distinct cell conditions. Generated TS1 kernels are compared directly with generic evaluation and integration through the same library.

## TS1 accuracy gate

Every run below includes one warm-up and 100 timed steps of 30 seconds. The supplied comparison tool checks every species and condition in every cell at `atol=1e-30`, `rtol=1e-10`.

| Input | Compiler | Cells × output columns | Maximum relative difference | Gate |
| --- | --- | ---: | ---: | --- |
| supplied seed 1 | clang | 100 × 211 | 3.48e-14 | PASS |
| supplied seed 1 | gcc-15 | 100 × 211 | 3.43e-14 | PASS |
| generated seed 314159 | clang | 1,000 × 211 | 6.46e-14 | PASS |
| generated seed 314159 | gcc-15 | 1,000 × 211 | 5.67e-14 | PASS |
| generated seed 2026 | clang | 10,000 × 211 | 1.26e-13 | PASS |
| generated seed 2026 | gcc-15 | 10,000 × 211 | 1.16e-13 | PASS |

There are 209 variable species plus temperature and pressure in each output row; the model has 210 species including its parameterized third body, and 547 reactions. The 10,000-cell comparison therefore checks 2,110,000 values per run.

Additional passing gates cover groups of 1, 17, 1000, partial final groups, 1/2/3/4/default workers, static/dynamic/guided scheduling, and generated kernel chunks of 7 and 32. Large clang runs with groups of 16, 32, and 64 all pass. The one-worker and all-worker MICM 10,000-cell results agree exactly for this input.

The 10,000-cell input SHA-256 is `43febf31a0a1f8e0acfe02d93ac96d95ecc024df9572547dd1000926a85ad8b5`. Its independently generated one-worker MICM output SHA-256 is `133290093c64317d520dff78e89622db36dc467110eb0f3d0e4bc6b24f16ba7d`.

## Performance

These are local elapsed times reported by the native programs, excluding initialization, warm-up, and CSV output. The main comparison uses 10,000 distinct cells, seed 2026, 100 timed steps of 30 seconds. Group size is 32 for the C port. Reference measurements are individual runs; the main clang result is the median of five runs. They are measurements of this machine, rather than a performance guarantee.

| Program | Workers | Elapsed ms | Speedup over matching MICM run |
| --- | --- | ---: | ---: |
| MICM, clang++ | 1 | 136,187.244 | 1.00× |
| C port, clang | 1 | 36,954.738 | 3.69× |
| C port, gcc-15 | 1 | 40,535.136 | 3.36× |
| MICM, clang++ | all 10 | 26,010.763 | 1.00× |
| C port, clang, median of 5 | all 10 | 6,935.493 | 3.75× |
| C port, gcc-15 | all 10 | 7,329.128 | 3.55× |

The five clang times were 6,813.886, 7,171.708, 6,855.476, 6,935.493, and 7,091.086 ms; every run passed the complete gate. The 100-cell, one-worker runs took 382.289 ms with clang and 407.043 ms with GCC. A prior same-flags MICM 100-cell measurement took 1,197.123 ms.

A coarse all-worker group-size sweep gave 6,688.210/7,088.742 ms for group 16, 6,813.886/7,171.708 ms for group 32, and 7,408.583/7,663.285 ms for group 64. Groups 16 and 32 were close relative to run variation. The default remains 32, which halves the number of adaptive controllers relative to group 16 while retaining a modest working set. The setting remains a build/run parameter.

Timed solver totals for the 10,000-cell C run were 281,700 internal/accepted steps, zero rejections, 563,400 forcing calls, 281,700 Jacobian updates/decompositions, and 845,100 solves. There are 313 independent groups. MICM's one-worker run has one group and reports 900 accepted steps, 1,800 forcing calls, 900 Jacobian updates/decompositions, and 2,700 solves; its ten-worker run reports ten times those counts. Both implementations take nine internal steps per group per timed call. Grouping accounts for the different totals; coefficients, tolerances, and step-size rules are unchanged.

## Reproduction

```sh
make env
mkdir -p results
make test-cpu CC=clang > results/clang-tests.log
make test-cpu CC=gcc-15 > results/gcc-tests.log
make micm-ref
.conda/bin/python ../reference/make_cells_csv.py 10000 2026 --output results/cells.csv
mkdir -p results/reference-run
(cd results/reference-run && ../../bin/micm_ts1_ref --threads 1 10000 100 30 ../cells.csv)
make bench CC=clang
.conda/bin/python tools/check_benchmark.py 10000 100 30 \
  results/cells.csv results/reference-run/ts1_final.csv results/clang-run \
  --compiler clang
```

Omit `--threads 1` on the reference command for its all-worker timing. Add `--threads 1`, `--group N`, or `--schedule dynamic,1` to the check tool for other settings. Select `CC=gcc-15` and label that compiler in the check tool to repeat with GCC. The tool retains stdout, the complete final CSV, gate output, input/reference hashes, and command metadata under the chosen output directory. Generated CSVs, build products, environment, caches, and local run logs are ignored by git.

The portability audit found only standard C headers and `omp.h` in runtime sources, and only CPU OpenMP pragmas. The native executable links the system C runtime and OpenMP; it has no Python, MICM, or MUSICA runtime dependency.
