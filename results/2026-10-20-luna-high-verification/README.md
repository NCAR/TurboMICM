# Verification of luna-high, 2026-10-20

luna-high's rerun was rebuilt from its branch, retested, and audited with both
compilers. The SSE deferred its timing session on 2026-10-20 at 18:33; it was
timed with all twelve entries on 2026-11-17 (see [Timing](#timing)).

| File | Contents |
|---|---|
| `tests.md` | `make test-cpu` results for luna-high with each compiler |

## Builds, tests, and audit

luna-high (`d3c2480`, 128 commits, Codex CLI 0.163.0) was built with `make
env`, `make cpu bench OPTFLAGS="-O3 -DNDEBUG"`, and `make test-cpu`, once with
`CC=clang` and once with `CC=gcc-15`, by `admin/verify_tier3.py --stages
build,tests,ts1,audit`, the method of
[`2026-10-20-luna-xhigh-verification`](../2026-10-20-luna-xhigh-verification/README.md).

- **Tests** (`tests.md`): 662 PASS, 0 FAIL, 18 N/A with each compiler, every id
  exactly once, as the agent reported. Ten of the 18 N/A cases are N/A in no
  other entry: four `State` copy and move cases and six incompatible-ordering
  cases of MICM's vector sparse matrices. The PASS lines come from a fixed
  list of the cases the agent counts as ported, printed by `tools/test_cases.py`
  once `bin/test_core` and `tests/test_codegen.py` have both passed
  (`Makefile:48-51`; `tools/test_cases.py:640-778`): they stand for one
  passing C test program, and the mapping of its checks to MICM's cases is the
  agent's.
- **TS1 gate:** passes with both compilers; largest relative difference from
  MICM 5.34e-14 with clang (BCARYO2VBS, cell 9034) and 1.43e-13 with gcc-15
  (BCARYO2VBS, cell 4054).
- **Flags, libraries, sources:** only `-std=c11`, warnings, `-O3 -DNDEBUG`,
  include paths, and OpenMP flags; both binaries link only the OpenMP runtime
  and `libSystem`. The scan of the 15 C sources for banned constructs matched
  75 lines, all a parameter or variable named `vector_size`, not GCC's vector
  attribute.
- **Warnings:** none, from the library, the benchmark, or the tests, with
  either compiler.
- **Same work, same clock:** `src/main.c` makes one group per thread
  (`src/main.c:291-292`), allocates each group's solver workspace before the
  timer (`src/main.c:294-307`), runs one untimed warm-up interval
  (`src/main.c:308-320`), then times the STEPS intervals with C11's
  `timespec_get` (`src/main.c:187-191`, `329-349`), one `omp parallel for
  schedule(static)` over the groups calling the solver once per group and
  interval (`src/main.c:330-348`). It does MICM's default solver work at ten
  threads: 9,000 internal steps, 18,000 function calls, 9,000 Jacobian updates
  and decompositions, and 27,000 solves.

## Found while writing Tier 1 Implementations

- **A generator defect.** For a configuration in MICM's CAMP format,
  `tools/read_config.py` gives each term the index of its species' input
  position (`tools/read_config.py:83-86`), and `tools/generate_c.py` emits that
  index (`tools/generate_c.py:126`) beside a species table in its Markowitz
  order (`tools/generate_c.py:100-102,130-137`). When the order differs from
  the input's, the generated reactions name the wrong species. Reproduced
  with copies of the two scripts on a four-species configuration (A + B → C,
  A + C → D, A + D → B, D → A): the species table came out B, A, C, D, and
  D → A was emitted as D → B. TS1's JSON carries no indices
  (`tools/import_ts1.py:56-61`), and the tests check the generated text only
  for the two tutorial configurations, whose order is unchanged
  (`tests/test_codegen.py:108-114`).

## Timing

Timed on 2026-11-17, 20:49 to 23:39, with every other entry and MICM's three
programs, in
[`2026-11-17-tier1-all-twelve-timing`](../2026-11-17-tier1-all-twelve-timing/README.md):
5 rounds at 10 and 1 threads, every run passing the gate.

| | 10 threads, median ms (range) | vs MICM's default | vs MICM's fastest | 1 thread, median ms (range) | vs MICM's default | vs MICM's fastest |
|---|---:|---:|---:|---:|---:|---:|
| luna-high | 50,942 (50,542–51,727) | 0.54× | 0.19× | 311,252 (308,706–312,341) | 0.47× | 0.16× |

It is the last of the twelve entries at both thread counts, and gains 6.11× from
one thread to ten, the most of any program timed. Its self-report, 48.5 s
against 27.3 s for MICM's default in one ten-thread run, compares with 50,942
and 27,265 ms here.
