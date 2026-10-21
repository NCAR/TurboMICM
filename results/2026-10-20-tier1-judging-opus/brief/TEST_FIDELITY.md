# Test-port fidelity across the nine Tier 1 entries

You are helping the first judge of TurboMICM's Tier 1 assess
**maintainability**, one of three secondary criteria. Part of it is how
faithfully each entry ported MICM's tests: a test suite protects a later
change only if each case still checks what MICM's case checks. Every entry
reports MICM's 680 case ids as PASS, FAIL, or N/A from `make test-cpu`, and
all nine pass. This task looks behind those lines, on a fixed sample of cases
that is the same for all nine entries.

## The contest, in short

Each entry is one OpenAI model at one reasoning effort that ported MICM's CPU
chemistry solver (C++) to ISO C11, with libm and OpenMP, and ported MICM's 680
test cases (672 GoogleTest cases in 59 executables and MICM's 8 tutorial
programs), in C or in Python through ctypes. `make test-cpu` prints one line
per case: `<id> PASS|FAIL|N/A <reason>`, with the ids of
`~/NSF-NCAR/ODE/benchmarks/reference/micm_cases.txt`. N/A is
allowed only where a case tests C++ mechanics with no C equivalent. Rules:
`~/NSF-NCAR/TurboMICM.wiki/Tier-1-Design.md`.

## Where things are (read-only)

- MICM's tests: `~/NSF-NCAR/ODE/benchmarks/micm/test/` (unit,
  integration, regression, tutorial). Find a case with, for example,
  `grep -rn "Robertson" micm/test`; a test may be a `TEST(Suite, Name)` that
  calls a templated helper in a `*_policy.hpp` file: follow it to the
  assertions. A program case's expected output is in MICM's tutorial sources
  or the entry's stored reference output.
- The entries: git worktrees at each final commit,
  `~/NSF-NCAR/tier1-work/verify/<run>/clang/<run>/`, for `<run>`
  in sol-high, sol-xhigh, sol-max, astra-high, astra-xhigh, astra-max,
  luna-max, luna-xhigh, luna-high. Start from each entry's `make test-cpu`
  target in its Makefile and its case runner (for example `tests/run.py`,
  `tests/run_cases.py`, `tests/runner.c`, `tools/run_tests.py`,
  `tools/test_cases.py`, `tests/main.c`) to find how an id reaches a check.
- `~/NSF-NCAR/TurboMICM.wiki/Tier-1-Implementations.md` and the
  Notes on `~/NSF-NCAR/TurboMICM.wiki/Tier-1-Results.md` describe
  each entry's tests briefly; check against the code.

Read-only: do not edit, build, or run anything. Light commands (`grep`,
`sed -n`, `wc`, reading files) are fine.

## For each case in your sample

1. Read MICM's test: what scenario it sets up and exactly what it asserts
   (values and tolerances, exceptions, counts). Summarize in two to four
   lines, citing MICM's file and lines.
2. For each of the nine entries, find the code that decides this id's
   verdict, and grade it:
   - **F, faithful:** the same scenario and the same assertions (values,
     tolerances, the error raised), through the entry's C API;
   - **P, partial:** the same scenario with fewer or weaker assertions, or
     with different inputs or looser tolerances (say which);
   - **S, surrogate:** a different or shared check stands in (for example one
     smoke test reported for many ids);
   - **M, mapped only:** the id's verdict comes from a list or mapping, with
     no check specific to it, or with a check you cannot find;
   - **N/A:** the entry marks it N/A (quote its reason).
   Cite `path:line` for the check and for the place that maps the id to it.

## Also, for each entry

One paragraph on the runner: how `make test-cpu` turns checks into the 680
lines. Does each id run its own check (a C function or Python test per case,
selected by name), does one program report many ids, or does a list print
PASS lines once a program has passed? Cite the code.

## What to hand back

A report at the path your task names: a summary table (rows: the sample's
cases; columns: the nine entries; cells: F, P, S, M, or N/A), a count of each
grade per entry, the runner paragraphs, then per case MICM's assertions and,
per entry, the grade with its citation and a line of explanation. Facts only,
no words about the agents or models beyond the run names. Your final message:
the summary table, the per-entry counts, and anything that looks like a case
reported PASS without a check behind it.
