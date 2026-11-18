# Tier 1 judging: Judge 1 on the three Terra entries, 2026-11-18

Judge 1, Claude Opus 5.5 at max effort, scored the three Terra entries (terra-max, terra-xhigh,
terra-high) on Tier 1's secondary criteria, readability, maintainability, and extensibility, to
the standard of its judgment of the first nine entries
([`2026-10-20-tier1-judging-opus`](../2026-10-20-tier1-judging-opus/README.md)), and ranked them
among the twelve. Judge 2 (Claude Fable 5.1) and Judge 3 (the SSE) have not scored the Terra
entries. The criteria, the exercises, and the scale are on the wiki's Tier-1-Judging page; the
report is Tier-1-Judging-Opus, which now covers all twelve.

## Decisions

- **The same judge and rubric** as for the nine: three scores of 1 to 10; extensibility is twice
  the mean of the five exercise ratings, kept to one decimal; equal weights; ties go to
  extensibility, then maintainability.
- **In the open**, as before: the session that judged is the one that orchestrates the
  experiment, and it had verified and timed the Terra entries and written their entrant pages and
  their sections of Tier-1-Implementations.
- **No helper agents.** For the nine, read-only helpers gathered the evidence to Judge 1's briefs.
  For the Terra entries Judge 1 gathered it itself, to the same briefs (`brief/` in the earlier
  folder) and the standards the helpers' reports state, so that the three are rated as the nine
  were. Where a case or an exercise needed a standard the earlier reports did not state, the
  standard chosen is stated in the report (below).
- **The nine entries' scores stand** as committed on 2026-10-21.
- **Exercises on paper**, nothing built. For E5, copies of the three readers and generators were
  run in a scratch folder, as for the nine.

## Procedure

1. **Counts** by Judge 1's rule (`tools/metrics_terra.py`, which runs
   `../2026-10-20-tier1-judging-opus/tools/funcs.py`): `metrics.json`, `metrics_table.md`.
2. **Reading**, at each final commit, in the read-only worktrees
   `~/NSF-NCAR/tier1-work/verify/<run>/clang/<run>/`: the documents, the public headers, the
   benchmark driver, the integrators and parameter sets, the step control, the LU and its plan, the
   rates, forcing and Jacobian, the reader and generator, and the case runner (`notes/`).
3. **The five exercises** for the three entries (`exercises/E1-E5-terra.md`). For E5, 24 runs of
   the copied generators, with Xcode's Python 3.9.6 in isolated mode (`tools/run_copy.py`), on
   Judge 1's inputs (`../2026-10-20-tier1-judging-opus/exercises/e5-inputs/`: the two MICM
   configurations, `permute`, and the third-body probe, here with M marked by
   `parameterization: air density` for terra-max and `third body: true` for the other two),
   checked by `tools/terra_terms_check.py` (`exercises/e5-check/run.log`, `terms.log`). The
   copies and their outputs are not kept here, since they are the entries' code.
4. **The test sample**: the twelve cases of Tier-1-Judging traced in the three entries
   (`tests/T-terra.md`), with MICM's assertions as summarized in the earlier folder's
   `tests/T1.md` and `T2.md`.
5. **Scores and report**: `scores.json` and the wiki pages.

Standards set for the Terra entries, in the report: in the test sample, a check that reaches a
different implementation from the one the solver path uses (a searching Mozart where the solver
replays a plan, a C ordering where the tables come from a Python one, a dense helper where the
solver has grouped kernels) is P, after Judge 1's grade for astra-high's case 6; for E2, an
entry with no hook whose controller is written once is rated 4, the general scale's "a few edits
in one module", where Judge 1's 3 was for no hook and the logic in more than one integrator.

## Scores

| Rank of 12 | Entrant | Readability | Maintainability | Extensibility (E1–E5) | Total of 30 |
|---:|---|---:|---:|---|---:|
| 2 | terra-max | 7 | 7 | 8.8 (5, 5, 5, 3, 4) | 22.8 |
| 8 | terra-high | 7 | 5 | 7.2 (5, 4, 2, 3, 4) | 19.2 |
| 9 | terra-xhigh | 6 | 6 | 6.8 (5, 3, 2, 3, 4) | 18.8 |

Among the twelve: astra-xhigh 23.8, terra-max 22.8, sol-xhigh 22.6, astra-high and astra-max
22.4, sol-high 20.6, sol-max and terra-high 19.2, terra-xhigh 18.8, luna-xhigh 15.0, luna-max
13.8, luna-high 9.4. sol-max and terra-high tie on total and extensibility; sol-max ranks seventh
on its higher maintainability.

Test sample, faithful of 12: terra-max 7, terra-xhigh 9, terra-high 8.

## Findings for the SSE

Bearing on the gate or the record, not on the scores; the report gives the citations.

- terra-high's 675 PASS lines are the exit statuses of the programs a table names for the ids;
  every id of a program shares its verdict, and all 32 LU ids rest on one program whose matrices
  are diagonally dominant (with the 36 vector-ordered ids of its verification record).
- terra-max's LU and linear-solver ids run the library's searching Mozart, and its forcing and
  Jacobian ids a per-cell dense helper, not the plan and grouped kernels the solver runs; no test
  calls those directly. terra-xhigh's and terra-high's Markowitz ids test a C ordering, where TS1's
  tables come from a Python copy of the rule.
- For other mechanisms, the three Terra readers convert MICM's molecule cm-3 units by reaction
  order, as the six Sol and Astra entries do, leave user-defined rates unscaled, and kept every
  species through their orderings in the E5 runs: no wrong-result defect was found.

## Files

- `notes/<run>.md`: Judge 1's reading notes, three files.
- `exercises/E1-E5-terra.md`: the five exercises; `exercises/e5-check/`: the E5 run log and the
  checker's log.
- `tests/T-terra.md`: the twelve traces in each entry, with the runners.
- `metrics.json`, `metrics_table.md`: the counts.
- `scores.json`: the three entries' scores and the rank of the twelve.
- `tools/`: `metrics_terra.py` (counts), `run_copy.py` (runs a copied generator in isolated
  mode), `terra_terms_check.py` (maps generated species indices to names and compares the
  reactions with the configuration's).

`<scratch>` in the logs stands for the session's scratch folder.
