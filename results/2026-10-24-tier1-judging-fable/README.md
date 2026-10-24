# Tier 1 judging: Judge 2, Claude Fable 5.1, 2026-10-24

The second of three judgments of Tier 1's secondary criteria, readability,
maintainability and extensibility, by Claude Fable 5.1 at max effort, after
Claude Opus 5.5 (Judge 1, 2026-10-21) and before the SSE (Judge 3). The
criteria, the exercises and the scale are on the wiki's Tier-1-Judging page;
the report, with every score and its reasons, is Tier-1-Judging-Fable.

## Independence

- The judge scored before reading Judge 1's report page, its helpers' notes,
  exercise reports and test traces, or the "Extension points" sections of the
  entrant pages; it read them afterwards for the comparison section of the
  report, which was added in a later commit.
- The judge had seen Judge 1's totals and ranks (in the project's instructions
  and on Tier-1-Judging) and the README of Judge 1's results folder, which
  states its procedure, its scores and its findings for the SSE.
- The helpers were fresh sessions with the briefs in `brief/`, which are Judge
  1's briefs with "first judge" changed to "second judge", plus
  `brief/INDEPENDENCE.md`, which forbids reading Judge 1's materials.

## Procedure

1. **The judge's own reading** of each entry at its final commit, in the
   read-only worktrees of `~/NSF-NCAR/tier1-work/verify/<run>/clang/<run>/`:
   the documents, the Makefile, the public header, the benchmark driver, the
   integrator and parameter sets, the step control, the sparse LU and plan,
   the rates, the forcing and Jacobian, the reader and generator, the test
   runner, and the generated TS1 source. Notes with citations: `notes/`.
2. **Five helpers** (Claude, fresh general-purpose sessions in this Claude
   Code session), each to a brief in `brief/`: E1 and E2, E3, E4 and E5
   (`exercises/`), and the test sample in two halves, cases 1 to 6 and 7 to 12
   (`tests/T1.md`, `tests/T2.md`). The E5 helper copied each entry's Python
   tools into a scratch folder and ran them with `python3 -I` on MICM's
   tutorial configurations, on the same with the species in reverse order,
   and on a units probe, with its own checker (`exercises/e5-check/`: the
   checker, a stand-in `yaml.py` that reads JSON, the three added
   configurations, and the log of the 45 runs). Nothing in C was built or run.
3. **Reconciliation** by the judge: every helper rating and finding compared
   with the judge's notes; the code behind every finding that moves a score
   re-read; one standard chosen per exercise and applied to all nine. The
   report states each standard and the places where the judge departed from a
   helper's rating (E2 for astra-xhigh, E4 for luna-max).
4. **Scores and report**: `scores.json` and the wiki page.

## Scores

| Rank | Entrant | Readability | Maintainability | Extensibility (E1–E5) | Total of 30 |
|---:|---|---:|---:|---|---:|
| 1 | astra-high | 8 | 8 | 8.8 (5, 4, 5, 4, 4) | 24.8 |
| 2 | sol-xhigh | 7 | 8 | 8.4 (5, 4, 5, 3, 4) | 23.4 |
| 3 | astra-xhigh | 7 | 7 | 9.2 (5, 5, 5, 4, 4) | 23.2 |
| 4 | sol-high | 7 | 8 | 8.0 (5, 4, 4, 3, 4) | 23.0 |
| 5 | astra-max | 6 | 7 | 8.4 (5, 4, 5, 3, 4) | 21.4 |
| 6 | sol-max | 6 | 6 | 8.0 (5, 4, 4, 3, 4) | 20.0 |
| 7 | luna-max | 5 | 4 | 7.2 (4, 3, 5, 3, 3) | 16.2 |
| 8 | luna-xhigh | 4 | 5 | 6.4 (5, 2, 3, 3, 3) | 15.4 |
| 9 | luna-high | 3 | 2 | 4.4 (2, 2, 2, 3, 2) | 9.4 |

## Findings for the SSE

Bearing on the gate or the record, not on the scores; the report gives the
citations.

- luna-high's PASS lines are printed from a list of 510 ids once one program
  has exited 0; of the twelve sampled ids one has no corresponding check.
- luna-max reports all 52 LU ids from one dense check, the four Troe ids from
  one fixed case without MICM's parameters, and `DAESolveEnforcesConservation`
  from an unrelated constraint step.
- luna-xhigh takes one binary's exit status for 30 rate-constant ids and the
  reaction-rate store ids.
- All three Luna entries accept a step whose error norm is NaN, through the
  `fmax` floor; sol-high floors the same way but inspects the candidate.
- luna-high's generator emits wrong reactions for a CAMP mechanism whose
  species order its Markowitz pass changes; verified on the species-reversed
  user-defined tutorial (all ten reactions wrong). TS1 is unaffected.
- The three Luna readers convert units with the tutorial's fixed exponents; a
  bimolecular Arrhenius in molecule cm-3 units is left unconverted.
- luna-max's backward Euler has no step reductions.
- astra-xhigh installs generated kernels with no check against the run-time
  plan; sol-max checks only the nonzero count.
- astra-high recomputes the rate constants inside the timed loop; luna-xhigh
  builds its symbolic plan and allocates inside the timer.

## Files

- `brief/`: `EXERCISES.md`, `TEST_FIDELITY.md` (Judge 1's briefs, reworded for
  the second judge), `INDEPENDENCE.md`.
- `notes/<run>.md`: the judge's reading notes, nine files.
- `exercises/E1-E2.md`, `E3.md`, `E4-E5.md`: the helpers' exercise reports;
  `exercises/e5-check/`: `e5_check.py`, `yaml.py`, `run.log`, `configs/`.
- `tests/T1.md`, `T2.md`: the twelve traces.
- `draft/`: the sections of the report as written, before assembly.
- `scores.json`: the scores and ranks.
- `tools/spotcheck_page.py`: samples citations from the report page and prints
  the cited code.
