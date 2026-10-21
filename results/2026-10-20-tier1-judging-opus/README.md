# Tier 1 judging: Judge 1, Claude Opus 5.5, 2026-10-20 and 21

The first of three judgments of Tier 1's secondary criteria, readability,
maintainability, and extensibility, by Claude Opus 5.5 at max effort; Claude
Fable 5.1 at max effort and the SSE follow. The criteria, the exercises, and
the scale are on the wiki's Tier-1-Judging page; the report, with every score
and its reasons, is Tier-1-Judging-Opus; each entrant's design is on its
Entrant-<run> page.

## Decisions

- **Three judges in turn** (the SSE, 2026-10-20), in place of the blind review
  by fresh sessions that Tier-1-Design had planned.
- **Judge 1 in the open** (the SSE, 2026-10-20): the session that judged is the
  one that orchestrates the experiment and had read the run records, so it knew
  each entry's model, effort, and speed.
- **Exercises on paper** (the SSE, 2026-10-20): the changes are worked out from
  the code as it stands, and nothing in the entries is built or run.
- **The rubric** (proposed by Judge 1 and stated to the SSE before the evidence
  was gathered): three scores of 1 to 10; extensibility is twice the mean of
  five exercise ratings of 1 to 5; equal weights; ties go to extensibility,
  then maintainability. The extensibility score is kept to one decimal.

## Procedure

1. **One counting rule for all nine** (`tools/metrics2.py`, `tools/funcs.py`):
   line counts by category, C function sizes by brace depth, and the share of
   lines with a comment, at each final commit from tracked files
   (`metrics.json`, `metrics_table.md`).
2. **Read-only helper agents** (Claude, in this session), each to a brief in
   `brief/`:
   - nine page agents wrote the entrant pages to `brief/PAGE_TEMPLATE.md`, with
     notes for the judge (`notes/`);
   - three exercise agents worked E1 and E2, E3, and E4 and E5 for all nine
     entries to `brief/EXERCISES.md` (`exercises/`); the E5 agent ran copies of
     each entry's reader and generator on the configurations in
     `exercises/e5-inputs/`, in its own scratch folder, with the checkers in
     `exercises/e5-check/`; the copies and their outputs are not kept here,
     since they are the entries' code;
   - two test agents traced twelve of MICM's cases in all nine entries to
     `brief/TEST_FIDELITY.md` (`tests/`).

   The page agents read condensed timelines of each session, made from the
   transcripts by `tools/transcript_digest.py`; the digests stay off `main`,
   like the transcripts, and the script remakes them from the `transcripts`
   branch.
3. **Checks by Judge 1.** 54 citations, six at random from each entrant page
   (`tools/spotcheck.py`), all held; every diagram rendered with mmdc
   (`tools/checkmmd.py`); and each finding that moves a score by more than a
   point confirmed in the code. One rating changed: luna-xhigh's E4 from 4 to 3
   (its reader has no order-dependent unit rule for a new type in MICM's
   units). One fact about C was checked by compiling a three-line program in a
   scratch folder: `fmax(NaN, 1e-10)` returns 1e-10.
4. **Scores and report** by Judge 1: `scores.json` and the wiki page.

`<scratch>` in the briefs and reports stands for the session's scratch folder.

## Scores

| Rank | Entrant | Readability | Maintainability | Extensibility (E1–E5) | Total of 30 |
|---:|---|---:|---:|---|---:|
| 1 | astra-xhigh | 7 | 8 | 8.8 (5, 5, 5, 4, 3) | 23.8 |
| 2 | sol-xhigh | 8 | 7 | 7.6 (5, 4, 4, 3, 3) | 22.6 |
| 3 | astra-high | 6 | 8 | 8.4 (5, 4, 5, 3, 4) | 22.4 |
| 4 | astra-max | 7 | 7 | 8.4 (5, 4, 5, 3, 4) | 22.4 |
| 5 | sol-high | 6 | 7 | 7.6 (5, 4, 3, 3, 4) | 20.6 |
| 6 | sol-max | 6 | 6 | 7.2 (5, 4, 3, 3, 3) | 19.2 |
| 7 | luna-xhigh | 5 | 4 | 6.0 (5, 3, 2, 3, 2) | 15.0 |
| 8 | luna-max | 4 | 3 | 6.8 (4, 4, 4, 3, 2) | 13.8 |
| 9 | luna-high | 3 | 2 | 4.4 (2, 2, 2, 3, 2) | 9.4 |

astra-high and astra-max tie on total and extensibility; astra-high ranks
third on maintainability.

## Findings for the SSE

Bearing on the primary gate or the record, not on the scores; the report gives
the citations.

- luna-high's 662 PASS lines are printed from a list once one test program
  passes; of the twelve sampled cases, one has no corresponding check.
- luna-max reports 363 of the 680 ids from seven shared checks; in luna-max and
  luna-xhigh the LU and linear-solver ids test dense routines, not the sparse
  LU the solver uses, and the Markowitz cases run no Markowitz ordering.
- All three Luna entries floor a NaN error norm to 1e-10 with `fmax`;
  luna-xhigh's three-stage `e[1]` differs from MICM's by 2.7e-14 relative;
  luna-max's backward Euler has no step-size control.
- The three Luna readers convert MICM's molecule cm-3 units with fixed
  exponents (a bimolecular Arrhenius rate comes out 6.0e17 times too small),
  and luna-high's generator names the wrong species when its Markowitz order
  permutes the input; TS1 is unaffected.
- astra-xhigh's context was compacted four times, not twice (corrected in
  `admin/runs.md`); astra-high was first timed in
  `results/2026-10-15-micm-configurations`.

## Files

- `brief/`: the three briefs.
- `notes/`: the page agents' notes on readability and maintainability, one per
  entrant.
- `exercises/`: the reports for E1–E2, E3, and E4–E5, with E5's inputs and
  checkers.
- `tests/`: the two test-sample reports, T1 (cases 1–6, with the runners) and T2
  (cases 7–12).
- `metrics.json`, `metrics_table.md`: the common counts.
- `scores.json`: Judge 1's scores and rank.
- `tools/`: `metrics2.py` and `funcs.py` (counts), `transcript_digest.py`
  (session timelines), `checkmmd.py` (renders a page's Mermaid diagrams),
  `spotcheck.py` (samples a page's citations).
