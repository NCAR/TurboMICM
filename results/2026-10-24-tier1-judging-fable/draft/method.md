# Tier 1 Judging: Judge 2, Claude Fable 5.1

*The second of three judgments of Tier 1's secondary criteria, readability,
maintainability and extensibility, by Claude Fable 5.1 at max effort, on
2026-10-24. Judge 1 was Claude Opus 5.5 ([Tier-1-Judging-Opus]); Judge 3 is
the SSE. The criteria, the exercises and the scale are on
[[Tier 1 Judging|Tier-1-Judging]]; each entrant's design is on its
`Entrant-<run>` page. Nothing in this report bears on the primary ranking,
TS1 wall time, or on the correctness gate, which all nine pass.*

## How this judgment was made

**In the open, like Judge 1.** The session that judged is the one that
orchestrates the experiment. It knew each entry's model, effort and speed.

**What the judge had and had not read.** Before scoring, the judge had read
the project's instructions and the [[Tier 1 Judging|Tier-1-Judging]] page,
both of which show Judge 1's totals and ranks, and the README of Judge 1's
results folder, which lists its procedure, its scores and its findings for
the SSE. The judge had not read Judge 1's report page, its helpers' notes,
exercise reports or test traces, or the "Extension points" sections of the
entrant pages, which summarize Judge 1's exercises; those were read only after
the scores below were fixed and committed, for the comparison at the end of
this page. The scores were formed from the code.

**Procedure.**

1. **The judge read every entry.** For each of the nine, at its final commit
   in a read-only worktree: README.md, DESIGN.md, DECISIONS.md, the Makefile,
   the public header, the benchmark driver, the integrator and its parameter
   sets, the step control, the sparse LU and its plan, the rate evaluation,
   the forcing and Jacobian, the configuration reader and generator, and the
   test runner, plus the generated TS1 source where the entry has one. Notes
   with citations are in `notes/<run>.md` of the results folder.
2. **Five helpers, fresh sessions, the same briefs as Judge 1's.** Three
   worked the exercises (E1 and E2; E3; E4 and E5) and two traced the twelve
   sampled test cases (cases 1 to 6; 7 to 12) in all nine entries, each to the
   briefs in `brief/`, which are Judge 1's with "first judge" changed to
   "second judge" and an added instruction not to read Judge 1's materials.
   The helpers did not know the judge's own ratings. For E5 only, the helper
   copied each entry's Python tools into a scratch folder and ran them on
   MICM's tutorial configurations, on the same configurations with the species
   listed in reverse order, and on a units probe, with a checker of its own
   (`exercises/e5-check/`); nothing in C was built or run, by the helpers or by
   the judge.
3. **Reconciliation.** The judge compared the helpers' ratings and findings
   with its own notes and re-read the code behind every finding that moves a
   score; where a helper's standard differed from the judge's first reading
   the judge chose one standard and applied it to all nine. The differences
   are noted under each exercise.
4. **Scores and report** by the judge: `scores.json` and this page. Compiler
   warnings come from the verification records, which note them only for the
   three Luna entries; nothing was rebuilt for this judgment.

**The rubric** is Judge 1's, stated on [[Tier 1 Judging|Tier-1-Judging]]:
readability and maintainability 1 to 10 each; extensibility twice the mean of
five exercise ratings of 1 to 5, kept to one decimal; equal weights, a total
of 30; ties go to extensibility, then maintainability.
