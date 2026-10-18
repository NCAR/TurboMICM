# Blind review of four Tier 3 entries

You are one of two reviewers in the blind review of a controlled experiment that
compares LLM coding agents. Four agents received the same task, PROMPT.md in
this folder: starting from their own earlier port of MICM's CPU chemistry solver
to portable C, research, choose, design, and test numerical methods for the
fastest solution of atmospheric chemistry that is accurate enough for an
atmospheric model. Their final repositories are in `A/`, `B/`, `C/`, and `D/`,
in an order and with labels chosen at random. Each has had its git history,
local paths, and run names removed; `ROOT/entry` in a file stands for the
entry's own directory. `scenarios/` holds the scenarios' README and scoring
tool, as the agents saw them.

Speed and accuracy are ranked separately, by measurement. This review is the
secondary ranking, of the research and of the code, as PROMPT.md defines it:

- **Research**
  1. *Use of the literature:* the right sources found and read, used correctly
     and in proportion; every claim that rests on a source tied to it;
     REFERENCES.md complete and accurate about what was taken from each.
  2. *Rigor of the scans:* controlled series that change one choice at a time;
     data logged in `scans/` and summarized in SCANS.md; methods and structural
     choices covered widely enough to support the choices made; negative results
     kept; prototypes that measure accuracy and operation counts, each C method
     checked against its prototype.
  3. *Soundness of the accuracy argument:* ACCURACY.md justifies a target for an
     atmospheric model; the evidence supports each declared operating point and
     the recommended one; margins and limits are stated honestly, and the claims
     do not go beyond the evidence (for example, to undisclosed conditions).
- **Code**
  4. *Readability:* structure, naming, comments, and documentation that let a
     scientific software engineer follow the methods.
  5. *Maintainability:* tests, build, separation of concerns, error handling,
     and how well the documentation matches the code.
  6. *Extensibility:* how hard it is to add a method, as DESIGN.md says it
     should be done, and to use another mechanism.

## How to work

- Read only inside this folder. Do not look anywhere else on this machine or
  online for these entries, their authors, or other solutions; do not try to
  work out which agent wrote which entry.
- Do not build or run the code, or anything heavy: the machine may be timing
  programs. Read the files; `grep`, `wc`, and similar light commands are fine.
- Judge each entry on what is in its folder. Verify claims where you can from
  its own data and code (for example, that a scan's CSV supports SCANS.md, or
  that the C method matches its description), and say when you could not.
- Review the entries in the order you are given, then compare them.

## Report

Write your report as Markdown, in this order:

1. For each entry, a section with: a summary of its methods and operating
   points (three to five sentences); a score from 1 to 10 for each of the six
   criteria, each with two to four sentences of justification citing files;
   its three main strengths and three main weaknesses.
2. A comparison: the four entries ranked for research (criteria 1–3) and for
   code (criteria 4–6), with the reasons for the order, and any ties.
3. At the end, one JSON block, exactly in this form:

```json
{"scores": {"A": [1, 2, 3, 4, 5, 6], "B": [...], "C": [...], "D": [...]},
 "research_rank": ["?", "?", "?", "?"], "code_rank": ["?", "?", "?", "?"]}
```

with the six scores in the order of the criteria above, and the ranks from
best to worst.
