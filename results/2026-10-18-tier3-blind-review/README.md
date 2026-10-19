# Tier 3 blind review, 2026-10-18

The secondary ranking of the four Tier 3 entries: a blind review of the
research (use of the literature, rigor of the scans, soundness of the
accuracy argument) and of the code (readability, maintainability,
extensibility), by two reviewers who did not know which agent wrote which
entry. The procedure follows `tier3/HANDOFF.md` and the reviewers'
instructions, `tier3/review/REVIEW.md`.

## Procedure

- **Anonymized copies.** `admin/anonymize_tier3.py` exported each entry's
  final commit (`git archive`, no history) without PROMPT.md, replaced the
  local paths, run names, and model and harness names in every text file,
  and searched every file for leftovers; its one hit was bytes inside a PNG
  (`anonymize.log`). The letters A to D were drawn at random
  (`secrets.SystemRandom().shuffle`) and the key kept in a private file
  until both reports were in. The review folder held the four copies
  (141, 201, 158, and 159 files), REVIEW.md, the Tier 3 prompt with its
  paths and run names abbreviated, the scenarios' README and `score.py`,
  and the reference README.
- **Fresh sessions.** Each reviewer was a new Claude Code session started
  from the review folder (so no project instructions or memory loaded) at
  effort max, with only the Read, Grep, Glob, and Bash tools, Bash limited
  to `wc`, `grep`, `head`, `tail`, and `cat`, and the web, agents, and
  editing tools disabled (`run-opus.sh`, `run-fable.sh`). The prompt told
  each to follow REVIEW.md, to read only inside the folder, not to run
  code, and to end with the JSON block. Opus read the entries A to D and
  Fable D to A. Both started at 19:43:37; Opus finished at 20:03:21
  (159 turns) and Fable at 20:01:06 (180 turns). The machine was running
  the entries' re-verification builds at the time and nothing was being
  timed.
- **Reports.** `report-opus.md` and `report-fable.md` are the sessions'
  final messages, unedited; `session-opus.json` and `session-fable.json`
  are the sessions' result records with their usage.

## The key

| Letter | Entry |
|---|---|
| A | sol-xhigh |
| B | astra-xhigh |
| C | astra-max |
| D | sol-max |

## Scores

Criteria: 1 literature, 2 scans, 3 accuracy argument (research); 4
readability, 5 maintainability, 6 extensibility (code); each 1 to 10.

| Entry | Reviewer | 1 | 2 | 3 | Research | 4 | 5 | 6 | Code |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| sol-xhigh (A) | Opus | 8 | 8 | 8 | 24 | 6 | 8 | 7 | 21 |
| | Fable | 9 | 8 | 8 | 25 | 8 | 8 | 7 | 23 |
| astra-xhigh (B) | Opus | 7 | 9 | 8 | 24 | 6 | 8 | 6 | 20 |
| | Fable | 8 | 9 | 9 | 26 | 7 | 8 | 7 | 22 |
| astra-max (C) | Opus | 7 | 8 | 7 | 22 | 7 | 8 | 7 | 22 |
| | Fable | 8 | 9 | 9 | 26 | 8 | 9 | 8 | 25 |
| sol-max (D) | Opus | 7 | 7 | 8 | 22 | 7 | 7 | 7 | 21 |
| | Fable | 7 | 8 | 9 | 24 | 6 | 8 | 7 | 21 |

## Ranks

| Ranking | Opus | Fable |
|---|---|---|
| Research | astra-xhigh, sol-xhigh, astra-max, sol-max | astra-max, astra-xhigh, sol-xhigh, sol-max |
| Code | astra-max, sol-xhigh, sol-max, astra-xhigh | astra-max, sol-xhigh, astra-xhigh, sol-max |

- **Where they agree:** astra-max has the best code (first for both), and
  sol-max the weakest research (last for both); sol-xhigh is second for
  code with both; astra-xhigh is first or second for research with both.
- **Where they differ:** the research winner. Opus put astra-xhigh first
  (with sol-xhigh tied on points) for the breadth of its scans and the
  depth of its accuracy analysis, and astra-max third; Fable put astra-max
  first (tied on points with astra-xhigh) for the rigor of its experimental
  record (named parents, hashed rows, a defect found and regression-tested)
  and sol-xhigh third. On code, they swap sol-max and astra-xhigh in the
  last two places.
- **Spreads are small:** 22 to 26 points of 30 for research and 20 to 25
  for code. Both reviewers note that no entry used its Python prototypes
  to explore methods on TS1; the prototypes verify the C on small problems.

Each report's per-entry sections give the justification for every score
and the three main strengths and weaknesses of each entry.
