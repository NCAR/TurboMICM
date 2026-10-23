# Tier 2 blind review, 2026-10-23

Tier 2's secondary ranking, as Tier-2-Design defines it: a blind review of
each entry's generator for readability, maintainability, and extensibility.
The reviewers did not know which agent wrote which entry. Their instructions
are [`tier2/review/REVIEW.md`](../../tier2/review/REVIEW.md), and the three
exercises they did on paper were also measured on the machine
([`2026-10-23-tier2-exercises`](../2026-10-23-tier2-exercises/README.md)).

## Procedure

- **Anonymized copies.** `admin/anonymize_tier2.py` exported each entry's final
  commit with `git archive` (no history), without PROMPT.md.
  - In every text file it replaced the local paths, run names, and model and
    harness names; five files changed in all.
  - Then it searched every file for anything left, and found nothing
    (`anonymize.log`).
  - At 00:26 Claude Code drew the letters A to I at random (a Fisher-Yates
    shuffle with Python's `secrets`) and kept the key outside the review
    folder until both reports were in.
- **The review folder** held:
  - the nine copies, of 30 to 128 files;
  - REVIEW.md;
  - the Tier 2 prompt, with its paths and run names abbreviated;
  - the scenarios' README and `score.py`, and the reference README;
  - MUSICA's v0 Chapman and Carbon Bond 5 configurations, for exercise (c).
- **Fresh sessions.** Each reviewer was a new Claude Code session (`claude -p`)
  started from the review folder, so no project instructions or memory
  loaded, at effort max (`tier2/review/run-review.sh`).
  - Its tools were Read, Grep, Glob, and Bash limited to `wc`, `grep`, `head`,
    `tail`, and `cat`; the web, agents, and editing were disabled.
  - The prompt told each to follow REVIEW.md, to read only inside the folder,
    to run no code, to do the exercises on paper, and to end with the JSON
    block.
  - Their transcripts show no path outside the folder. The only refused
    commands were a few `ls` and `cd`.
- **The reviewers: a deviation.** The design names Opus and Fable at max.
  - Fable's session, started at 01:12:16 beside Opus's, ended at once: "You've
    reached your Fable limit" (HTTP 429; `times-fable.txt`). Probes at 01:39
    and 02:09 found the limit still in place, and usage credits are not spent
    without the SSE's permission.
  - So two fresh Opus 5.5 sessions reviewed instead. Reviewer 1 read the
    entries A to I, from 01:12:16 to 01:37:15 (187 turns). Reviewer 2 read
    them I to A, from 01:40:53 to 02:08:15 (182 turns).
  - Fable reviews when its limit allows, as a third report, with the same
    folder, prompt, and tools.
  - During the reviews the machine ran the measured exercises' builds; nothing
    was timed.
- **Reports:**
  - `report-opus-1.md` and `report-opus-2.md` are the sessions' final
    messages, unedited.
  - `session-opus-1.json` and `session-opus-2.json` are their result records,
    with usage: 139,000 and 142,000 output tokens, $16.09 and $14.64 at API
    list prices, on the subscription.

## The key

| Letter | Entry |
|---|---|
| A | astra-max |
| B | sol-xhigh |
| C | astra-xhigh |
| D | astra-high |
| E | sol-high |
| F | luna-max |
| G | sol-max |
| H | luna-xhigh |
| I | luna-high |

## Scores

Criteria: 1 readability, 2 maintainability, 3 extensibility, each 1 to 10.

| Entry | Reviewer | 1 | 2 | 3 | Sum | Both |
|---|---|---:|---:|---:|---:|---:|
| sol-xhigh (B) | 1 | 7 | 9 | 8 | 24 | **47** |
| | 2 | 8 | 8 | 7 | 23 | |
| astra-high (D) | 1 | 7 | 8 | 8 | 23 | **46** |
| | 2 | 7 | 9 | 7 | 23 | |
| sol-max (G) | 1 | 6 | 9 | 8 | 23 | **44** |
| | 2 | 6 | 8 | 7 | 21 | |
| astra-max (A) | 1 | 6 | 8 | 7 | 21 | **44** |
| | 2 | 7 | 8 | 8 | 23 | |
| astra-xhigh (C) | 1 | 6 | 8 | 7 | 21 | **41** |
| | 2 | 6 | 7 | 7 | 20 | |
| sol-high (E) | 1 | 6 | 7 | 6 | 19 | **39** |
| | 2 | 7 | 7 | 6 | 20 | |
| luna-max (F) | 1 | 5 | 4 | 4 | 13 | **31** |
| | 2 | 7 | 5 | 6 | 18 | |
| luna-xhigh (H) | 1 | 4 | 3 | 3 | 10 | **24** |
| | 2 | 5 | 4 | 5 | 14 | |
| luna-high (I) | 1 | 4 | 4 | 3 | 11 | **24** |
| | 2 | 4 | 5 | 4 | 13 | |

## Ranks

| Rank | Reviewer 1 (read A to I) | Reviewer 2 (read I to A) | Both, by total |
|---:|---|---|---|
| 1 | sol-xhigh, 24 | sol-xhigh, 23 | sol-xhigh, 47 |
| 2 | sol-max, 23 | astra-high, 23 | astra-high, 46 |
| 3 | astra-high, 23 | astra-max, 23 | sol-max and astra-max, 44 |
| 4 | astra-max, 21 | sol-max, 21 | |
| 5 | astra-xhigh, 21 | sol-high, 20 | astra-xhigh, 41 |
| 6 | sol-high, 19 | astra-xhigh, 20 | sol-high, 39 |
| 7 | luna-max, 13 | luna-max, 18 | luna-max, 31 |
| 8 | luna-high, 11 | luna-xhigh, 14 | luna-xhigh and luna-high, 24 |
| 9 | luna-xhigh, 10 | luna-high, 13 | |

Each reviewer broke its own ties:
- Reviewer 1 put sol-max above astra-high for its tests: every kernel
  compared at every cell count, and a built-in test of adding a reaction and
  changing a coefficient.
- Reviewer 2 broke ties by the criteria in order: readability, then
  maintainability.

## Where they agree

- **sol-xhigh first.** Reviewer 1 cited a readable generator (named formulas,
  explicit schedules) behind a clean descriptor API, with the most complete
  verification. Reviewer 2 cited the most readable and best-explained
  generator, robust and tested input validation, and strict tests with exact
  solver statistics.
- **The same four on top, and the Luna entries at the bottom.** Both put
  sol-xhigh, astra-high, sol-max, and astra-max in their top four, and
  astra-xhigh and sol-high fifth and sixth in some order. They put luna-max
  seventh and luna-xhigh and luna-high last.
  - The six Astra and Sol entries scored 19 to 24, the three Luna entries 10
    to 18.
  - On the Luna entries, both find weak or near-vacuous tests of other
    mechanisms. Reviewer 1 also calls all three generators single-file or
    monolithic; reviewer 2 says that of luna-high's alone.
  - Both found that the Luna readers leave Arrhenius A unconverted from
    molecule cm-3, against MICM's convention. Reviewer 2 adds that their
    other factors fit only the tutorial's Troe and ternary reactions.
- **Weaknesses all nine share:**
  - **Third body and Ea:** no reader recognizes `"tracer type": "THIRD_BODY"`,
    so Chapman's M would become an integrated species. None reads `"Ea"`; it
    is zero throughout Carbon Bond 5, but a nonzero value would be silently
    ignored.
  - **Label namespaces:** none prefixes photolysis and emission labels with
    MICM's `PHOTO.` and `EMIS.`, so Carbon Bond 5's `NO2` and `ALD2` labels
    collide. By reviewer 1's count, five entries merge them silently into one
    parameter, sol-xhigh, astra-xhigh, and astra-high keep two slots under one
    name, and sol-high rejects the file.
  - **Annotation:** no generated C names the reactions or species it
    implements.
  - **A shared reader:** in every entry the generic library and the generator
    share the reader, so comparing generic and generated results cannot catch
    a reader error.

## Where they differ

- **sol-max:** second for reviewer 1, for its tests. Fourth for reviewer 2,
  which judged its density and many options the hardest generator of the
  strong group to follow, with templates that repeat the control flow.
- **astra-max:** reviewer 2 found it the easiest path to a new mechanism and
  rated exercise (c) easy, the only "easy" either reviewer gave for (c).
  Reviewer 1 rated it moderate and noted terse code (numeric rate codes,
  positional parameters).
- **Exercise (c) for the Luna entries:** hard for reviewer 1 (their readers
  misread general mechanisms), moderate for reviewer 2.
- **luna-max:** 13 against 18, the widest gap. Reviewer 2 credited a compact,
  readable generator and a clean optional-kernel interface.

## Against the measured exercises

- **(a) add a reaction and (b) change a rate parameter:** both reviewers rate
  them easy for all nine. Measured, all nine pass both.
- **(c) another mechanism:** measured, every generator took Chapman and
  Carbon Bond 5 through the entry's own benchmark programs, within 3.4e-14 of
  MICM. The reviewers rated (c) moderate or hard. Both findings stand, because
  they test different things:
  - The measured exercise gave each entry a mechanism resolved in advance
    (units, Ea, third body, label namespaces, as MICM's former parser
    resolved them), written by a converter into the entry's own format. It
    tested the generators.
  - The reviewers judged the path from MUSICA's v0 files through each entry's
    own reader, and found the gaps above.
  - So every generator handles a new mechanism given correct input, but no
    entry's reader takes MUSICA's v0 files as MICM would.
- **A fallback seen both ways:** reviewer 1 read that astra-high's and
  sol-high's generated kernels serve only full 32-cell groups, falling back to
  generic kernels otherwise, so their tests exercise the generated code only
  at 32 cells. The measured exercise found the same fallback in its 4-cell
  runs.

The reports' per-entry sections give the reasons for every score, the
exercises on paper, and each entry's main strengths and weaknesses.
