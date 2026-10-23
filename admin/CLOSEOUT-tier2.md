# Tier 2 closeout

**State on 2026-10-23 at 07:25.** Tier 2 is done except one item: Fable's
blind review, which waits for Fable's usage reset at 22:00 (the SSE's word).
Everything else is verified, written up, committed, and pushed.
- **This repository:** `main`; the nine branches `tier2/<model>/<effort>`; the
  nine transcripts on `transcripts`.
- **The wiki:** Tier-2-Results, Home, Running-the-Experiment, and
  Role-of-Claude-in-the-Design.

The SSE will resume the orchestrating session to finish it. This file says
what is done, what is running, and what to do when Fable's report arrives.

## What is done

| Item | Where |
|---|---|
| Nine runs, 2026-10-22 09:58 to 2026-10-23 00:09, in the SSE's order | `admin/runs.md` (Tier 2), branches `tier2/<model>/<effort>`, `transcripts/tier2/<run>.jsonl.gz` |
| Verification of every entry, both compilers: the 680 cases, both TS1 gates (seed 2026 and the undisclosed input), `make test-specialized`, ten day-night scenarios, the audit | `results/2026-10-22-tier2-<run>-verification`, `results/2026-10-23-tier2-luna-high-verification` |
| The held-out scenarios and the ranking seed, published after the last run | `tier2/heldout/` |
| The ranking at 18 threads (00:34 to 01:11) and at one thread (02:11 to 04:30), 255 runs, all passing the gate | `results/2026-10-23-tier2-ranking` |
| The review's three exercises, measured on all nine | `results/2026-10-23-tier2-exercises` |
| The blind review by two fresh Opus sessions, decoded | `results/2026-10-23-tier2-blind-review` |
| Plan, status, and tooling | `tier2/README.md`, CLAUDE.md (Tier 2) |

The results, in brief:

| Entry | Speed rank, 18 threads (s) | Speed rank, 1 thread (s) | Blind review, two Opus reviewers (of 60) | Exercises |
|---|---|---|---|---|
| astra-xhigh | 1 (2.05) | 1 (26.5) | 5 (41) | all pass |
| sol-max | 2 (2.40) | 2 (30.5) | 3= (44) | all pass |
| sol-xhigh | 3 (2.71) | 4 (35.3) | 1 (47) | all pass |
| sol-high | 4 (2.80) | 3 (34.7) | 6 (39) | all pass |
| astra-high | 5= (3.16) | 6 (41.0) | 2 (46) | all pass |
| astra-max | 5= (3.22) | 5 (36.5) | 3= (44) | all pass |
| luna-max | 7 (5.38) | 8 (52.2) | 7 (31) | all pass |
| luna-high | 8 (7.57) | 9 (89.8) | 8= (24) | all pass |
| luna-xhigh | 9 (10.06) | 7 (47.1) | 8= (24) | all pass |

## Deviations from the design, all recorded

- **The blind review's second reviewer.** The design names Opus and Fable.
  Fable stopped at its usage limit at 01:12 and stayed there, so two fresh
  Opus sessions reviewed instead (A to I, then I to A). Fable's is to be a
  third report.
- **Exercise (c)'s first input.** It had 4 identical cells, which leave
  astra-high's and sol-high's generated kernels unused, since those kernels
  serve only full 32-cell groups. sol-high's agent found this. A 64-cell input
  of distinct cells was added, and every entry passes it.
- **The exercise agents' setup.** Their instructions wrongly said `make env`
  uses an existing environment: it creates one, downloading from conda-forge.
  Three agents let it download; the others reused an environment the same
  entry's `make env` had built.

## What is running

- **The watcher:** tmux `tier2-review-fable` on socket `turbomicm`, running
  `~/NSF-NCAR/tier2-work/blind-review/fable-when-available.sh`.
  - It waits until 22:00:30, checking the clock every minute. Then it probes
    Fable every 5 minutes; a refused probe costs nothing.
  - Its log is `~/NSF-NCAR/tier2-work/blind-review/fable-probe.log`.
- **When Fable answers,** the watcher runs `tier2/review/run-review.sh fable
  "I H G F E D C B A" REVIEW OUT fable`. That is the same anonymized folder,
  prompt, and read-only tools as the two Opus reviews, at max effort.
  - It writes `report-fable.json`, `report-fable.err`, and then
    `report-fable.done` in `~/NSF-NCAR/tier2-work/blind-review/`.
  - The Opus reviews took 25 and 27 minutes.
- **The laptop must be awake** around 22:00 for the review to run. If it
  sleeps, the watcher carries on when it wakes.
- **Nothing else runs.** No Codex session, no timing, no agents.

Check it with:

    tail -3 ~/NSF-NCAR/tier2-work/blind-review/fable-probe.log
    ~/homebrew/bin/tmux -L turbomicm ls
    ls ~/NSF-NCAR/tier2-work/blind-review/report-fable.*

## When Fable's report is in

1. **Confirm it ran.** `report-fable.done` exists, `report-fable.err` ends
   with `exit 0`, and `report-fable.json` has `"is_error": false` and a JSON
   block at the end of its `result`.
   - If Fable was refused again, the watcher's log says so; rerun the
     watcher, or ask the SSE.
2. **Check that it stayed blind.** Run
   `python3 admin/check_review_session.py <session_id from report-fable.json>`.
   It must report 0 paths outside the folder. The review folder's parent holds
   `key.txt`, and the key is now public in
   `results/2026-10-23-tier2-blind-review/key.txt`, so this check matters.
3. **Decode it.** Run `python3 admin/review_tables.py`, which prints all three
   reviewers' scores, totals of 90, ranks, and exercises, decoded.
4. **Save it** to `results/2026-10-23-tier2-blind-review/`:
   - `report-fable.md`, the session's final message, unedited;
   - `session-fable.json`, the whole record;
   - `times-fable.txt`, from `report-fable.err`.

   First rename the existing `times-fable.txt`, from the failed 01:12
   attempt, to `times-fable-attempt1.txt`.
5. **Update that folder's README.md:**
   - Procedure: Fable's session, its times, turns, and usage, and the
     transcript check.
   - Scores: a third row per entry, and totals of 90.
   - Ranks: a third column.
   - Where the three agree and differ, and whether Fable changes the order.

   Keep the two Opus reviewers' results as they are.
6. **Update the wiki:**
   - Tier-2-Results, "Blind review": a reviewer-3 column, the total of 90, and
     the deviation note (Fable is now the third reviewer, the design's two
     plus one).
   - Tier-2-Results, the status line and the summary table.
   - Home's Tier 2 status.
7. **Update the records:**
   - `admin/runs.md`, Tier 2 "blind review and exercises": the Fable session.
   - CLAUDE.md's Tier 2 status: remove the open item.
   - `tier2/README.md` status: Tier 2 complete.
8. **Close the watcher** (`tmux -L turbomicm kill-session -t
   tier2-review-fable`, if it has not exited), then commit and push
   TurboMICM and the wiki. ODE has no Tier 2 change to push: the mirror
   matches, `diff` of the committed files.

## Left in place, removable later

These are all in `~/NSF-NCAR/tier2-work`, outside the run root, and all 46 are
registered as worktrees of this repository:
- `verify/<run>/{clang,gcc-15}/<run>`: the verification builds, whose clang
  binaries were timed (11 GB);
- `tier1-builds/<run>`: the Tier 1 commits built here (1.9 GB);
- `variants/astra-high-nocache`: astra-high without its rate cache, from its
  verification;
- `exercises/ab/<run>/<run>` and `exercises/c/<run>/<run>`: the exercises
  (3.1 GB with their outputs).

Keep `verify/` and `tier1-builds/` while a retiming is possible. The exercises'
results are copied into `results/2026-10-23-tier2-exercises`. When the
worktrees go, remove them with `git worktree remove`, then `git worktree
prune`.
