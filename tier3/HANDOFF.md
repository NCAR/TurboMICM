# Tier 3 handoff, 2026-10-18

**State on 2026-10-18 at 21:30:** the SSE decided at 19:32 to finish on the
same laptop, and steps 1 to 6 below are done or running there: the inputs
checked (step 2), the entries re-verified with identical results into
`~/NSF-NCAR/tier3-work/verify/` (step 3), the MICM scan run and its baselines
chosen (step 4; `micm/L2` rtol 1e-2 atol 1e-15, `micm/L3` 3e-4 and 1e-15,
plus `micm/L3b` 2e-4 and 1e-15 as insurance for L3's thin isoprene margin),
the programs file written (step 5; `~/NSF-NCAR/tier3-work/programs.txt`, 22
programs), the blind review done (step 8; `results/2026-10-18-tier3-blind-review/`),
and the timing session started at 21:19 in tmux session `tier3-timing` on
socket `turbomicm` (step 6; `~/NSF-NCAR/tier3-work/run-timing.sh`, output in
`~/NSF-NCAR/tier3-work/timing/`, a load audit beside it). Remaining: step 7
(rank) when it ends, and step 9 (publish). The text below is the original
handoff, kept for the record.

All four Tier 3 runs are done: collected, archived, and verified on the SSE's
Apple M5 laptop. The rest was to move to another machine:
- the MICM baselines;
- the timing session;
- the ranking;
- the blind review;
- publishing.

Everything needed is on GitHub. Regenerate results on the new machine rather than
copying them. Read CLAUDE.md first; this file adds only what Tier 3 still needs.

## Where things stand

| Run | Branch, final commit | Points | Verified (1,000 cells, 6 disclosed + 5 held-out scenarios) |
|---|---|---|---|
| sol-xhigh | `tier3/sol/xhigh`, `63d246f` | 5 | all pass the gate; each reaches its declared level on all 11 |
| sol-max | `tier3/sol/max`, `c122d12` | 5 | the same |
| astra-xhigh | `tier3/astra/xhigh`, `db44838` | 5 | the same |
| astra-max | `tier3/astra/max`, `3cab602` | 4 | all pass the gate; fast3, declared for level 3, reaches it on 9 of 11 (misses two held-out scenarios) |

- **Run records:** each run's facts are in `admin/runs.md`, its verification in
  `results/2026-10-18-tier3-<run>-verification/`, and its redacted transcript at
  `transcripts/tier3/<run>.jsonl.gz` on the `transcripts` branch.
- **Held-out scenarios:** published in `tier3/heldout/` once the last run had
  ended. That folder has the generator, a driver that makes all the timing
  inputs, and checksums.
- **Tools:**
  - `admin/verify_tier3.py`.
  - `admin/time_entries.py --scenario`: score once, by SHA-256. It was tested on
    2026-10-18: 2 programs on the 6 disclosed scenarios, 2 rounds.
  - `admin/time_tier3.sh` and the MICM scan list `admin/tier3_micm_scan.txt`.
  - `admin/rank_tier3.py`, tested on synthetic timings.
  - `admin/anonymize_tier3.py` and the reviewers' instructions in `tier3/review/`.

## Decisions already made

- **Timing:** one interleaved session on one machine, which must be:
  - idle: no agent running and no builds;
  - on AC power;
  - kept awake with `caffeinate`.

  It covers all 19 declared points plus MICM's two baselines, run on the
  10,000-cell versions of the 11 scenarios, 5 rounds each, at all cores. All
  cores is the prompt's default: 10 on the M5; use the new machine's count.
- **Accuracy for the ranking** comes from the same session:
  - each 10,000-cell output's first 1,000 cells are scored against the
    1,000-cell reference;
  - the gate applies to every cell.
- **Level L:** L digits overall and L - 1 for every species, on every scenario
  (`score.py --level`).
- **Rank:** for each entry and level, its fastest declared point that reaches the
  level on all 11 scenarios. "Fastest" is the sum over the 11 scenarios of the
  point's median elapsed_ms. The sums over the disclosed and held-out sets are
  shown too. This rule was put to the SSE, who has chosen no other.
- **MICM baselines:**
  - MICM's default three-stage Rosenbrock with 32-cell groups.
  - One setting per level, used on all 11 scenarios, as the entries' points
    must be: the fastest setting in the scan that reaches the level on all 11 at
    1,000 cells.
  - Confirm it in the timing session.
- **Blind review:**
  - Opus and Fable at max effort, in fresh sessions, on anonymized copies with
    random letters.
  - Keep the letter key private until both reports are in.
  - Run it while nothing is being timed.

## Steps on the new machine

1. **Set up** as CLAUDE.md describes, at the same paths under `$HOME`:
   - Clone `NCAR/TurboMICM` to `~/NSF-NCAR/TurboMICM` and its wiki to
     `~/NSF-NCAR/TurboMICM.wiki`.
   - Clone ODE to `~/NSF-NCAR/ODE`, on branch `micm-c-port-experiment`.
   - Build MICM e15898ce in `ODE/benchmarks/micm`: Release, double precision,
     OpenMP; see "MICM build" in Pitfalls.
   - Install Homebrew's libomp and gcc-15, and conda (the entries' `make env`).
   - Set `SDKROOT` for Apple clang.
   - Make the run root's links, which ODE does not track:
     `cd ~/NSF-NCAR/ODE/benchmarks/tier3 && ln -sn ../micm micm && ln -sn ../reference reference`.
   - Keep the working files in `~/NSF-NCAR/tier3-work/`. Don't use `/private/tmp`:
     the timing session needs its binaries for hours.

2. **Timing inputs** (about 30 minutes):
   - Run `zsh tier3/heldout/make_data.sh`.
   - Then check: `cd ~/NSF-NCAR/.tier3-heldout && shasum -a 256 -c ~/NSF-NCAR/TurboMICM/tier3/heldout/SHA256SUMS`.
   - The scenario files should match exactly. The references may differ in the
     last digits.

3. **Re-verify the four entries** (about 5 minutes each). This also builds the
   clang binaries the timing uses, and shows whether the new machine reproduces
   the committed reports, digit for digit with the same compiler versions.

       W=~/NSF-NCAR/tier3-work/verify
       for spec in sol-xhigh:tier3/sol/xhigh sol-max:tier3/sol/max astra-xhigh:tier3/astra/xhigh astra-max:tier3/astra/max; do
         r=${spec%%:*} b=${spec#*:} p=results/2026-10-18-tier3-$r-verification/points.txt
         python3 admin/verify_tier3.py $r $b $W/$r --points $p
         cp $W/$r/report.md $W/$r/report-disclosed.md
         python3 admin/verify_tier3.py $r $b $W/$r --points $p --scenarios ~/NSF-NCAR/.tier3-heldout/heldout --stages points
       done

4. **MICM scan:**
   - Run `admin/time_tier3.sh ~/NSF-NCAR/tier3-work/micm-scan admin/tier3_micm_scan.txt disclosed,heldout --threads N --rounds 1 --pause 1`.
   - Then `admin/rank_tier3.py ~/NSF-NCAR/tier3-work/micm-scan --threads N`.
   - There, each setting counts as its own entry. Pick the fastest that reaches
     level 2 on all 11, and the fastest that reaches level 3.

5. **Programs file:** one line per declared point, named ENTRY/POINT, plus the
   baselines as `micm/L2` and `micm/L3`. The entry name `micm` gives the ranking
   its speedup column.

       W=~/NSF-NCAR/tier3-work/verify
       for r in sol-xhigh sol-max astra-xhigh astra-max; do
         grep -v '^#' results/2026-10-18-tier3-$r-verification/points.txt |
           while IFS=$'\t' read -r name opts; do echo "$r/$name=$W/$r/clang/$r/bin/ts1_scenario $opts"; done
       done > ~/NSF-NCAR/tier3-work/programs.txt
       echo 'micm/L2=$HOME/NSF-NCAR/.tier3-heldout/bin/micm_scenario --group-size 32 --rtol R2 --atol A2' >> ~/NSF-NCAR/tier3-work/programs.txt
       echo 'micm/L3=$HOME/NSF-NCAR/.tier3-heldout/bin/micm_scenario --group-size 32 --rtol R3 --atol A3' >> ~/NSF-NCAR/tier3-work/programs.txt

6. **Timing session** (about 8 hours on the M5, 1,155 runs):
   - Run `caffeinate -i admin/time_tier3.sh ~/NSF-NCAR/tier3-work/timing ~/NSF-NCAR/tier3-work/programs.txt disclosed-10000,heldout-10000 --threads N --rounds 5`.
   - Each output is 200–400 MB, and `time_entries.py` deletes it after scoring.

7. **Rank:**
   - Run `admin/rank_tier3.py ~/NSF-NCAR/tier3-work/timing --threads N`.
   - Rescore exactly, with `score.py --level`, any point it flags.
   - Keep the ranking and each scenario's `timings.csv`, `summary.md`, and
     `environment.txt` in `results/<date>-tier3-ranking/`.

8. **Blind review:**
   - Draw a random order of the four runs for the letters A–D, and record the
     key privately.
   - Make the copies: `admin/anonymize_tier3.py REVIEW sol-xhigh=tier3/sol/xhigh=<letter> ...`.
     It exits 1 when it finds possible leftovers; check each one. A PNG false
     positive is known.
   - Copy into REVIEW: `tier3/review/REVIEW.md` and `PROMPT.md`;
     `tier3/scenarios/README.md` and `score.py` into `scenarios/`;
     `reference/README.md` into `reference/`.
   - Start two fresh subagents at effort max: model opus reading A to D, and
     model fable reading D to A. Each follows REVIEW.md, reads only REVIEW, and
     ends with the JSON block.
   - Put both reports and the key in `results/<date>-tier3-blind-review/`.

9. **Publish:**
   - The ranking's README.
   - The wiki's Tier-3-Results: each entry's fastest point at 2 and 3 digits in
     the matrix, the held-out scenarios named, and their tables. Also Home.
   - `admin/runs.md`.
   - CLAUDE.md: the Tier 3 section and Open items.
   - `tier3/README.md`; delete this file.

## Caveats

- **astra-max's fast3:** if the 10,000-cell scores agree, its level-3 rank falls
  to balanced, its recommended point.
- **Narrow margins:** two level-2 points sit close to the per-species threshold
  of 1. sol-max's fast has 1.06 disclosed and 1.09 held out; astra-xhigh's fast
  has 1.035 and 1.13. At 10,000 cells, the last cell group within the first
  1,000 cells is full rather than partial, so scores can move slightly.
- **Rounding:** `timings.csv` records digits rounded to three decimals, and
  `rank_tier3.py` flags values that round exactly to a threshold.
- **Tuning to the M5:** the entries tuned threads, grouping, and scheduling on
  the M5 (10 cores), so another chip may shift their relative times. Name the
  machine in the results, and compare times only within one session.
- **Held-out scenarios are now public:** any later Tier 3 run, such as Luna's,
  would need new ones.
