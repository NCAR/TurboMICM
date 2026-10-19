# Tier 3: numerical algorithm design

Numbered Tier 2, in `tier2/`, until 2026-10-16, when specialization for a frozen
mechanism became Tier 2 and this tier became Tier 3.

Tier 1 asked agents to port MICM's CPU solver to portable C and make it fast while
keeping MICM's numerical methods. Tier 3 frees the methods. Each agent starts
from its own Tier 1 solution, researches the literature, prototypes and scans
method structures, and implements the fastest methods it can justify as accurate
enough for an atmospheric model.

Status: designed and approved on 2026-10-15; tooling built; launched on 2026-10-18,
one run at a time in the order of the Tier 1 runs (sol-xhigh, sol-max,
astra-xhigh, astra-max; see `admin/runs.md`). All four ended and were verified
that day, and the rest followed on the same laptop by the morning of
2026-10-19: MICM's baselines, the blind review
(`results/2026-10-18-tier3-blind-review/`), the timing session, and the
ranking (`results/2026-10-19-tier3-ranking/`; the wiki's Tier-3-Results has
the tables). Luna's arm is open; it would need new held-out scenarios, since
`heldout/` is public.
The prompt template is `prompt.v1`. This directory mirrors the committed files of
the Tier 3 run root, `ODE/benchmarks/tier3`; keep the two identical, except for
this README, `heldout/` (the held-out scenarios, published after the last
run), and `review/` (the blind review's instructions).

## Plan

- **Arms:** Sol and Astra at xhigh and max effort, a 2 x 2 matrix, each starting
  from its own Tier 1 run. Luna gets one run at max effort if its Tier 1 runs pass
  comfortably. The top two Claude models follow at the same two effort levels,
  after their Tier 1 runs.
- **Methods:** any mechanism-generic numerical method. Mechanism reduction and
  learned surrogates are excluded. MICM's methods and the Tier 1 gate are kept.
- **Research:**
  - Web search is allowed for the literature (`codex --search`). No code may be
    copied, and every source used is cited in REFERENCES.md.
  - Methods may be prototyped in Python (numpy and scipy, or driving the agent's
    C pieces through ctypes), measuring accuracy and operation counts. Speed
    claims come only from C, and a test checks each C method against its
    prototype.
  - Structure scans are logged in scans/ and summarized in SCANS.md.
- **Accuracy:** each agent decides what is good enough and justifies it in
  ACCURACY.md. The SSE measures every entry with one yardstick:
  - significant digits (Sandu et al., 1997) against high-accuracy MICM
    references, over species at or above 1e-14 mol m-3;
  - a gate of no NaN or infinity and no value below -1e-20 mol m-3.

  The threshold covers 177 of TS1's 209 species in most cells, including OH, HO2,
  NO, NO2, NO3, O3, and N2O5. MICM itself produces no negative values on TS1.
- **Ranking:** an entry's result at level L is its fastest declared operating
  point that reaches L digits overall, and L - 1 digits for every species, on
  every scenario. Entries are ranked at L = 2 and L = 3, on 10,000 cells, median
  of 5 runs, idle machine.
- **Time:** no limit. Overnight runs are expected; an agent continues while it
  is making progress, and the SSE may stop a run that stalls.

## Prompt placeholders

`<model>`, `<run>`, `<micm_commit>`, as in Tier 1, plus:
- `<root>`: the Tier 3 run root, planned as `/Users/fillmore/NSF-NCAR/ODE/benchmarks/tier3`;
- `<tier1_commit>`: the run's final Tier 1 commit, for example `001fb64` for
  sol-max.

## Tooling

| Path | Purpose | Status |
|---|---|---|
| `prompts/make_prompts.py` | fills the template; `<tier1_commit>` is the Tier 1 run's last commit | built; all four prompts made |
| `admin/setup_runs.py` | clones each Tier 1 run repository into the Tier 3 root (no shared files, no remote), writes PROMPT.md, links `micm/` and `reference/` | built and tested; all four runs set up |
| `scenarios/micm_scenario.cpp` | MICM on a scenario: f(t) photolysis at the start of each step, any MICM Rosenbrock method and tolerances, snapshots | built; bit-identical to `micm_ts1_ref` on Tier 1's scenario |
| `admin/make_scenarios.sh`, `scenarios/make_scenario.py` | the six disclosed scenarios, and 10,000-cell versions for timing | built; regenerates the disclosed set byte for byte |
| `admin/make_references.sh` | references: MICM's six-stage Rosenbrock at 1e-10 relative and 1e-23 absolute, checked against 1e-9 and 1e-22 | all six done on 2026-10-15 (10 threads, 16:31–17:00): converged to 10.18–10.79 digits overall, no species below 8.86 |
| `scenarios/score.py` | significant digits overall and per species, and the gate | built |
| `scenarios/README.md` | the scenario, output, and scoring formats, for the agents | written |

Done after the launch (the agents needed neither):
- the undisclosed scenarios, made the same way and kept outside the run root:
  five, approved by the SSE and fixed on 2026-10-18 at 01:50, before any result
  (the script's SHA-256 is in `admin/runs.md`); their files and references were
  made on 2026-10-18 between sol-xhigh and sol-max, converged to 10.16–10.69
  digits overall, and published in `heldout/` when the last run ended;
- a test of the scenario mode of `admin/time_entries.py`, on 2026-10-18.

## How the ranking was made

The rules, fixed before the runs (in `HANDOFF.md` until the ranking was
published on 2026-10-19):
- **Timing:** one interleaved session on one idle machine on AC power under
  `caffeinate`: all 19 declared points and MICM's baselines on the 10,000-cell
  versions of the 11 scenarios, 5 rounds in rotating order, all cores (10 on
  the M5). It ran from 21:19 on 2026-10-18 to 05:46 on 2026-10-19.
- **Accuracy for the ranking** from the same session: each 10,000-cell output's
  first 1,000 cells scored against the 1,000-cell reference, the gate applied
  to every cell.
- **Level L:** L digits overall and L - 1 for every species, on every scenario
  (`scenarios/score.py --level`).
- **Rank:** for each entry and level, its fastest declared point that reaches
  the level on all 11 scenarios, by the sum over the scenarios of the point's
  median elapsed time; the disclosed and held-out sums shown beside it.
- **MICM's baselines:** MICM's default three-stage Rosenbrock with 32-cell
  groups, one tolerance setting per level, the fastest of an 18-setting scan
  that reaches the level on all 11 scenarios at 1,000 cells, confirmed in the
  session (relative 1e-2 with absolute 1e-15 for level 2, 3e-4 with 1e-15 for
  level 3).
- **Blind review:** Opus and Fable at max effort, in fresh sessions, on
  anonymized copies with random letters; the key private until both reports
  were in; run while nothing was timed.
- **Tools:** `admin/verify_tier3.py`, `admin/time_entries.py --scenario`,
  `admin/time_tier3.sh`, `admin/rank_tier3.py`, `admin/anonymize_tier3.py`,
  and `review/`.

## MICM's default on the scenarios

MICM's default three-stage Rosenbrock (relative 1e-6, absolute 1e-3 mol m-3 for
every species; MUSICA's TS1 and TS4 configurations set no absolute tolerances)
is accurate on Tier 1's benchmark but not over a day:

| Scenario | Reference convergence | MICM default |
|---|---:|---|
| tier1-30s | 10.79 | 4.63 digits |
| diurnal-30min | 10.31 | fails: NaN in step 1 |
| diurnal-10min | 10.22 | fails: step size too small in step 3 |
| diurnal-5min | 10.18 | 1.21 digits (worst species 0.21) |
| polluted-10min | 10.25 | fails: NaN in step 1 |
| clean-10min | 10.20 | 2.96 digits (worst species 2.16) |

The default takes 9 internal steps per scenario step whatever the step length,
which suggests its error control is inactive: an absolute tolerance of 1e-3 mol
m-3 (about 24 ppm) is above every trace gas. Tightening only the absolute
tolerance fixes it. These are single runs on `diurnal-30min`, 1,000 cells, 10
threads:

| Three-stage Rosenbrock | Time | Overall digits | Worst species |
|---|---:|---:|---|
| default | — | fails | NaN in step 1 |
| absolute 1e-12 | 1.6 s | 2.76 | ISOP 1.50 |
| absolute 1e-15 | 4.5 s | 4.77 | ISOP 3.63 |
| relative 1e-3, absolute 1e-15 | 2.0 s | 2.97 | ISOP 1.71 |
| relative 1e-2, absolute 1e-15 | 1.6 s | 2.49 | ISOP 1.24 |

The agents' scenarios README gives the first table, which matches what
`references.log` already shows them. The tolerance experiment stays here: choosing
tolerances is part of their research.

## Launch

Claude Code launches each run from TurboMICM with

    admin/launch_codex_run.sh --tier3 <run> <model> <effort>

which checks that the run directory is an untouched clone at its Tier 1 commit,
then types into a detached tmux session (`tmux -L turbomicm attach -t tier3-<run>`)

    cd ~/NSF-NCAR/ODE/benchmarks/tier3/<run>
    codex --search -m <model> -c model_reasoning_effort=<effort> -c features.memories=false

and then

    /goal read PROMPT.md and work until done, commit
