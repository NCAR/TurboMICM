# Tier 2: numerical algorithm design

Tier 1 asked agents to port MICM's CPU solver to portable C and make it fast while
keeping MICM's numerical methods. Tier 2 frees the methods. Each agent starts
from its own Tier 1 solution, researches the literature, prototypes and scans
method structures, and implements the fastest methods it can justify as accurate
enough for an atmospheric model.

Status: designed and approved on 2026-10-15; tooling built; not yet launched.
The prompt template is `prompt.v1`. This directory mirrors the committed files of
the Tier 2 run root, `ODE/benchmarks/tier2`; keep the two identical.

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
- `<root>`: the Tier 2 run root, planned as `/Users/fillmore/NSF-NCAR/ODE/benchmarks/tier2`;
- `<tier1_commit>`: the run's final Tier 1 commit, for example `001fb64` for
  sol-max.

## Tooling

| Path | Purpose | Status |
|---|---|---|
| `prompts/make_prompts.py` | fills the template; `<tier1_commit>` is the Tier 1 run's last commit | built; sol-xhigh and sol-max prompts made |
| `admin/setup_runs.py` | clones each Tier 1 run repository into the Tier 2 root (no shared files, no remote), writes PROMPT.md, links `micm/` and `reference/` | built and tested; sol-xhigh and sol-max set up |
| `scenarios/micm_scenario.cpp` | MICM on a scenario: f(t) photolysis at the start of each step, any MICM Rosenbrock method and tolerances, snapshots | built; bit-identical to `micm_ts1_ref` on Tier 1's scenario |
| `admin/make_scenarios.sh`, `scenarios/make_scenario.py` | the six disclosed scenarios, and 10,000-cell versions for timing | built; regenerates the disclosed set byte for byte |
| `admin/make_references.sh` | references: MICM's six-stage Rosenbrock at 1e-10 relative and 1e-23 absolute, checked against 1e-9 and 1e-22 | tier1-30s done (converged to 10.8 digits; MICM's default reaches 4.6); the other five suspended while astra-high runs |
| `scenarios/score.py` | significant digits overall and per species, and the gate | built |
| `scenarios/README.md` | the scenario, output, and scoring formats, for the agents | written |

Still to do before launch:
- the five remaining references, on an idle machine;
- the undisclosed scenarios, made the same way and kept outside the run root;
- scenario timing in `admin/time_entries.py`, then MICM's time on every
  scenario, for context.

## Launch

    cd ~/NSF-NCAR/ODE/benchmarks/tier2/<run>
    codex --search -m <model> -c model_reasoning_effort=<effort> -c features.memories=false

then

    /goal read PROMPT.md and work until done, commit
