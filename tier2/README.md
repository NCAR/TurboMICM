# Tier 2: numerical algorithm design

Tier 1 asked agents to port MICM's CPU solver to portable C and make it fast while
keeping MICM's numerical methods. Tier 2 frees the methods. Each agent starts
from its own Tier 1 solution, researches the literature, prototypes and scans
method structures, and implements the fastest methods it can justify as accurate
enough for an atmospheric model.

Status: designed and approved on 2026-10-15; not yet launched. The prompt template
is `prompt.v1`.

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

## To build before launch

1. The Tier 2 run root, `ODE/benchmarks/tier2`, with `micm/` and `reference/`
   linked to Tier 1's, so the entries' Makefiles still work, and a new
   `scenarios/`.
2. A MICM scenario driver: photolysis rates updated at the start of each step
   from f(t), tight-tolerance settings, and snapshot output. Check it against
   `micm_ts1_ref` on Tier 1's scenario.
3. The scenarios:
   - Tier 1's benchmark;
   - 24-hour day-night runs, with f(t) from the solar zenith angle, at 5, 10,
     and 30 minute steps;
   - clean and polluted cell sets;
   - undisclosed variants: other seeds, latitudes, seasons, and 48-hour runs.
4. High-accuracy references: MICM's Rosenbrock at very tight tolerances, checked
   for convergence.
5. `scenarios/score.py` (significant digits overall and per species, the NaN and
   negativity checks) and `scenarios/README.md` (the scenario file format and
   the ts1_snapshots.csv format).
6. Prompt generation and setup for Tier 2: the `<tier1_commit>` placeholder, and
   run directories cloned from `runs/<model>/<effort>`.
7. Scenario support in the timing harness.
8. MICM's own accuracy and time on every scenario, for context.

## Launch

    cd ~/NSF-NCAR/ODE/benchmarks/tier2/<run>
    codex --search -m <model> -c model_reasoning_effort=<effort> -c features.memories=false

then

    /goal read PROMPT.md and work until done, commit
