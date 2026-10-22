# Tier 2: specialization for a frozen mechanism

Tier 1 asked agents to port MICM's CPU solver to portable C and make it fast while
keeping MICM's numerical methods, in a library for any mechanism. Tier 2 keeps the
methods and specializes the solver to one frozen mechanism, TS1, through a code
generator that the agent keeps and that works for any mechanism. Tier 3 (in
`tier3/`) frees the methods.

Status: adopted by the SSE on 2026-10-16, with the prompt's choices settled the
same day (`prompt.v1`). On 2026-10-22 the SSE chose the arms, all nine Tier 1
entrants, and their order, and the runs moved to a second machine (below);
`prompt.v2` is the template for that machine. The run log is
`admin/runs.md`, under "Tier 2"; the wiki has Tier-2-Design, Tier-2-Prompt, and
Tier-2-Results.

This directory mirrors the committed files of the Tier 2 run root,
`ODE/benchmarks/tier2`; keep the two identical (`diff -r`), except for this
README and, after the last run, `heldout/`.

## The machine

Tier 1 and Tier 3 ran on the SSE's other Mac, an Apple M5 laptop (Gemini.local:
4 Super and 6 Efficiency cores, 16 GB, macOS 27.0.1, Xcode, Homebrew in
`/opt/homebrew`). Tier 2 runs on an Apple M5 Pro laptop (LALAL0726071122: 6
Super and 12 Performance cores, 24 GB, macOS 26.7.1, Command Line Tools 27
without Xcode, Homebrew in `~/homebrew`, installed without administrator
rights). What changes:
- **Paths.** Every Tier 1 Makefile names `/opt/homebrew/opt/libomp`, which does
  not exist here, so each agent must repoint its Makefile; `prompt.v2` says so,
  and gives this machine's paths.
- **Toolchain.** Apple clang 21.0.0 as before; GCC 15.3.0 (Homebrew `gcc@15`,
  built from source here on 2026-10-22) instead of 15.2.0; libomp 23.1.3
  instead of 22.1.8; cmake 4.4.4; conda 26.5.3. Apple clang links here with or
  without `SDKROOT`.
- **MICM** is the same commit, e15898ce, built here on 2026-10-22 with Apple
  clang, Release, double precision, OpenMP: all 67 ctest tests pass, the case
  list regenerated from this build is identical to `reference/micm_cases.txt`
  (680 cases), and `micm_ts1_ref` reproduces
  `reference/micm_ts1_final_100_seed1.csv` byte for byte.
- **Timing.** Times from this machine are never compared with Tier 1's. Every
  Tier 2 entry is timed here, with its own Tier 1 benchmark and MICM in the same
  session, at 18 threads and at 1.
- **Codex** is the same version, 0.163.0, and the same subscription; its
  `config.toml` here does not grant full access, which every contest run has had
  (approval `never`, sandbox `danger-full-access`), so the launch must set it.

## Plan

- **Arms:** all nine Tier 1 entrants, one run each, each from its own Tier 1
  repository at its final commit.
- **Order** (the SSE's, 2026-10-22, "in reverse order of capability"): astra-max,
  astra-xhigh, astra-high, sol-max, sol-xhigh, sol-high, luna-max, luna-xhigh,
  luna-high. One at a time: each is collected, archived, verified, and timed
  before the next launches, and the wiki is updated as results come in.
- **Time:** no limit; the SSE may stop a run that flounders.
- **Web:** none, as in Tier 1.

## Prompt placeholders

`<model>`, `<run>`, `<micm_commit>`, as in Tier 1, plus:
- `<root>`: the Tier 2 run root, `/Users/dfillmor/NSF-NCAR/ODE/benchmarks/tier2`;
- `<tier1_commit>`: the run's final Tier 1 commit, from the TurboMICM branch
  `runs/<model>/<effort>`, for example `0dc7f67` for astra-max.

## Tooling

| Path | Purpose |
|---|---|
| `prompt.v1` | the template as settled on 2026-10-16, for the M5 laptop; never used |
| `prompt.v2` | the template for this machine; the changes are listed in `admin/runs.md` |
| `prompts/make_prompts.py` | fills the template; `<tier1_commit>` comes from TurboMICM's `runs/<model>/<effort>` |
| `admin/setup_runs.py` | copies each Tier 1 history into `ROOT/<run>` (branch `main`, no remote, no tags), writes PROMPT.md, links `micm/` and `reference/` |
| `admin/make_scenarios.sh`, `scenarios/make_scenario.py` | the five disclosed day-night scenarios, Tier 3's five of six that run through a day, byte for byte |
| `admin/make_scenario_data.sh` | MICM's results at the gate's settings (`micm_*.csv`) and the references (`reference_*.csv`), with `micm.log` |
| `scenarios/micm_scenario.cpp` | MICM on a scenario, as in Tier 3 (one comment line differs) |
| `scenarios/score.py` | significant digits and the gate, as in Tier 3 |
| `scenarios/README.md` | the scenario, output, gate, and scoring formats, for the agents |
| `../admin/launch_codex_run.sh --tier2` | the launch: no `--search`, tmux session `tier2-<run>` |
| `../admin/verify_tier2.py` | the verification: builds, the 680 cases, both TS1 gates, the generator's tests, the day-night gate, the audit |

## Undisclosed inputs

Both fixed on 2026-10-22, before any Tier 2 run, in `~/NSF-NCAR/.tier2-heldout/`
(a repository with no remote, outside the run root); their hashes are in
`admin/runs.md`.
- **Five held-out day-night scenarios**, varying seeds, latitudes, seasons,
  start times, and initial NOx, with steps of 5 to 30 minutes and runs of up to
  48 hours, as `scenarios/README.md` announces. `make_heldout.sh` will be
  published in `heldout/` when the last run ends.
- **The TS1 ranking input:** 10,000 cells from `reference/make_cells_csv.py`
  with a 63-bit seed from Python's `secrets`.
