# Tier 2 exercises, measured, 2026-10-23

The prompt told the agents that reviewers would "add a reaction to TS1 and
regenerate, change a rate parameter, and run the generator on another
mechanism". The blind reviewers do these three exercises on paper
([`2026-10-23-tier2-blind-review`](../2026-10-23-tier2-blind-review/README.md));
here Claude Code (Opus 5.5) did them on the machine, for all nine entries, on
the Apple M5 Pro.

| Folder | Contents |
|---|---|
| `ab/<run>/` | exercises (a) and (b): `exercises.md` and `exercises.json`, from `admin/tier2_exercises.py` |
| `c/<run>/` | exercise (c): each agent's `exercise-c.md` (commands, findings) and `exercise-c.json`, its converter `convert.py`, and any helper scripts |
| `c/canonical/<mechanism>/` | the shared inputs: `canonical.json`, `cells.csv` (4 identical cells), `cells64.csv` (64 distinct cells), MICM's header for the mechanism, and MICM's final states |
| `c/mechanisms/` | MUSICA's v0 Chapman and Carbon Bond 5, as the blind reviewers received them |

## Exercises (a) and (b): change TS1 and regenerate

Every entry builds its generic and its specialized TS1 benchmarks from one TS1
configuration file, so these two exercises were the same for all nine.
`admin/tier2_exercises.py` ran them in a fresh worktree of each entry's Tier 2
branch, with clang and OPTFLAGS "-O3 -DNDEBUG":
- **baseline:** `make env` and `make cpu bench`; both benchmarks on
  `reference/ts1_cells_100_seed1.csv` (1 thread, 100 steps of 30 s);
- **(a) add a reaction:** a copy of TS1's NO + O3 → NO2 + O2 Arrhenius
  reaction appended to the configuration, a second channel at the same rate;
  then `make bench` and both benchmarks;
- **(b) change a rate parameter:** the configuration restored, then that
  reaction's A multiplied by 1.5; `make bench` and both benchmarks.

An exercise passes when `make bench` succeeds and rebuilds the specialized
benchmark, both benchmarks' results change from the baseline, and the
specialized result agrees with the generic one within 1e-10
(`reference/compare_final.py`). The script edits each entry's configuration
in its own format: CAMP-style records for five entries, sol-max's variant of
them (`reactant sequence`, `product sequence`), and the Luna entries' own
top-level `reactions` lists.

| Entry | TS1 configuration | Clean build (s) | (a) rebuild (s) | (a) specialized vs generic | (b) rebuild (s) | (b) specialized vs generic | Verdict |
|---|---|---:|---:|---:|---:|---:|---|
| astra-max | `configs/ts1.json` | 19.2 | 14.9 | 0 | 14.8 | 0 | pass, pass |
| astra-xhigh | `configs/ts1.json` | 75.0 | 68.6 | 4.4e-14 | 67.6 | 4.7e-14 | pass, pass |
| astra-high | `mechanisms/ts1.json` | 28.0 | 23.1 | 0 | 23.4 | 0 | pass, pass |
| sol-max | `configs/ts1.json` | 35.6 | 30.2 | 0 | 29.0 | 0 | pass, pass |
| sol-xhigh | `mechanisms/ts1.json` | 20.8 | 17.8 | 0 | 17.6 | 0 | pass, pass |
| sol-high | `configs/ts1.json` | 19.3 | 16.1 | 4.4e-14 | 15.8 | 4.7e-14 | pass, pass |
| luna-max | `mechanisms/ts1.json` | 7.0 | 2.6 | 0 | 2.7 | 0 | pass, pass |
| luna-xhigh | `config/ts1.json` | 71.1 | 67.2 | 3.1e-14 | 62.3 | 3.4e-14 | pass, pass |
| luna-high | `mechanisms/ts1.json` | 15.1 | 13.0 | 0 | 13.0 | 0 | pass, pass |

- **All nine pass both.** Each build noticed the edited configuration,
  regenerated, and rebuilt; the generic and specialized results moved
  together and stayed within 4.7e-14 of each other.
- **The cost of a mechanism change is a rebuild of 3 to 69 s.** It is
  dominated by compiling the generated C; luna-max's is the smallest.
  The entries ran three or four at a time, so these times run high next to
  the verification's single builds and compare only roughly.
- **The script's first pass** knew only the CAMP-style format and stopped on
  sol-max and the three Luna entries (01:13 to 01:20). Taught their formats,
  it reran those four from fresh worktrees (01:23 to 01:28). The entries were
  never at fault.

## Exercise (c): another mechanism

The reviewers' mechanisms were Chapman (5 species, 7 reactions) and Carbon Bond
5 (67 species, 200 reactions, 14 of them emissions without reactants), from
MUSICA's `configs/v0`. A third variant, Carbon Bond 5 without its emissions,
was there in case reactions without reactants were the obstacle; none of the
nine entries needed it.

**The shared problem.** So that every entry solved the same problem,
`admin/v0_to_canonical.py` settled the configurations' meaning once, as MICM's
former v0 parser did (MICM `6755647d^`, `solver_config.hpp`):
- rate constants converted to mol m-3, and Ea to C;
- photolysis and emissions made user-defined rates, a label times a scaling
  factor;
- M a third body parameterized on the air density, as in TS1 (Carbon Bond 5
  declares M without the tracer type that would say so).

Labels became C identifiers (`PHOTO_*`, `EMIS_*`). Five Carbon Bond 5
photolysis rates have no initial value and are zero, as in MICM. MICM's
reference is `reference/micm_ts1_ref.cpp`, built unchanged against a header
that `admin/micm_header.py` writes for each mechanism. Its default settings
(relative tolerance 1e-6, absolute 1e-3) run one warm-up step and then 60 steps
of 60 s, on two inputs:
- `cells.csv`, 4 identical cells;
- `cells64.csv`, 64 distinct cells from `admin/vary_cells.py`, with
  temperature ±10 %, pressure ±20 %, lognormal concentrations, and a sun factor
  on the photolysis rates.

MICM gives the same final state on the 64 cells with its default grouping and
with 32-cell groups. All 540 steps per run are accepted, so the step sequence
does not depend on the grouping.

**The method.** One Claude Code agent (Opus 5.5) per entry worked in its own
folder with a fresh worktree of the entry's Tier 2 branch, under these rules:
- no change to the entry;
- the agent writes its own converter from the canonical form into the entry's
  input format, modeled on the entry's TS1 file;
- the agent puts the converted file where the build reads TS1's configuration,
  or passes it through a Makefile variable;
- the agent runs `make bench` with clang and OPTFLAGS "-O3 -DNDEBUG", then
  both benchmark programs (`--threads 1 4 60 60 cells.csv`, and 64 cells for
  the second input);
- every final state is compared with the other program's and with MICM's
  (`compare_final.py`, 1e-10).

At the end, each agent restored the TS1 configuration. Each entry then rebuilt
TS1, which passed again, and `git status --porcelain` showed no change. Claude
Code reran `compare_final.py` on every final state the agents wrote. Every
protocol run passes; the only mismatches were the agents' own extra runs on
other inputs: astra-xhigh's longer timing runs, astra-max's 64 identical cells,
and a sol-xhigh pass with stale binaries, which its agent caught and redid.

| Entry | How the mechanism reached the build | Converter (lines) | Generated C: Chapman, Carbon Bond 5 (lines) | `make bench`, Carbon Bond 5 (s) | 4 cells: specialized vs MICM | Generated kernels at 4 cells | 64 cells: specialized vs generic | 64 cells: specialized vs MICM |
|---|---|---:|---:|---:|---:|---|---:|---:|
| astra-max | swap of `configs/ts1.json` | 124 | 863; 36,844 | 5.2 | 5.7e-15 | yes, the runtime-stride form | 0 | 2.7e-14 |
| astra-xhigh | swap of `configs/ts1.json` | 119 | 825; 25,494 | 7.5 | 5.7e-15 | yes, a padded 4-cell tile | 2.6e-14 | 2.7e-14 |
| astra-high | swap of `mechanisms/ts1.json` | 76 | 569; 15,649 | 7.0 | 0 | no: generic kernels below 32 cells | 0 | 0 |
| sol-max | swap of `configs/ts1.json` | 90 | 902; 22,867 | 6.9 | 4.4e-15 | yes | 0 | 2.9e-14 |
| sol-xhigh | swap of `mechanisms/ts1.json` | 118 | 670; 36,218 | 4.7 | 0 | yes, the runtime-stride form | 0 | 0 |
| sol-high | swap of `configs/ts1.json` | 84 | 452; 9,685 | 6.2 | 4.6e-15 | no: generic kernels below 32 cells | 2.6e-14 | 2.7e-14 |
| luna-max | `TS1_JSON=` on the make command line | 124 | 304; 9,682 | 1.4 | 4.5e-15 | yes | 0 | 3.0e-14 |
| luna-xhigh | `CONFIG=` on the make command line | 109 | 602; 16,456 | 7.1 | 5.0e-15 | yes | 2.7e-14 | 3.4e-14 |
| luna-high | swap of `mechanisms/ts1.json` | 102 | 301; 10,610 | 5.2 | 0 | yes | 0 | 2.3e-14 |

Each difference is the largest over the three mechanisms. Every generic
program agrees with MICM as closely as its specialized one. The generator ran
in 0.03 to 0.13 s.

- **All nine generators took both mechanisms with no change to the entry and
  no driver:** only a converter of 76 to 124 lines into the entry's own
  format, and every result agrees with MICM within 3.4e-14.
  - The 14 emissions, reactions without reactants, which TS1 lacks, worked in
    every entry.
  - astra-high and sol-xhigh reproduce MICM bit for bit on all three
    mechanisms, and luna-high does at 4 cells.
- **The generated kernels ran.** In astra-high and sol-high they serve only
  full 32-cell groups, and narrower groups take the generic kernels. In those
  two entries the 4-cell runs therefore tested only the generic path; sol-high's
  agent found this. The 64-cell input fills two groups.
  - Every agent showed from the code path, and most also from the binaries or
    from call counters, that all nine specialized programs ran their generated
    kernels on 64 cells.
  - sol-high and astra-high also pass at 4 cells when rebuilt for that width,
    with the Makefile variables `SPECIALIZE_WIDTH=4` and `GROUP_SIZE=4`.
- **Every TS1 benchmark was already mechanism-generic.** Each entry's
  benchmark programs find the cells file's columns by name and take the
  mechanism from the generated code, so each ran a new mechanism once the
  build was pointed at it. Two Makefiles take the configuration as a
  variable; the other seven name `ts1.json` directly.
- **Traps a user would meet,** found by the agents:
  - **Units:** the astra-xhigh and astra-high readers default to molecule cm-3
    and silently rescale A, k0_A, and kinf_A unless each reaction declares mol
    m-3. luna-xhigh's CAMP reader applies molecule-cm-3 factors to Troe
    parameters, so its agent wrote the entry's own descriptor format instead.
  - **Configuration checks:** astra-xhigh fails on an initial rate parameter
    that no reaction uses; sol-max rejects a species that no reaction uses.
  - **Undocumented options:** luna-max's `TS1_JSON` override is undocumented;
    its README documents `make mechanism CONFIG=…`, which builds only a
    generic library.
  - **Photolysis naming:** the scenario programs of sol-max and sol-xhigh treat
    only parameters whose names begin with `j` as photolysis, TS1's naming.
    This was outside the exercise.
  - **A mechanism without user-defined rates** makes luna-high's `src/main.c`
    divide by zero: undefined behavior, with four `-Wpedantic` warnings. It
    was tried as a side check.
  - **Its own test on Carbon Bond 5:** astra-max's generator test, pointed at
    Carbon Bond 5, makes both of its solvers exceed 1,000 steps on its
    synthetic inputs.
- **macOS's make (3.81) compares times to the second.** A configuration
  swapped in within a second of the previous build can leave stale binaries.
  Three agents hit this in extra passes and reran them after pausing before
  each swap; every reported build rebuilt fully. In exercises (a) and (b),
  the benchmark runs separated each build from the next edit, and the script
  checked that every build rebuilt the specialized program.
- **Setup, a deviation in the instructions:** they said `make env` uses an
  existing environment. In fact it creates one inside the worktree from
  conda-forge.
  - astra-max's, sol-max's, and sol-high's agents ran it as prescribed; it
    downloads any package not yet cached (about 25 MB for sol-high).
  - The others reused an environment that the same entry's `make env` had
    built for verification, from the same `environment.yml`, or used the
    Makefile's fallback to the system `python3`.
- **Times:** the build and generator times above come from nine agents
  working at once. They show scale, not a ranking.
