# Tier 2 ranking, 2026-10-23

Tier 2's primary ranking, by the rules on Tier-2-Design. Each of the nine
entries' specialized TS1 benchmark, built with clang, was timed on the Apple M5
Pro beside its own generic solver, its Tier 1 commit, and three MICM programs.
All ran in one session on the undisclosed ranking input: 10,000 cells from
`reference/make_cells_csv.py` with the seed in
[`tier2/heldout/SEED`](../../tier2/heldout/SEED), 100 steps of 30 s. Each time
is the median of 5 runs, on an idle machine.

| Folder | Session |
|---|---|
| `18-threads/` | 18 threads, all cores, from 00:34:39 to 01:11:50: the ranking |
| `1-thread/` | 1 thread: to follow, after the blind review and the exercises |

Each folder holds:
- `timings.csv`: every run;
- `summary.md`: each program's median and range;
- `ranking.md`: the table below, from `admin/rank_tier2.py`;
- `programs.txt`: the programs, as given to `admin/time_tier2.sh`;
- `environment.txt`: the machine and the input;
- `load-audit.log`: the busiest processes, once a minute;
- `wait.log`: the idle checks before the session.

## The session at 18 threads

First, `admin/time_tier2.sh` waited for 10 idle checks a minute apart, from
00:25:38 to 00:34:39. Each check required:
- no Codex session record written in 10 minutes;
- CPU below one core;
- AC power.

Then it ran `admin/time_entries.py` under `caffeinate`: 30 programs in 5
rotating rounds, 150 runs. Each run's final state was checked against MICM's
for the input with `compare_final.py` (1e-10), and all 150 pass.

The programs:
- **MICM** e15898ce's TS1 benchmark, `micm_ts1_ref`, run three ways:
  - with its default grouping, one group of cells per thread;
  - with 32-cell groups;
  - in MICM's fastest Tier 1 configuration: vector-ordered matrices of width 8,
    the in-place Mozart decomposition, and 32-cell groups.
- **Each entry,** three programs:
  - `bin/ts1_bench_specialized`, from the clang build in its verification
    (`results/2026-10-2x-tier2-<run>-verification`);
  - `bin/ts1_bench_cpu`, its generic solver, from the same build;
  - `bin/ts1_bench_cpu` from its Tier 1 commit, built here with its OpenMP flags
    on the make command line.

Medians, ms (MICM: default 12,192.1; 32-cell groups 11,492.9; fastest
configuration 4,705.3):

| Rank | Entry | Specialized (ms) | Spread | Generic, Tier 2 (ms) | Tier 1 build (ms) | Gain from specialization | Over MICM's default | Over MICM, 32-cell groups | Over MICM's fastest |
|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | astra-xhigh | 2,048.9 | 8.0 % | 2,991.7 | 2,972.7 | 1.46 | 5.95 | 5.61 | 2.30 |
| 2 | sol-max | 2,403.4 | 2.2 % | 2,979.7 | 2,993.6 | 1.24 | 5.07 | 4.78 | 1.96 |
| 3 | sol-xhigh | 2,706.9 | 5.0 % | 3,544.9 | 3,581.9 | 1.31 | 4.50 | 4.25 | 1.74 |
| 4 | sol-high | 2,796.9 | 4.5 % | 3,780.5 | 3,838.7 | 1.35 | 4.36 | 4.11 | 1.68 |
| 5 | astra-high | 3,156.8 | 1.7 % | 4,028.9 | 4,047.7 | 1.28 | 3.86 | 3.64 | 1.49 |
| 6 | astra-max | 3,216.4 | 8.8 % | 3,929.2 | 3,906.6 | 1.22 | 3.79 | 3.57 | 1.46 |
| 7 | luna-max | 5,381.0 | 7.9 % | 15,110.5 | 15,131.5 | 2.81 | 2.27 | 2.14 | 0.87 |
| 8 | luna-high | 7,567.9 | 2.8 % | 22,954.3 | 22,795.6 | 3.03 | 1.61 | 1.52 | 0.62 |
| 9 | luna-xhigh | 10,059.6 | 1.2 % | 10,756.3 | 10,760.6 | 1.07 | 1.21 | 1.14 | 0.47 |

The spread is the range of the five runs over their median. The gain from
specialization is the generic median over the specialized one. The last three
columns divide MICM's medians by the entry's.

## How firm the places are

- **Every adjacent pair but one is consistent across all five rounds:** the
  higher-placed entry's specialized benchmark was the faster in each round.
  - sol-xhigh and sol-high are 3.3 % apart, and their ranges overlap, but
    sol-xhigh was faster in all five rounds.
  - Every other gap is 12.6 % or more.
- **Places 5 and 6 are a tie.** astra-high was faster than astra-max in 3 of the
  5 rounds, and their medians differ by 1.9 %, within astra-max's spread.
  astra-max's own provisional sessions measured it at about 2,930 ms, faster
  than astra-high.
- **Ranking by the ratio instead gives the same order,** but for places 5 and 6
  swapped. The ratio is the median, over the rounds, of MICM's default over the
  entry; it cancels a round that ran slow for both.

## Spread

- **The entries' 27 programs varied a median 3.7 % across the rounds.** Their
  specialized benchmarks ranged from 1.2 % (luna-xhigh) to 8.8 % (astra-max).
- **MICM's programs varied more, 12 to 20 %.** The six programs at the head of
  the rotation (MICM's three, then astra-max's three) ran 8.5 to 15 % slower in
  rounds 2 to 4 than in rounds 1 and 5.
  - In every round those six ran back to back, in a block of about 80 s.
  - The other 24 programs show no such pattern. Their times averaged 1.000 of
    their medians in the first ten places of a round and 1.001 elsewhere.
  - The load audit shows no process above 11 % of a core during those blocks.
    The cause is not known.
- **The speedups over MICM therefore carry about 10 %;** the order of the
  entries does not depend on them.

## Load audit

Once a minute, the audit recorded the timed program and the machine's other
busy processes, always on AC power.
- **Usually under 5 % of a core each:**
  - this Claude Code session;
  - the managed laptop's agents (BESAgent, Jamf, a printer installer, Druva
    inSync);
  - WindowServer, Terminal, trustd.
- **Briefly more:**
  - Spotlight's indexer: 26.5 % at 00:36:39 and 15.5 % at 01:04:42;
  - logd: 14.1 % at 00:56:41;
  - the harness's own gate check, single-threaded Python at up to 99.5 %,
    between runs.

During the session, Claude Code wrote `admin/rank_tier2.py`, tried it on an
earlier session's file, edited CLAUDE.md, and otherwise waited. The exercises
and the blind reviews started after the session ended.

## Findings

- **astra-xhigh is first,** 5.95 times faster than MICM's default and 2.30
  times faster than MICM's fastest configuration. Its gain from
  specialization, 1.46, is the largest of the six entries whose generic
  solvers were already fast.
  - It and sol-max were level before specialization: their Tier 1 builds took
    2,972.7 and 2,993.6 ms here.
- **Specialization reordered only the bottom of the field.**
  - Before it, the Tier 1 builds ran in the order astra-xhigh, sol-max,
    sol-xhigh, sol-high, astra-max, astra-high, luna-xhigh, luna-max,
    luna-high.
  - luna-max and luna-high gained 2.81 and 3.03 times and passed luna-xhigh,
    which gained 1.07.
  - astra-high and astra-max swapped, within noise.
- **Six entries beat MICM's fastest configuration and three don't.** The six
  are the Sol and Astra entries, by 1.46 to 2.30 times; all nine beat MICM's
  default, by 1.21 to 5.95.
- **Every generic solver kept its Tier 1 speed.** Each entry's generic path in
  its Tier 2 commit is within 1.6 % of its Tier 1 build, run in the same
  session.
