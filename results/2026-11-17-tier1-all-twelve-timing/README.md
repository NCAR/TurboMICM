# Timing of all twelve Tier 1 entries, 2026-11-17

All twelve finished Tier 1 entries were timed in one interleaved session
alongside MICM's default, MICM with 32-cell groups, and MICM's fastest
configuration, at 10 and 1 threads: the nine timed before or verified
earlier, and luna-high, terra-max, terra-xhigh, and terra-high, timed here for
the first time. The method, machine, input, reference, and builds are those of
[`2026-10-20-luna-xhigh-verification`](../2026-10-20-luna-xhigh-verification/README.md):
the clang builds of the entries and MICM's three programs in
`~/NSF-NCAR/tier1-work`, and the seed-2026 input and MICM final state, whose
SHA-256 sums (`43febf31…ad8b5`, `13329009…6ba7d`) were checked before the
session. The four new entries' builds are those of their verifications
(`2026-10-20-luna-high-verification`, `2026-11-17-terra-max-verification`,
`2026-11-17-terra-xhigh-verification`, `2026-11-17-terra-high-verification`),
with the same Apple clang 21.0.0 as the others.

| File | Contents |
|---|---|
| `timings.csv` | one row per timed run: elapsed_ms, solver statistics, process wall time, gate |
| `summary.md` | median and range of elapsed_ms per program and thread count |
| `environment.txt` | the machine and the exact commands timed |
| `load-audit.log` | one sample a minute: load, AC power, every process above 5 % of a core |
| `leaderboard.md` | `admin/leaderboard.py` at 10 threads over all seven timing sessions |
| `leaderboard-1-thread.md` | the same at 1 thread, over the five sessions that ran it |

## The wait

The session was first armed on 2026-11-17 at 08:49 and deferred by the SSE at
11:43, before it started: two macOS daemons, `suggestd` and `maild`, held the
machine near 1.1 cores, above the gate's one-core limit (`admin/runs.md`,
"terra-high"). After a reboot, the session was armed again at 20:03:58 with
the gate's default limit, one core (`IDLE_CPU_LIMIT=100`): the two daemons were
gone. `time-when-idle.sh` waited for 10 checks in a row, a minute apart, with no
Codex session record written for 10 minutes, total CPU below one core, and AC
power. macOS's own work after the reboot kept resetting it until 20:39:
XProtect's remediators (`XProtectRemediatorPirrit`, `XProtectRemediatorMRTv3`,
about one core each) from 20:07, the media analysis daemon (up to 1.9 cores)
and Spotlight's `spotlightknowledged` (one core, until about 20:33) after
them, and two one-minute spikes, at 20:36 and 20:39. The tenth quiet check
came at 20:49:01, 45 minutes after the session was armed
(`timing-all-twelve-wait.log`, kept in `~/NSF-NCAR/tier1-work`).

## Timing

`admin/time_entries.py` from 20:49 to 23:39: 5 rounds, each running all
fifteen programs at 10 threads and then at 1 thread, in rotating order, with a
5 s pause between runs, under `caffeinate`, each run's final state checked
against MICM's.

Every one of the 150 runs passed the gate. The largest relative differences
from MICM's default were those of each entry's verification: 2.96e-13
(astra-max), 2.78e-13 (luna-max), 2.51e-13 (terra-max), 2.03e-13 (astra-xhigh
and sol-high), 1.77e-13 (luna-xhigh and terra-xhigh), 1.26e-13 (sol-max),
1.25e-13 (terra-high), 5.34e-14 (luna-high), and 0 (astra-high, sol-xhigh, and
MICM's three programs).

| Program | 10 threads: median ms (range) | vs MICM | vs MICM's fastest | 1 thread: median ms (range) | vs MICM | vs MICM's fastest |
|---|---:|---:|---:|---:|---:|---:|
| MICM, one group per thread | 27,265 (26,777–27,933) | 1.00× |  | 146,347 (145,102–146,770) | 1.00× |  |
| MICM, 32-cell groups | 22,838 (22,817–23,514) | 1.19× |  | 136,013 (134,466–139,733) | 1.08× |  |
| MICM's fastest: vector8, in-place Mozart, 32-cell groups | 9,593 (9,513–9,783) | 2.84× | 1.00× | 50,981 (50,592–51,092) | 2.87× | 1.00× |
| sol-max | 7,114 (7,055–7,232) | 3.83× | 1.35× | 37,704 (37,662–37,820) | 3.88× | 1.35× |
| astra-xhigh | 7,370 (7,280–7,486) | 3.70× | 1.30× | 38,275 (38,175–38,457) | 3.82× | 1.33× |
| sol-xhigh | 8,193 (8,044–8,220) | 3.33× | 1.17× | 44,330 (44,189–44,441) | 3.30× | 1.15× |
| astra-max | 8,433 (8,336–8,608) | 3.23× | 1.14× | 45,514 (45,490–45,629) | 3.22× | 1.12× |
| sol-high | 8,568 (8,383–8,687) | 3.18× | 1.12× | 46,880 (46,730–46,983) | 3.12× | 1.09× |
| terra-max | 8,668 (8,591–8,707) | 3.15× | 1.11× | 43,691 (43,676–43,762) | 3.35× | 1.17× |
| astra-high | 9,193 (9,017–9,299) | 2.97× | 1.04× | 49,447 (49,305–49,687) | 2.96× | 1.03× |
| luna-xhigh | 22,057 (21,864–22,201) | 1.24× | 0.43× | 60,338 (60,313–60,344) | 2.43× | 0.84× |
| terra-high | 25,343 (25,173–26,021) | 1.08× | 0.38× | 144,516 (143,201–144,696) | 1.01× | 0.35× |
| terra-xhigh | 29,239 (28,957–30,459) | 0.93× | 0.33× | 168,884 (166,591–169,035) | 0.87× | 0.30× |
| luna-max | 31,541 (31,199–31,926) | 0.86× | 0.30× | 182,483 (182,291–182,738) | 0.80× | 0.28× |
| luna-high | 50,942 (50,542–51,727) | 0.54× | 0.19× | 311,252 (308,706–312,341) | 0.47× | 0.16× |

## Findings

- **terra-max is sixth of the twelve at ten threads and third on one.** At ten
  threads it takes 8,668 ms, 3.15× MICM's default and 1.11× MICM's fastest
  configuration, between sol-high (8,568 ms) and astra-high (9,193 ms); on one
  thread 43,691 ms, 3.35× and 1.17×, behind only sol-max and astra-xhigh. With
  its 64-cell groups it is 2.63× faster than MICM with 32-cell groups at ten
  threads and 3.11× on one; MICM with 64-cell groups was not timed here. It
  gains 5.04× from one thread to ten, less than the Sol and Astra entries,
  whose layout it shares (5.19–5.47×), which is why it places lower at ten
  threads than on one.
- **terra-high and terra-xhigh, one cell at a time, run near MICM's default.**
  terra-high takes 25,343 ms at ten threads (1.08× MICM's default) and 144,516
  ms on one (1.01×); terra-xhigh 29,239 and 168,884 ms (0.93× and 0.87×).
  terra-high is 1.15× faster than terra-xhigh at ten threads and 1.17× on one,
  although its LU makes 6,986 updates per factorization to terra-xhigh's 6,633,
  in MICM's species order, and its timer includes a million `calloc` and `free`
  pairs of about 366 KB (`src/rosenbrock.c:167-175,269-270` in its tree). Each
  takes 9,000,000 internal steps at ten threads, MICM's 9 per interval in every
  cell, where MICM's default takes 9,000 in ten groups. Where each spends its
  time is for the analysis pass to measure.
- **luna-high is the slowest of the twelve:** 50,942 ms at ten threads (0.54×
  MICM's default) and 311,252 ms on one (0.47×), with MICM's default grouping
  and solver work (9,000 internal steps at ten threads). luna-max is 1.62×
  faster at ten threads and 1.71× on one, luna-xhigh 2.31× and 5.16×.
- **Scaling.** From one thread to ten: terra-max 5.04×, terra-high 5.70×,
  terra-xhigh 5.78×, and luna-high 6.11×, the most of any program timed;
  luna-xhigh's 2.74× is still the least. MICM's three programs gain 5.31–5.96×.
- **The other eight keep their order** at both thread counts: sol-max,
  astra-xhigh, sol-xhigh, astra-max, sol-high, astra-high, luna-xhigh,
  luna-max. Among the twelve, terra-max falls between sol-high and astra-high
  at ten threads and between astra-xhigh and sol-xhigh on one; terra-high and
  terra-xhigh fall between luna-xhigh and luna-max at both, and luna-high last.
- **The measurements reproduce.** The eight entries' medians are within 2.9 %
  of the 2026-10-20 session at ten threads (luna-xhigh, +2.9 %) and within
  0.7 % on one thread; MICM's three programs within 3.3 % (its default, −3.3 %)
  and 1.3 %.
- **Self-reports** (single runs at ten threads): terra-max 9,632 ms against
  23,835 ms for MICM with 64-cell groups; terra-xhigh 30,496 against 29,617 ms
  for MICM's default; terra-high 25,689 against 26,324 ms; luna-high 48.5 s
  against 27.3 s. Here they take 8,668, 29,239, 25,343, and 50,942 ms, against
  27,265 ms for MICM's default.

Other load (`load-audit.log`, 170 samples): AC attached throughout. Processes
above 5 % of a core other than the benchmark and the timing script: macOS's
media analysis daemon in 10 samples (mean 0.5 core, at most 1.0), Spotlight's
knowledge updater in 5 (mean 0.4, at most 0.6), `contactsd` in 1 (1.15 cores,
at 21:27, with `AddressBookSourceSync` at 0.4), `opendirectoryd` in 1 (0.8, at
22:26), WindowServer in 13 (mean 0.07, at most 0.19), Terminal in 8 (0.07), and
Claude Code in 20 (mean 0.07, at most 0.10). Four of the heavy samples fell in
10-thread runs: 20:50 and 20:55 (round 1, media analysis at 0.5 and 0.7 core),
21:27 (round 2, `contactsd`), and 22:35 (round 4, media analysis at 0.9 core);
the rotation spreads them over the programs, and each median is of five rounds.
No Codex session ran. This Claude Code session waited, in a foreground wait,
through every 10-thread phase, and during the 1-thread runs read and wrote
files for the analysis of the Terra entries; its samples above 5 % fall at
phase boundaries and during the waits, at 5–10 % of a core.

These use seed 2026, which the agents knew. The ranking itself will use an
undisclosed seed.
