# MICM's 24 CPU configurations, timed with the entries, 2026-10-15

MICM's TS1 reference (`reference/micm_ts1_ref.cpp`) built in all 24 of
`micm_bench`'s CPU configurations, and timed alongside the Tier 1 entries:
6 matrix layouts (standard; vector ordering with 1, 2, 4, 8, or 128 cells
interleaved) × LU storage (in-place or separate) × LU algorithm (Mozart or
Doolittle), each with one cell group per thread and with 32-cell groups.

| File | Contents |
|---|---|
| `micm_ts1_ref_configs.cpp` | the configurable copy of the reference, made by `admin/micm_configs_source.py` |
| `check.csv` | every configuration on the 100-cell example: Mozart configurations reproduce MICM's default exactly, Doolittle ones within 3.6e-14 |
| `inplace-mozart/` | session 1 (14:38–14:58): in-place Mozart in all six layouts, both groupings, and the four entries; 3 rounds |
| `others/` | session 2 (15:02–16:05): the other 18 configurations, with the default and vector8 in-place Mozart again as baselines; 3 rounds |
| `leaderboard.md`, `matrices.md` | `admin/leaderboard.py` over both sessions: every program ranked, and the configurations as layout-by-LU matrices |

Each session is `admin/time_entries.py` output (`timings.csv`, `summary.md`,
`environment.txt`), run by `admin/micm_configs.sh` with `THREADS=10 ROUNDS=3`:
10,000 cells (seed 2026), 100 steps of 30 s, clang `-O3 -DNDEBUG`, on AC power.
Every one of the 168 runs passed the gate against MICM's default result.

Findings, at 10 threads:
- **Fastest:** vector8 ordering with 32-cell groups and in-place Mozart LU,
  9,752 ms, 2.85 times faster than the default (27,751 ms).
- **Layout dominates:** vector8 gains 2.2–2.4 times with 32-cell groups and
  about 1.5 with one group per thread.
- **Storage and algorithm:** separate LU storage costs 2–14 % more (28–36 % with
  vector128 and 32-cell groups). Mozart and Doolittle differ by a median of
  0.9 %.
- **The entries:** 1.05 to 1.36 times faster than MICM's fastest configuration,
  and 3.00 to 3.88 times faster than its default.

Other load:
- **Web-research agents** ran from about 14:45 to 15:35, covering session 1's
  rounds 2 and 3 and the start of session 2. With them, Claude Code used about
  16 % of one core while the benchmark used 9.6 of the 10.
- **Checking of their sources** (network requests and small scripts) continued
  through session 2.
- **Two runs in session 2** were slow outliers: vector4 with separate storage
  and one group per thread, at 25.3 s and 27.8 s against a usual 20 s. The
  medians of 3 are unaffected.

The rotating rounds spread any effect across all the programs.
