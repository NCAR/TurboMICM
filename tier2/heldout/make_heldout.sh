#!/bin/zsh
# Write Tier 2's undisclosed (held-out) day-night scenarios: cells.csv,
# photolysis.csv, scenario.txt.
#
#   make_heldout.sh MAKE_SCENARIO DIR [CELLS]
#
# MAKE_SCENARIO is tier2/scenarios/make_scenario.py. With the default 1,000 cells
# this makes the held-out set, and tier2/admin/make_scenario_data.sh then adds
# MICM's results at the gate's settings and the reference solutions.
#
# Fixed on 2026-10-22 by Claude Code (Opus 5.5), orchestrating Tier 2, before any
# Tier 2 run; kept outside the run root until the last Tier 2 run ends. Tier 3's
# held-out scenarios are public, so Tier 2 has its own. They vary what the
# agents' scenarios README announces: seeds, latitudes, seasons, start times, and
# initial NOx, steps of 5 to 30 minutes, and runs of up to 48 hours, recording a
# snapshot every 6 hours. arctic-summer-10min has the midnight sun (photolysis
# never reaches zero), south-48h-30min the southern summer, and night-start-5min
# starts at dusk, with photolysis off for its first steps.
set -eu
make=${1:A}
dir=$2
cells=${3:-1000}
mkdir -p $dir
cd $dir
python3 $make spring45-20min --cells $cells --seed 301 --dt 1200 --steps 72 \
    --diurnal --lat 45 --doy 110 --start 6 --snapshot-every 18
python3 $make tropics-polluted-15min --cells $cells --seed 302 --dt 900 --steps 96 \
    --diurnal --lat 5 --doy 200 --start 6 --nox-scale 3 --snapshot-every 24
python3 $make arctic-summer-10min --cells $cells --seed 303 --dt 600 --steps 144 \
    --diurnal --lat 70 --doy 172 --start 6 --nox-scale 0.5 --snapshot-every 36
python3 $make south-48h-30min --cells $cells --seed 304 --dt 1800 --steps 96 \
    --diurnal --lat -35 --doy 15 --start 6 --nox-scale 2 --snapshot-every 12
python3 $make night-start-5min --cells $cells --seed 305 --dt 300 --steps 288 \
    --diurnal --lat 40 --doy 300 --start 18 --nox-scale 0.1 --snapshot-every 72
