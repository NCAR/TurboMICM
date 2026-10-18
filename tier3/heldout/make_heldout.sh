#!/bin/zsh
# Write Tier 3's undisclosed (held-out) scenarios: cells.csv, photolysis.csv, scenario.txt.
#
#   make_heldout.sh MAKE_SCENARIO DIR [CELLS]
#
# MAKE_SCENARIO is tier3/scenarios/make_scenario.py. With the default 1,000 cells
# this makes the held-out set, and tier3/admin/make_references.sh then adds the
# reference solutions. With CELLS = 10000 it makes the 10,000-cell versions used
# for timing, whose first 1,000 cells are the held-out ones.
#
# The SSE approved these parameters on 2026-10-18, before any Tier 3 result. They
# vary what the disclosed README announces: seeds, latitudes, seasons, and 48-hour
# runs. Steps stay within 5 to 30 minutes and initial NOx within x0.1 to x10, as in
# the disclosed set; every scenario starts at 06:00 local solar time and records a
# snapshot every 6 hours.
set -eu
make=${1:A}
dir=$2
cells=${3:-1000}
mkdir -p $dir
cd $dir
python3 $make summer-15min --cells $cells --seed 201 --dt 900 --steps 96 \
    --diurnal --lat 40 --doy 172 --start 6 --snapshot-every 24
python3 $make tropics-20min --cells $cells --seed 202 --dt 1200 --steps 72 \
    --diurnal --lat 10 --doy 80 --start 6 --snapshot-every 18
python3 $make winter60-10min --cells $cells --seed 203 --dt 600 --steps 144 \
    --diurnal --lat 60 --doy 355 --start 6 --snapshot-every 36
python3 $make polluted-48h-30min --cells $cells --seed 204 --dt 1800 --steps 96 \
    --diurnal --lat 30 --doy 200 --start 6 --nox-scale 5 --snapshot-every 12
python3 $make clean-48h-5min --cells $cells --seed 205 --dt 300 --steps 576 \
    --diurnal --lat 50 --doy 260 --start 6 --nox-scale 0.2 --snapshot-every 72
