#!/bin/bash
# MICM's leading configurations at 256 threads (two per core), 3 rounds: the
# vector 8, 4, and 2 orderings with in-place Mozart and Doolittle LU, vector
# 128 and vector 8 with separate factors, plus the default, each with one cell
# group per thread and with 32-cell groups.
#   derecho/tools/submit.sh micm-configs-t256 128 01:00:00 pbs/micm_t256.sh
set -eu
here=$(cd "$(dirname "$0")/.." && pwd)
cd "$here"
export THREADS=256 ROUNDS=3 LABEL=t256
export ONLY='^vector(8|4|2)-in-place-(mozart|doolittle)$|^vector128-in-place-mozart$|^vector8-separate-mozart$'
tools/micm_configs.sh work/micm-configs work/cells_2026.csv work/micm-reference/micm_final.csv
