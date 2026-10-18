#!/bin/bash
# Verify the Derecho build on a held-out input: a seed the tuning never used.
#
#   derecho/pbs/verify_heldout.sh OUT SEED
#
# Writes the 10,000-cell input for SEED, MICM's final state for it (the
# default configuration at 128 threads: MICM's result does not depend on the
# thread count), then runs derecho/tools/verify_port.sh against it: the ported
# test suite and the gate at the default, one thread, groups of 16 and 64, and
# the generic kernels. Run with a whole node:
#   derecho/tools/submit.sh verify-heldout 128 01:00:00 pbs/verify_heldout.sh OUT 7
set -eu
here=$(cd "$(dirname "$0")/.." && pwd)
repo=$(cd "$here/.." && pwd)
out=$(realpath -m "$1"); seed=$2
cd "$here"
mkdir -p "$out"
cells=$here/work/cells_$seed.csv
[[ -f $cells ]] || python3 "$repo/reference/make_cells_csv.py" 10000 "$seed" --output "$cells"
echo "== MICM's final state for seed $seed ($(date -Is))"
THREADS=128 tools/micm_reference.sh "work/micm-reference-$seed" "$cells"
cp "work/micm-reference-$seed/micm_run.txt" "$out/micm_run_seed_$seed.txt"
echo "== the port against it ($(date -Is))"
SETTINGS=$'default|\n1 thread|--threads 1\ngroups of 16|--group-size 16\ngroups of 64|--group-size 64\ngeneric kernels|--generic' \
    tools/verify_port.sh "$out" "work/micm-reference-$seed/micm_final.csv" "$cells"
echo "== done ($(date -Is))"
