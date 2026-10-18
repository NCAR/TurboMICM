#!/bin/zsh
# Run admin/time_entries.py --scenario on every scenario of the given sets, one
# after another, with the same programs, writing OUT/<set>/<scenario>/.
#
#   time_tier3.sh OUT PROGRAMS_FILE SET[,SET...] [time_entries.py options ...]
#
# PROGRAMS_FILE: one NAME=COMMAND per line ('#' lines and blank lines skipped; $HOME is
# expanded), for example admin/tier3_micm_scan.txt.
# SETs: disclosed (1,000 cells), heldout (1,000 cells), disclosed-10000,
# heldout-10000. The 10,000-cell sets are scored against the 1,000-cell
# references of the same scenarios (their first 1,000 cells are the same).
# Between runs only: never while a contestant runs.
set -eu
TM=${0:A:h:h}
T3=$HOME/NSF-NCAR/ODE/benchmarks/tier3
H=$HOME/NSF-NCAR/.tier3-heldout
out=$1 programs_file=$2 sets=$3
shift 3
programs=()
while IFS= read -r line; do
  [[ -z $line || $line == \#* ]] && continue
  programs+=("${line//\$HOME/$HOME}")
done < $programs_file
(( ${#programs} )) || { echo "no programs in $programs_file" >&2; exit 1; }
for set in ${(s:,:)sets}; do
  case $set in
    disclosed) dirs=($T3/scenarios/*(/)) refs=$T3/scenarios ;;
    heldout) dirs=($H/heldout/*(/)) refs=$H/heldout ;;
    disclosed-10000) dirs=($H/disclosed-10000/*(/)) refs=$T3/scenarios ;;
    heldout-10000) dirs=($H/heldout-10000/*(/)) refs=$H/heldout ;;
    *) echo "unknown set $set" >&2; exit 1 ;;
  esac
  for d in $dirs; do
    name=${d:t}
    ref=$refs/$name/reference_snapshots.csv
    [[ -f $ref ]] || ref=$refs/$name/reference_final.csv
    [[ -f $ref ]] || { echo "no reference for $name in $refs" >&2; exit 1; }
    echo "== $set/$name: $(date '+%H:%M:%S')"
    python3 $TM/admin/time_entries.py --scenario $d/scenario.txt --reference $ref \
        --out $out/$set/$name "$@" $programs
  done
done
echo "== done: $(date '+%H:%M:%S')"
