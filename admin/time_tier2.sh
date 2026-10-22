#!/bin/zsh
# Time Tier 2 entries on the ranking input, once the machine is idle.
#
#   admin/time_tier2.sh OUT THREADS ROUNDS NAME=COMMAND...
#
# For example:
#   admin/time_tier2.sh ~/NSF-NCAR/tier2-work/timing/astra-max 18 5 \
#       "micm=~/NSF-NCAR/tier2-work/bin/micm_ts1_ref" \
#       "micm-g32=~/NSF-NCAR/tier2-work/bin/micm_ts1_ref --group-size 32" \
#       "astra-max/specialized=.../bin/ts1_bench_specialized" ...
#
# 1. Waits until the machine has been idle for 10 checks in a row, a minute apart:
#    no Codex session record written in the last 10 minutes, the processes' CPU
#    summing to less than one core, and AC power (log: OUT/wait.log). NOWAIT=1
#    skips the wait.
# 2. Runs admin/time_entries.py under caffeinate on the Tier 2 ranking input, the
#    10,000 cells of the undisclosed seed, 100 steps of 30 s, with MICM's final
#    state for it as the gate (TIER2_INPUT and TIER2_REFERENCE replace them),
#    THREADS (comma-separated) and ROUNDS as given, into OUT.
# 3. Records the load once a minute while it times (OUT/load-audit.log): load
#    averages, the power source, and the five busiest processes.
# 4. Touches OUT/DONE, and writes the exit status to OUT/exit.
set -u
here=${0:A:h}
out=${1:A}
threads=$2
rounds=$3
shift 3
input=${TIER2_INPUT:-$HOME/NSF-NCAR/.tier2-heldout/ts1/cells.csv}
reference=${TIER2_REFERENCE:-$HOME/NSF-NCAR/.tier2-heldout/ts1/micm_final.csv}
mkdir -p $out
rm -f $out/DONE $out/exit

idle_once() {
  local recent cpu power
  recent=$(find $HOME/.codex/sessions -name 'rollout-*.jsonl' -mmin -10 2>/dev/null | wc -l | tr -d ' ')
  cpu=$(ps -axo pcpu= | awk '{s += $1} END {printf "%d", s}')
  power=$(pmset -g batt | head -1 | command grep -c "AC Power")
  echo "$(date '+%F %T') records written in 10 min: $recent; CPU: $cpu %; AC: $power" >> $out/wait.log
  (( recent == 0 && cpu < 100 && power == 1 ))
}

if [[ ${NOWAIT:-0} != 1 ]]; then
  streak=0
  while (( streak < 10 )); do
    if idle_once; then (( streak++ )); else streak=0; fi
    (( streak < 10 )) && sleep 60
  done
  echo "$(date '+%F %T') idle for 10 checks: timing" >> $out/wait.log
fi

(
  while [[ ! -e $out/DONE ]]; do
    {
      echo "== $(date '+%F %T') $(sysctl -n vm.loadavg) $(pmset -g batt | head -1)"
      ps -axo pcpu=,comm= -r | head -5
    } >> $out/load-audit.log
    sleep 60
  done
) &
audit=$!

caffeinate -ims python3 $here/time_entries.py --input $input --reference $reference --out $out \
    --threads $threads --rounds $rounds "$@" > $out/time_entries.log 2>&1
rc=$?  # not "status", which zsh reserves
echo $rc > $out/exit
echo "exit $rc" >> $out/time_entries.log
touch $out/DONE
wait $audit 2>/dev/null
exit $rc
