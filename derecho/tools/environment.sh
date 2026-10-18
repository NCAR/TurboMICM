#!/bin/bash
# Record the Derecho node a timing session ran on: OUT/environment-derecho.txt.
#
#   derecho/tools/environment.sh OUT
#
# admin/time_entries.py's environment.txt knows only macOS commands; this adds
# the node, its CPU and NUMA layout, the kernel, huge-page setting, compiler,
# modules, OpenMP binding, and the PBS job.
set -eu
out=$1
mkdir -p "$out"
{
  echo "date: $(date -Is)"
  echo "host: $(hostname)  job: ${PBS_JOBID:-none}  cpus allowed: $(grep Cpus_allowed_list /proc/self/status | cut -f2)"
  echo "kernel: $(uname -r)"
  lscpu | grep -E 'Model name|^CPU\(s\)|Thread\(s\) per core|Core\(s\) per socket|Socket|NUMA node\(s\)|L1d|L1i|L2|L3|max MHz|min MHz|boost' | sed 's/  */ /g'
  echo "transparent hugepages: $(cat /sys/kernel/mm/transparent_hugepage/enabled)"
  echo "compiler: $(which gcc) $(gcc --version | head -1)"
  echo "modules: $(module -t list 2>&1 | tr '\n' ' ')"
  echo "OMP_PROC_BIND=${OMP_PROC_BIND:-}  OMP_PLACES=${OMP_PLACES:-}  OMP_NUM_THREADS=${OMP_NUM_THREADS:-}"
  echo "load: $(cat /proc/loadavg)"
} > "$out/environment-derecho.txt"
cat "$out/environment-derecho.txt"
