#!/bin/bash
# Submit a command to Derecho's develop queue (cpudev) as a PBS job.
#
#   [PLACE=excl|shared] [MEM=...] derecho/tools/submit.sh NAME NCPUS WALLTIME COMMAND...
#
# The job runs COMMAND with bash from the repository's derecho/ directory after
# sourcing env.sh, on one node with NCPUS cores for WALLTIME (HH:MM:SS). The
# dev nodes advertise their 128 cores as ncpus (the second hardware thread of
# each core is not scheduled separately), so NCPUS=128 takes a whole node.
# Output streams to derecho/logs/NAME.log (PBS -k oed), starting with the node,
# its CPU set, and the compiler; the job id is printed. The account is
# $TURBOMICM_PBS_ACCOUNT (env.sh).
set -eu
name=$1 ncpus=$2 walltime=$3; shift 3
here=$(cd "$(dirname "$0")/.." && pwd)
source "$here/env.sh"
mkdir -p "$here/logs"
command=$(printf '%q ' "$@")
place=${PLACE:+#PBS -l place=$PLACE}
mem=${MEM:+:mem=$MEM}
qsub <<JOB
#PBS -N $name
#PBS -A $TURBOMICM_PBS_ACCOUNT
#PBS -q develop
#PBS -l select=1:ncpus=$ncpus:ompthreads=$ncpus$mem
#PBS -l walltime=$walltime
#PBS -j oe
#PBS -k oed
#PBS -o $here/logs/$name.log
$place
cd "$here"
source env.sh 2>/dev/null
echo "host: \$(hostname)  start: \$(date -Is)  ncpus: $ncpus  jobid: \$PBS_JOBID"
echo "cpus allowed: \$(grep Cpus_allowed_list /proc/self/status | cut -f2)  nproc: \$(nproc)"
echo "compiler: \$(which gcc) \$(gcc --version | head -1)"
echo "modules: \$(module -t list 2>&1 | tr '\n' ' ')"
echo "command: $command"
$command
echo "end: \$(date -Is)  exit: \$?"
JOB
