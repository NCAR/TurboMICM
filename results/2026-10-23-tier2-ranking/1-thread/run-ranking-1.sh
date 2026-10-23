#!/bin/zsh
# Tier 2's ranking session at one thread: MICM's three programs and each entry's
# specialized and generic benchmarks (the Tier 1 builds ran at 18 threads only,
# within 1.6 % of the generic ones), 5 rotating rounds, after the idle wait.
W=$HOME/NSF-NCAR/tier2-work
args=("micm=$W/bin/micm_ts1_ref" "micm-g32=$W/bin/micm_ts1_ref --group-size 32"
      "micm-vector8-g32=$W/bin/micm-vector8-in-place-mozart --group-size 32")
for run in astra-max astra-xhigh astra-high sol-max sol-xhigh sol-high luna-max luna-xhigh luna-high; do
  args+=("$run/specialized=$W/verify/$run/clang/$run/bin/ts1_bench_specialized"
         "$run/generic=$W/verify/$run/clang/$run/bin/ts1_bench_cpu")
done
exec zsh $HOME/NSF-NCAR/TurboMICM/admin/time_tier2.sh $W/timing/ranking-1 1 5 $args
