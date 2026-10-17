# TS1 benchmark checks

10000 cells (seed 2026), 100 steps of 30 s; the gate is rtol 1e-10, atol 1e-30 against MICM.

| Variant | Arguments | Median ms | Runs | Gate | Largest relative difference | Statistics |
|---|---|---:|---:|---|---:|---|
| default (all threads, groups of 32) | `` | 4,577 | 1 | PASS | 1.04e-13 | internal_steps=281700 accepted=281700 rejected=0 function_calls=563400 jacobian_updates=281700 decompositions=281700 solves=845100 |
| 1 thread | `--threads 1` | 25,950 | 1 | PASS | 1.04e-13 | internal_steps=281700 accepted=281700 rejected=0 function_calls=563400 jacobian_updates=281700 decompositions=281700 solves=845100 |
| groups of 4 | `--group-size 4` | 6,395 | 1 | PASS | 1.04e-13 | internal_steps=2250000 accepted=2250000 rejected=0 function_calls=4500000 jacobian_updates=2250000 decompositions=2250000 solves=6750000 |
| groups of 8 | `--group-size 8` | 4,950 | 1 | PASS | 1.04e-13 | internal_steps=1125000 accepted=1125000 rejected=0 function_calls=2250000 jacobian_updates=1125000 decompositions=1125000 solves=3375000 |
| groups of 16 | `--group-size 16` | 4,876 | 1 | PASS | 1.04e-13 | internal_steps=562500 accepted=562500 rejected=0 function_calls=1125000 jacobian_updates=562500 decompositions=562500 solves=1687500 |
| groups of 64 | `--group-size 64` | 5,464 | 1 | PASS | 1.04e-13 | internal_steps=141300 accepted=141300 rejected=0 function_calls=282600 jacobian_updates=141300 decompositions=141300 solves=423900 |
| generic kernels | `--generic` | 27,829 | 1 | PASS | 0 | internal_steps=281700 accepted=281700 rejected=0 function_calls=563400 jacobian_updates=281700 decompositions=281700 solves=845100 |
| static schedule | `--schedule static` | 4,917 | 1 | PASS | 1.04e-13 | internal_steps=281700 accepted=281700 rejected=0 function_calls=563400 jacobian_updates=281700 decompositions=281700 solves=845100 |
| guided schedule | `--schedule guided` | 4,866 | 1 | PASS | 1.04e-13 | internal_steps=281700 accepted=281700 rejected=0 function_calls=563400 jacobian_updates=281700 decompositions=281700 solves=845100 |
| 2 threads, groups of 3 | `--threads 2 --group-size 3` | 49,481 | 1 | PASS | 1.04e-13 | internal_steps=3000600 accepted=3000600 rejected=0 function_calls=6001200 jacobian_updates=3000600 decompositions=3000600 solves=9001800 |
