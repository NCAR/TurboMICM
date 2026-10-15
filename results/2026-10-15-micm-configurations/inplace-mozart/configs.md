Median elapsed_ms (speed relative to the default), one cell group per thread:

| Matrix layout | in-place Mozart (default) | in-place Doolittle | separate Mozart | separate Doolittle |
|---|---:|---:|---:|---:|
| standard | 27,752 (1.00×) |  |  |  |
| vector1 | 27,357 (1.01×) |  |  |  |
| vector2 | 22,380 (1.24×) |  |  |  |
| vector4 | 19,354 (1.43×) |  |  |  |
| vector8 | 18,158 (1.53×) |  |  |  |
| vector128 | 19,627 (1.41×) |  |  |  |

Median elapsed_ms (speed relative to the default), 32-cell groups:

| Matrix layout | in-place Mozart (default) | in-place Doolittle | separate Mozart | separate Doolittle |
|---|---:|---:|---:|---:|
| standard | 23,485 (1.00×) |  |  |  |
| vector1 | 22,242 (1.06×) |  |  |  |
| vector2 | 15,311 (1.53×) |  |  |  |
| vector4 | 11,362 (2.07×) |  |  |  |
| vector8 | 9,826 (2.39×) |  |  |  |
| vector128 | 36,531 (0.64×) |  |  |  |

Largest relative difference of the final states from MICM_FINAL, over both groupings:

| Matrix layout | in-place Mozart (default) | in-place Doolittle | separate Mozart | separate Doolittle |
|---|---:|---:|---:|---:|
| standard | 0 |  |  |  |
| vector1 | 0 |  |  |  |
| vector2 | 0 |  |  |  |
| vector4 | 0 |  |  |  |
| vector8 | 0 |  |  |  |
| vector128 | 0 |  |  |  |

| Configuration | One group per thread: median ms (range) | vs default | 32-cell groups: median ms (range) | vs default | Internal steps (group per thread; 32-cell) | Gate | Largest relative difference |
|---|---:|---:|---:|---:|---:|---|---:|
| standard-in-place-mozart | 27,752 (27,750–28,144) | 1.00× | 23,485 (23,475–24,502) | 1.00× | 9000; 281700 | PASS | 0 |
| vector1-in-place-mozart | 27,357 (26,805–27,649) | 1.01× | 22,242 (22,221–22,894) | 1.06× | 9000; 281700 | PASS | 0 |
| vector2-in-place-mozart | 22,380 (22,327–22,543) | 1.24× | 15,311 (15,144–15,446) | 1.53× | 9000; 281700 | PASS | 0 |
| vector4-in-place-mozart | 19,354 (19,133–19,723) | 1.43× | 11,362 (11,272–11,692) | 2.07× | 9000; 281700 | PASS | 0 |
| vector8-in-place-mozart | 18,158 (18,010–18,405) | 1.53× | 9,826 (9,727–10,135) | 2.39× | 9000; 281700 | PASS | 0 |
| vector128-in-place-mozart | 19,627 (19,431–19,803) | 1.41× | 36,531 (35,340–37,389) | 0.64× | 9000; 281700 | PASS | 0 |

| Program | Median ms (range) | vs default, one group per thread | vs default, 32-cell groups | Gate | Largest relative difference |
|---|---:|---:|---:|---|---:|
| sol-high | 8,681 (8,610–8,759) | 3.20× | 2.71× | PASS | 2e-13 |
| sol-xhigh | 8,196 (8,112–8,250) | 3.39× | 2.87× | PASS | 0 |
| sol-max | 7,160 (7,135–7,166) | 3.88× | 3.28× | PASS | 1.3e-13 |
| astra-high | 9,265 (9,200–9,299) | 3.00× | 2.53× | PASS | 0 |
