Median elapsed_ms (speed relative to the default), one cell group per thread:

| Matrix layout | in-place Mozart (default) | in-place Doolittle | separate Mozart | separate Doolittle |
|---|---:|---:|---:|---:|
| standard | 2,532 (1.00×) |  |  |  |
| vector1 |  |  |  |  |
| vector2 | 1,943 (1.30×) | 2,157 (1.17×) |  |  |
| vector4 | 1,526 (1.66×) | 1,590 (1.59×) |  |  |
| vector8 | 1,403 (1.80×) | 1,448 (1.75×) | 2,334 (1.08×) |  |
| vector128 | 11,058 (0.23×) |  |  |  |

Median elapsed_ms (speed relative to the default), 32-cell groups:

| Matrix layout | in-place Mozart (default) | in-place Doolittle | separate Mozart | separate Doolittle |
|---|---:|---:|---:|---:|
| standard | 3,076 (1.00×) |  |  |  |
| vector1 |  |  |  |  |
| vector2 | 2,286 (1.35×) | 2,641 (1.16×) |  |  |
| vector4 | 1,735 (1.77×) | 1,858 (1.66×) |  |  |
| vector8 | 1,574 (1.95×) | 1,645 (1.87×) | 1,798 (1.71×) |  |
| vector128 | 11,911 (0.26×) |  |  |  |

Largest relative difference of the final states from MICM_FINAL, over both groupings:

| Matrix layout | in-place Mozart (default) | in-place Doolittle | separate Mozart | separate Doolittle |
|---|---:|---:|---:|---:|
| standard | 0 |  |  |  |
| vector1 |  |  |  |  |
| vector2 | 1e-13 | 2e-13 |  |  |
| vector4 | 1e-13 | 1e-13 |  |  |
| vector8 | 1e-13 | 1e-13 | 1e-13 |  |
| vector128 | 1e-13 |  |  |  |

| Configuration | One group per thread: median ms (range) | vs default | 32-cell groups: median ms (range) | vs default | Internal steps (group per thread; 32-cell) | Gate | Largest relative difference |
|---|---:|---:|---:|---:|---:|---|---:|
| standard-in-place-mozart | 2,532 (2,527–2,598) | 1.00× | 3,076 (3,053–3,086) | 1.00× | 230400; 281700 | PASS | 0 |
| vector2-in-place-mozart | 1,943 (1,941–1,943) | 1.30× | 2,286 (2,273–2,287) | 1.35× | 230400; 281700 | PASS | 9.9e-14 |
| vector2-in-place-doolittle | 2,157 (2,150–2,208) | 1.17× | 2,641 (2,612–2,647) | 1.16× | 230400; 281700 | PASS | 1.6e-13 |
| vector4-in-place-mozart | 1,526 (1,525–1,535) | 1.66× | 1,735 (1,728–1,773) | 1.77× | 230400; 281700 | PASS | 1.4e-13 |
| vector4-in-place-doolittle | 1,590 (1,587–1,594) | 1.59× | 1,858 (1,850–1,882) | 1.66× | 230400; 281700 | PASS | 9.7e-14 |
| vector8-in-place-mozart | 1,403 (1,401–1,407) | 1.80× | 1,574 (1,572–1,593) | 1.95× | 230400; 281700 | PASS | 1e-13 |
| vector8-in-place-doolittle | 1,448 (1,443–1,457) | 1.75× | 1,645 (1,641–1,667) | 1.87× | 230400; 281700 | PASS | 1e-13 |
| vector8-separate-mozart | 2,334 (2,286–2,366) | 1.08× | 1,798 (1,759–1,829) | 1.71× | 230400; 281700 | PASS | 1e-13 |
| vector128-in-place-mozart | 11,058 (10,607–11,116) | 0.23× | 11,911 (11,905–12,083) | 0.26× | 230400; 281700 | PASS | 1.4e-13 |
