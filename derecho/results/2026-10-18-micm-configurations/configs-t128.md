Median elapsed_ms (speed relative to the default), one cell group per thread:

| Matrix layout | in-place Mozart (default) | in-place Doolittle | separate Mozart | separate Doolittle |
|---|---:|---:|---:|---:|
| standard | 3,528 (1.00×) | 3,494 (1.01×) | 4,158 (0.85×) | 3,987 (0.88×) |
| vector1 | 3,609 (0.98×) | 3,686 (0.96×) | 4,128 (0.85×) | 3,859 (0.91×) |
| vector2 | 3,006 (1.17×) | 3,360 (1.05×) | 3,204 (1.10×) | 3,715 (0.95×) |
| vector4 | 2,198 (1.61×) | 2,361 (1.49×) | 2,842 (1.24×) | 2,811 (1.26×) |
| vector8 | 1,803 (1.96×) | 1,954 (1.81×) | 2,555 (1.38×) | 2,849 (1.24×) |
| vector128 | 2,716 (1.30×) | 2,886 (1.22×) | 7,259 (0.49×) | 7,010 (0.50×) |

Median elapsed_ms (speed relative to the default), 32-cell groups:

| Matrix layout | in-place Mozart (default) | in-place Doolittle | separate Mozart | separate Doolittle |
|---|---:|---:|---:|---:|
| standard | 3,813 (1.00×) | 3,872 (0.98×) | 4,343 (0.88×) | 4,039 (0.94×) |
| vector1 | 3,508 (1.09×) | 3,562 (1.07×) | 3,881 (0.98×) | 4,085 (0.93×) |
| vector2 | 3,181 (1.20×) | 3,416 (1.12×) | 2,595 (1.47×) | 3,718 (1.03×) |
| vector4 | 2,245 (1.70×) | 2,258 (1.69×) | 2,204 (1.73×) | 2,478 (1.54×) |
| vector8 | 1,864 (2.05×) | 1,865 (2.04×) | 2,080 (1.83×) | 2,194 (1.74×) |
| vector128 | 5,679 (0.67×) | 5,856 (0.65×) | 14,238 (0.27×) | 14,329 (0.27×) |

Largest relative difference of the final states from MICM_FINAL, over both groupings:

| Matrix layout | in-place Mozart (default) | in-place Doolittle | separate Mozart | separate Doolittle |
|---|---:|---:|---:|---:|
| standard | 0 | 9e-14 | 0 | 9e-14 |
| vector1 | 0 | 9e-14 | 0 | 9e-14 |
| vector2 | 1e-13 | 2e-13 | 1e-13 | 2e-13 |
| vector4 | 1e-13 | 1e-13 | 1e-13 | 1e-13 |
| vector8 | 1e-13 | 1e-13 | 1e-13 | 1e-13 |
| vector128 | 1e-13 | 1e-13 | 1e-13 | 1e-13 |

| Configuration | One group per thread: median ms (range) | vs default | 32-cell groups: median ms (range) | vs default | Internal steps (group per thread; 32-cell) | Gate | Largest relative difference |
|---|---:|---:|---:|---:|---:|---|---:|
| standard-in-place-mozart | 3,528 (3,188–3,608) | 1.00× | 3,813 (3,791–3,921) | 1.00× | 115200; 281700 | PASS | 0 |
| standard-in-place-doolittle | 3,494 (3,486–3,494) | 1.01× | 3,872 (3,457–3,882) | 0.98× | 115200; 281700 | PASS | 9.3e-14 |
| standard-separate-mozart | 4,158 (3,893–4,318) | 0.85× | 4,343 (3,870–4,388) | 0.88× | 115200; 281700 | PASS | 0 |
| standard-separate-doolittle | 3,987 (3,930–4,242) | 0.88× | 4,039 (3,657–4,255) | 0.94× | 115200; 281700 | PASS | 9.3e-14 |
| vector1-in-place-mozart | 3,609 (3,289–3,629) | 0.98× | 3,508 (3,434–3,824) | 1.09× | 115200; 281700 | PASS | 0 |
| vector1-in-place-doolittle | 3,686 (2,989–3,699) | 0.96× | 3,562 (3,328–3,827) | 1.07× | 115200; 281700 | PASS | 9.3e-14 |
| vector1-separate-mozart | 4,128 (3,840–4,187) | 0.85× | 3,881 (3,805–4,192) | 0.98× | 115200; 281700 | PASS | 0 |
| vector1-separate-doolittle | 3,859 (3,793–4,252) | 0.91× | 4,085 (3,968–4,233) | 0.93× | 115200; 281700 | PASS | 9.3e-14 |
| vector2-in-place-mozart | 3,006 (2,974–3,029) | 1.17× | 3,181 (3,096–3,248) | 1.20× | 115200; 281700 | PASS | 9.9e-14 |
| vector2-in-place-doolittle | 3,360 (2,932–3,440) | 1.05× | 3,416 (3,281–3,702) | 1.12× | 115200; 281700 | PASS | 1.6e-13 |
| vector2-separate-mozart | 3,204 (3,023–3,221) | 1.10× | 2,595 (2,553–2,933) | 1.47× | 115200; 281700 | PASS | 1e-13 |
| vector2-separate-doolittle | 3,715 (3,671–3,809) | 0.95× | 3,718 (3,322–3,810) | 1.03× | 115200; 281700 | PASS | 1.6e-13 |
| vector4-in-place-mozart | 2,198 (2,180–2,208) | 1.61× | 2,245 (2,240–2,264) | 1.70× | 115200; 281700 | PASS | 1.4e-13 |
| vector4-in-place-doolittle | 2,361 (2,318–2,377) | 1.49× | 2,258 (2,061–2,468) | 1.69× | 115200; 281700 | PASS | 9.7e-14 |
| vector4-separate-mozart | 2,842 (2,698–2,911) | 1.24× | 2,204 (2,191–2,210) | 1.73× | 115200; 281700 | PASS | 1.5e-13 |
| vector4-separate-doolittle | 2,811 (2,697–2,959) | 1.26× | 2,478 (2,477–2,580) | 1.54× | 115200; 281700 | PASS | 9.7e-14 |
| vector8-in-place-mozart | 1,803 (1,659–1,929) | 1.96× | 1,864 (1,841–1,985) | 2.05× | 115200; 281700 | PASS | 1e-13 |
| vector8-in-place-doolittle | 1,954 (1,895–2,051) | 1.81× | 1,865 (1,723–1,961) | 2.04× | 115200; 281700 | PASS | 1e-13 |
| vector8-separate-mozart | 2,555 (2,516–2,646) | 1.38× | 2,080 (2,047–2,102) | 1.83× | 115200; 281700 | PASS | 1e-13 |
| vector8-separate-doolittle | 2,849 (2,644–3,071) | 1.24× | 2,194 (2,165–2,197) | 1.74× | 115200; 281700 | PASS | 1e-13 |
| vector128-in-place-mozart | 2,716 (2,677–2,801) | 1.30× | 5,679 (5,463–5,816) | 0.67× | 115200; 281700 | PASS | 1.4e-13 |
| vector128-in-place-doolittle | 2,886 (2,817–3,079) | 1.22× | 5,856 (5,346–6,204) | 0.65× | 115200; 281700 | PASS | 9.7e-14 |
| vector128-separate-mozart | 7,259 (6,865–7,383) | 0.49× | 14,238 (14,177–14,260) | 0.27× | 115200; 281700 | PASS | 1.4e-13 |
| vector128-separate-doolittle | 7,010 (6,915–7,364) | 0.50× | 14,329 (14,239–14,499) | 0.27× | 115200; 281700 | PASS | 9.7e-14 |
