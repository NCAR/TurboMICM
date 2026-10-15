Median elapsed_ms (speed relative to the default), one cell group per thread:

| Matrix layout | in-place Mozart (default) | in-place Doolittle | separate Mozart | separate Doolittle |
|---|---:|---:|---:|---:|
| standard | 27,559 (1.00×) | 27,880 (0.99×) | 28,361 (0.97×) | 28,487 (0.97×) |
| vector1 |  | 27,358 (1.01×) | 27,910 (0.99×) | 28,316 (0.97×) |
| vector2 |  | 22,115 (1.25×) | 22,908 (1.20×) | 23,531 (1.17×) |
| vector4 |  | 19,076 (1.44×) | 20,311 (1.36×) | 20,241 (1.36×) |
| vector8 | 18,226 (1.51×) | 18,014 (1.53×) | 19,123 (1.44×) | 19,027 (1.45×) |
| vector128 |  | 19,672 (1.40×) | 22,085 (1.25×) | 20,914 (1.32×) |

Median elapsed_ms (speed relative to the default), 32-cell groups:

| Matrix layout | in-place Mozart (default) | in-place Doolittle | separate Mozart | separate Doolittle |
|---|---:|---:|---:|---:|
| standard | 22,787 (1.00×) | 22,759 (1.00×) | 24,694 (0.92×) | 24,231 (0.94×) |
| vector1 |  | 22,443 (1.02×) | 24,114 (0.94×) | 24,145 (0.94×) |
| vector2 |  | 15,361 (1.48×) | 16,640 (1.37×) | 17,056 (1.34×) |
| vector4 |  | 11,309 (2.01×) | 12,815 (1.78×) | 12,873 (1.77×) |
| vector8 | 9,704 (2.35×) | 9,819 (2.32×) | 10,890 (2.09×) | 10,990 (2.07×) |
| vector128 |  | 36,919 (0.62×) | 49,635 (0.46×) | 47,392 (0.48×) |

Largest relative difference of the final states from MICM_FINAL, over both groupings:

| Matrix layout | in-place Mozart (default) | in-place Doolittle | separate Mozart | separate Doolittle |
|---|---:|---:|---:|---:|
| standard | 0 | 1e-13 | 0 | 1e-13 |
| vector1 |  | 1e-13 | 0 | 1e-13 |
| vector2 |  | 1e-13 | 0 | 1e-13 |
| vector4 |  | 1e-13 | 0 | 1e-13 |
| vector8 | 0 | 1e-13 | 0 | 1e-13 |
| vector128 |  | 1e-13 | 0 | 1e-13 |

| Configuration | One group per thread: median ms (range) | vs default | 32-cell groups: median ms (range) | vs default | Internal steps (group per thread; 32-cell) | Gate | Largest relative difference |
|---|---:|---:|---:|---:|---:|---|---:|
| standard-in-place-mozart | 27,559 (27,013–28,304) | 1.00× | 22,787 (22,585–23,722) | 1.00× | 9000; 281700 | PASS | 0 |
| standard-in-place-doolittle | 27,880 (26,910–28,158) | 0.99× | 22,759 (22,682–23,522) | 1.00× | 9000; 281700 | PASS | 1.1e-13 |
| standard-separate-mozart | 28,361 (28,036–29,072) | 0.97× | 24,694 (23,960–24,803) | 0.92× | 9000; 281700 | PASS | 0 |
| standard-separate-doolittle | 28,487 (27,647–29,135) | 0.97× | 24,231 (23,813–25,295) | 0.94× | 9000; 281700 | PASS | 1.1e-13 |
| vector1-in-place-doolittle | 27,358 (27,023–30,968) | 1.01× | 22,443 (22,224–23,607) | 1.02× | 9000; 281700 | PASS | 1.1e-13 |
| vector1-separate-mozart | 27,910 (27,788–29,573) | 0.99× | 24,114 (23,924–24,464) | 0.94× | 9000; 281700 | PASS | 0 |
| vector1-separate-doolittle | 28,316 (27,890–28,695) | 0.97× | 24,145 (23,885–28,070) | 0.94× | 9000; 281700 | PASS | 1.1e-13 |
| vector2-in-place-doolittle | 22,115 (22,098–23,668) | 1.25× | 15,361 (15,060–15,529) | 1.48× | 9000; 281700 | PASS | 1.1e-13 |
| vector2-separate-mozart | 22,908 (22,723–28,814) | 1.20× | 16,640 (16,561–21,171) | 1.37× | 9000; 281700 | PASS | 0 |
| vector2-separate-doolittle | 23,531 (23,294–30,200) | 1.17× | 17,056 (16,843–20,426) | 1.34× | 9000; 281700 | PASS | 1.1e-13 |
| vector4-in-place-doolittle | 19,076 (19,076–25,873) | 1.44× | 11,309 (11,112–15,230) | 2.01× | 9000; 281700 | PASS | 1.1e-13 |
| vector4-separate-mozart | 20,311 (20,125–27,834) | 1.36× | 12,815 (12,558–15,735) | 1.78× | 9000; 281700 | PASS | 0 |
| vector4-separate-doolittle | 20,241 (20,119–25,346) | 1.36× | 12,873 (12,822–15,256) | 1.77× | 9000; 281700 | PASS | 1.1e-13 |
| vector8-in-place-mozart | 18,226 (18,087–18,246) | 1.51× | 9,704 (9,701–9,776) | 2.35× | 9000; 281700 | PASS | 0 |
| vector8-in-place-doolittle | 18,014 (17,758–18,558) | 1.53× | 9,819 (9,713–9,999) | 2.32× | 9000; 281700 | PASS | 1.1e-13 |
| vector8-separate-mozart | 19,123 (19,067–19,466) | 1.44× | 10,890 (10,823–10,968) | 2.09× | 9000; 281700 | PASS | 0 |
| vector8-separate-doolittle | 19,027 (18,995–19,090) | 1.45× | 10,990 (10,844–11,020) | 2.07× | 9000; 281700 | PASS | 1.1e-13 |
| vector128-in-place-doolittle | 19,672 (19,635–20,508) | 1.40× | 36,919 (36,657–38,242) | 0.62× | 9000; 281700 | PASS | 1.1e-13 |
| vector128-separate-mozart | 22,085 (22,015–22,090) | 1.25× | 49,635 (49,204–50,984) | 0.46× | 9000; 281700 | PASS | 0 |
| vector128-separate-doolittle | 20,914 (20,795–20,935) | 1.32× | 47,392 (46,862–47,914) | 0.48× | 9000; 281700 | PASS | 1.1e-13 |
