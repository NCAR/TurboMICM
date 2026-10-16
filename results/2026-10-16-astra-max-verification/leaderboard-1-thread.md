| Rank | Entrant | Configuration | Cell groups | Median ms (range) | Runs | vs MICM default | vs fastest MICM | Gate | Largest relative difference |
|---:|---|---|---|---:|---:|---:|---:|---|---:|
| 1 | **GPT-6.1-Sol Max** | the agent's C port (`sol-max`) | 32 cells | 37,645 (37,598–37,698) | 10 | 🔶 **3.84×** | 1.35× | PASS | 1e-13 |
| 2 | **GPT-6-Astra Extra High** | the agent's C port (`astra-xhigh`) | 32 cells | 38,374 (38,095–38,514) | 10 | 🔶 **3.77×** | 1.32× | PASS | 2e-13 |
| 3 | **GPT-6.1-Sol Extra High** | the agent's C port (`sol-xhigh`) | 32 cells | 44,124 (44,026–44,298) | 10 | 🔶 **3.28×** | 1.15× | PASS | 0 |
| 4 | **GPT-6-Astra Max** | the agent's C port (`astra-max`) | 32 cells | 45,388 (45,296–45,457) | 5 | 🔶 **3.18×** | 1.12× | PASS | 3e-13 |
| 5 | **GPT-6.1-Sol High** | the agent's C port (`sol-high`) | 32 cells | 46,779 (46,663–46,911) | 10 | 🔶 **3.09×** | 1.08× | PASS | 2e-13 |
| 6 | **GPT-6-Astra High** | the agent's C port (`astra-high`) | 32 cells | 49,251 (49,180–49,561) | 10 | 🔶 **2.93×** | 1.03× | PASS | 0 |
| 7 | MICM | vector8, in-place Mozart | 32 cells | 50,681 (50,087–51,070) | 10 | 🔶 **2.85×** | 1.00× | PASS | 0 |
| 8 | MICM | standard, in-place Mozart | 32 cells | 135,880 (133,920–136,927) | 10 | 🔶 **1.06×** | 0.37× | PASS | 0 |
| 9 | MICM | standard, in-place Mozart | one per thread | 144,509 (143,714–145,461) | 10 | 🔶 **1.00×** | 0.35× | PASS | 0 |
