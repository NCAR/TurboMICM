| Rank | Entrant | Configuration | Cell groups | Median ms (range) | Runs | vs MICM default | vs fastest MICM | Gate | Largest relative difference |
|---:|---|---|---|---:|---:|---:|---:|---|---:|
| 1 | **GPT-6.1-Sol Max** | the agent's C port (`sol-max`) | 32 cells | 37,629 (37,573–37,698) | 15 | 🔶 **3.84×** | 1.35× | PASS | 1e-13 |
| 2 | **GPT-6-Astra Extra High** | the agent's C port (`astra-xhigh`) | 32 cells | 38,365 (38,095–38,514) | 15 | 🔶 **3.77×** | 1.32× | PASS | 2e-13 |
| 3 | **GPT-6.1-Sol Extra High** | the agent's C port (`sol-xhigh`) | 32 cells | 44,164 (44,026–44,298) | 15 | 🔶 **3.27×** | 1.15× | PASS | 0 |
| 4 | **GPT-6-Astra Max** | the agent's C port (`astra-max`) | 32 cells | 45,390 (45,294–45,488) | 10 | 🔶 **3.18×** | 1.12× | PASS | 3e-13 |
| 5 | **GPT-6.1-Sol High** | the agent's C port (`sol-high`) | 32 cells | 46,769 (46,640–46,911) | 15 | 🔶 **3.09×** | 1.08× | PASS | 2e-13 |
| 6 | **GPT-6-Astra High** | the agent's C port (`astra-high`) | 32 cells | 49,258 (49,180–49,611) | 15 | 🔶 **2.93×** | 1.03× | PASS | 0 |
| 7 | MICM | vector8, in-place Mozart | 32 cells | 50,673 (50,087–51,211) | 15 | 🔶 **2.85×** | 1.00× | PASS | 0 |
| 8 | MICM | standard, in-place Mozart | 32 cells | 135,739 (133,882–136,927) | 15 | 🔶 **1.06×** | 0.37× | PASS | 0 |
| 9 | MICM | standard, in-place Mozart | one per thread | 144,522 (143,714–145,461) | 15 | 🔶 **1.00×** | 0.35× | PASS | 0 |
| 10 | **GPT-6-Luna Max** | the agent's C port (`luna-max`) | one per thread | 183,244 (182,292–184,055) | 5 | 🔶 **0.79×** | 0.28× | PASS | 3e-13 |
