| Rank | Entrant | Configuration | Cell groups | Median ms (range) | Runs | vs MICM default | vs fastest MICM | Gate | Largest relative difference |
|---:|---|---|---|---:|---:|---:|---:|---|---:|
| 1 | **GPT-6.1-Sol Max** | the agent's C port (`sol-max`) | 32 cells | 37,633 (37,573–37,698) | 20 | 🔶 **3.84×** | 1.35× | PASS | 1e-13 |
| 2 | **GPT-6-Astra Extra High** | the agent's C port (`astra-xhigh`) | 32 cells | 38,332 (38,095–38,514) | 20 | 🔶 **3.77×** | 1.33× | PASS | 2e-13 |
| 3 | **GPT-6.1-Sol Extra High** | the agent's C port (`sol-xhigh`) | 32 cells | 44,177 (44,026–44,298) | 20 | 🔶 **3.27×** | 1.15× | PASS | 0 |
| 4 | **GPT-6-Astra Max** | the agent's C port (`astra-max`) | 32 cells | 45,376 (45,271–45,488) | 15 | 🔶 **3.18×** | 1.12× | PASS | 3e-13 |
| 5 | **GPT-6.1-Sol High** | the agent's C port (`sol-high`) | 32 cells | 46,790 (46,640–46,911) | 20 | 🔶 **3.09×** | 1.09× | PASS | 2e-13 |
| 6 | **GPT-6-Astra High** | the agent's C port (`astra-high`) | 32 cells | 49,285 (49,180–49,611) | 20 | 🔶 **2.93×** | 1.03× | PASS | 0 |
| 7 | MICM | vector8, in-place Mozart | 32 cells | 50,824 (50,064–51,211) | 20 | 🔶 **2.84×** | 1.00× | PASS | 0 |
| 8 | **GPT-6-Luna Extra High** | the agent's C port (`luna-xhigh`) | one per thread | 59,916 (59,903–59,993) | 5 | 🔶 **2.41×** | 0.85× | PASS | 2e-13 |
| 9 | MICM | standard, in-place Mozart | 32 cells | 135,764 (133,882–136,927) | 20 | 🔶 **1.06×** | 0.37× | PASS | 0 |
| 10 | MICM | standard, in-place Mozart | one per thread | 144,491 (143,680–145,585) | 20 | 🔶 **1.00×** | 0.35× | PASS | 0 |
| 11 | **GPT-6-Luna Max** | the agent's C port (`luna-max`) | one per thread | 183,273 (182,292–184,055) | 10 | 🔶 **0.79×** | 0.28× | PASS | 3e-13 |
