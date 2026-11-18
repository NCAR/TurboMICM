| Rank | Entrant | Configuration | Cell groups | Median ms (range) | Runs | vs MICM default | vs fastest MICM | Gate | Largest relative difference |
|---:|---|---|---|---:|---:|---:|---:|---|---:|
| 1 | **GPT-6.1-Sol Max** | the agent's C port (`sol-max`) | 32 cells | 37,661 (37,573–37,820) | 25 | 🔶 **3.84×** | 1.35× | PASS | 1e-13 |
| 2 | **GPT-6-Astra Extra High** | the agent's C port (`astra-xhigh`) | 32 cells | 38,330 (38,095–38,514) | 25 | 🔶 **3.77×** | 1.33× | PASS | 2e-13 |
| 3 | **GPT-6.1-Terra Max** | the agent's C port (`terra-max`) | 64 cells | 43,691 (43,676–43,762) | 5 | 🔶 **3.31×** | 1.17× | PASS | 3e-13 |
| 4 | **GPT-6.1-Sol Extra High** | the agent's C port (`sol-xhigh`) | 32 cells | 44,202 (44,026–44,441) | 25 | 🔶 **3.27×** | 1.15× | PASS | 0 |
| 5 | **GPT-6-Astra Max** | the agent's C port (`astra-max`) | 32 cells | 45,390 (45,271–45,629) | 20 | 🔶 **3.19×** | 1.12× | PASS | 3e-13 |
| 6 | **GPT-6.1-Sol High** | the agent's C port (`sol-high`) | 32 cells | 46,791 (46,640–46,983) | 25 | 🔶 **3.09×** | 1.09× | PASS | 2e-13 |
| 7 | **GPT-6-Astra High** | the agent's C port (`astra-high`) | 32 cells | 49,322 (49,180–49,687) | 25 | 🔶 **2.93×** | 1.03× | PASS | 0 |
| 8 | MICM | vector8, in-place Mozart | 32 cells | 50,950 (50,064–51,211) | 25 | 🔶 **2.84×** | 1.00× | PASS | 0 |
| 9 | **GPT-6-Luna Extra High** | the agent's C port (`luna-xhigh`) | one per thread | 60,153 (59,903–60,344) | 10 | 🔶 **2.41×** | 0.85× | PASS | 2e-13 |
| 10 | MICM | standard, in-place Mozart | 32 cells | 135,790 (133,882–139,733) | 25 | 🔶 **1.07×** | 0.38× | PASS | 0 |
| 11 | **GPT-6.1-Terra High** | the agent's C port (`terra-high`) | one cell | 144,516 (143,201–144,696) | 5 | 🔶 **1.00×** | 0.35× | PASS | 1e-13 |
| 12 | MICM | standard, in-place Mozart | one per thread | 144,668 (143,680–146,770) | 25 | 🔶 **1.00×** | 0.35× | PASS | 0 |
| 13 | **GPT-6.1-Terra Extra High** | the agent's C port (`terra-xhigh`) | one cell | 168,884 (166,591–169,035) | 5 | 🔶 **0.86×** | 0.30× | PASS | 2e-13 |
| 14 | **GPT-6-Luna Max** | the agent's C port (`luna-max`) | one per thread | 182,866 (182,291–184,055) | 15 | 🔶 **0.79×** | 0.28× | PASS | 3e-13 |
| 15 | **GPT-6-Luna High** | the agent's C port (`luna-high`) | one per thread | 311,252 (308,706–312,341) | 5 | 🔶 **0.46×** | 0.16× | PASS | 5e-14 |
