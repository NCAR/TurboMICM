| Rank | Entrant | Configuration | Cell groups | Median ms (range) | Runs | vs MICM default | vs fastest MICM | Gate | Largest relative difference |
|---:|---|---|---|---:|---:|---:|---:|---|---:|
| 1 | **GPT-6.1-Sol Max** | the agent's C port (`sol-max`) | 32 cells | 37,684 (37,610–37,698) | 5 | 🔶 **3.84×** | 1.34× | PASS | 1e-13 |
| 2 | **GPT-6-Astra Extra High** | the agent's C port (`astra-xhigh`) | 32 cells | 38,328 (38,095–38,365) | 5 | 🔶 **3.77×** | 1.32× | PASS | 2e-13 |
| 3 | **GPT-6.1-Sol Extra High** | the agent's C port (`sol-xhigh`) | 32 cells | 44,105 (44,079–44,252) | 5 | 🔶 **3.28×** | 1.15× | PASS | 0 |
| 4 | **GPT-6.1-Sol High** | the agent's C port (`sol-high`) | 32 cells | 46,702 (46,663–46,828) | 5 | 🔶 **3.10×** | 1.08× | PASS | 2e-13 |
| 5 | **GPT-6-Astra High** | the agent's C port (`astra-high`) | 32 cells | 49,247 (49,180–49,561) | 5 | 🔶 **2.94×** | 1.03× | PASS | 0 |
| 6 | MICM | vector8, in-place Mozart | 32 cells | 50,659 (50,087–50,959) | 5 | 🔶 **2.86×** | 1.00× | PASS | 0 |
| 7 | MICM | standard, in-place Mozart | 32 cells | 136,210 (134,031–136,927) | 5 | 🔶 **1.06×** | 0.37× | PASS | 0 |
| 8 | MICM | standard, in-place Mozart | one per thread | 144,668 (143,832–145,461) | 5 | 🔶 **1.00×** | 0.35× | PASS | 0 |
