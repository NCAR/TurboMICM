| Rank | Entrant | Configuration | Cell groups | Median ms (range) | Runs | vs MICM default | vs fastest MICM | Gate | Largest relative difference |
|---:|---|---|---|---:|---:|---:|---:|---|---:|
| 1 | **GPT-6.1-Sol Max** | the agent's C port (`sol-max`) | 32 cells | 7,136 (6,942–7,266) | 23 | 🔶 **3.89×** | 1.37× | PASS | 1e-13 |
| 2 | **GPT-6-Astra Extra High** | the agent's C port (`astra-xhigh`) | 32 cells | 7,296 (7,095–7,478) | 20 | 🔶 **3.81×** | 1.34× | PASS | 2e-13 |
| 3 | **GPT-6.1-Sol Extra High** | the agent's C port (`sol-xhigh`) | 32 cells | 8,196 (7,886–8,359) | 23 | 🔶 **3.39×** | 1.19× | PASS | 0 |
| 4 | **GPT-6-Astra Max** | the agent's C port (`astra-max`) | 32 cells | 8,408 (8,090–8,663) | 15 | 🔶 **3.31×** | 1.16× | PASS | 3e-13 |
| 5 | **GPT-6.1-Sol High** | the agent's C port (`sol-high`) | 32 cells | 8,681 (8,331–8,759) | 23 | 🔶 **3.20×** | 1.13× | PASS | 2e-13 |
| 6 | **GPT-6-Astra High** | the agent's C port (`astra-high`) | 32 cells | 9,101 (8,857–9,449) | 23 | 🔶 **3.05×** | 1.07× | PASS | 0 |
| 7 | MICM | vector8, in-place Mozart | 32 cells | 9,769 (9,195–10,135) | 26 | 🔶 **2.84×** | 1.00× | PASS | 0 |
| 8 | MICM | vector8, in-place Doolittle | 32 cells | 9,819 (9,713–9,999) | 3 | 🔶 **2.83×** | 0.99× | PASS | 1e-13 |
| 9 | MICM | vector8, separate Mozart | 32 cells | 10,890 (10,823–10,968) | 3 | 🔶 **2.55×** | 0.90× | PASS | 0 |
| 10 | MICM | vector8, separate Doolittle | 32 cells | 10,990 (10,844–11,020) | 3 | 🔶 **2.53×** | 0.89× | PASS | 1e-13 |
| 11 | MICM | vector4, in-place Doolittle | 32 cells | 11,309 (11,112–15,230) | 3 | 🔶 **2.46×** | 0.86× | PASS | 1e-13 |
| 12 | MICM | vector4, in-place Mozart | 32 cells | 11,362 (11,272–11,692) | 3 | 🔶 **2.45×** | 0.86× | PASS | 0 |
| 13 | MICM | vector4, separate Mozart | 32 cells | 12,815 (12,558–15,735) | 3 | 🔶 **2.17×** | 0.76× | PASS | 0 |
| 14 | MICM | vector4, separate Doolittle | 32 cells | 12,873 (12,822–15,256) | 3 | 🔶 **2.16×** | 0.76× | PASS | 1e-13 |
| 15 | MICM | vector2, in-place Mozart | 32 cells | 15,311 (15,144–15,446) | 3 | 🔶 **1.81×** | 0.64× | PASS | 0 |
| 16 | MICM | vector2, in-place Doolittle | 32 cells | 15,361 (15,060–15,529) | 3 | 🔶 **1.81×** | 0.64× | PASS | 1e-13 |
| 17 | MICM | vector2, separate Mozart | 32 cells | 16,640 (16,561–21,171) | 3 | 🔶 **1.67×** | 0.59× | PASS | 0 |
| 18 | MICM | vector2, separate Doolittle | 32 cells | 17,056 (16,843–20,426) | 3 | 🔶 **1.63×** | 0.57× | PASS | 1e-13 |
| 19 | MICM | vector8, in-place Doolittle | one per thread | 18,014 (17,758–18,558) | 3 | 🔶 **1.54×** | 0.54× | PASS | 1e-13 |
| 20 | MICM | vector8, in-place Mozart | one per thread | 18,192 (18,010–18,405) | 6 | 🔶 **1.53×** | 0.54× | PASS | 0 |
| 21 | MICM | vector8, separate Doolittle | one per thread | 19,027 (18,995–19,090) | 3 | 🔶 **1.46×** | 0.51× | PASS | 1e-13 |
| 22 | MICM | vector4, in-place Doolittle | one per thread | 19,076 (19,076–25,873) | 3 | 🔶 **1.46×** | 0.51× | PASS | 1e-13 |
| 23 | MICM | vector8, separate Mozart | one per thread | 19,123 (19,067–19,466) | 3 | 🔶 **1.45×** | 0.51× | PASS | 0 |
| 24 | MICM | vector4, in-place Mozart | one per thread | 19,354 (19,133–19,723) | 3 | 🔶 **1.44×** | 0.50× | PASS | 0 |
| 25 | MICM | vector128, in-place Mozart | one per thread | 19,627 (19,431–19,803) | 3 | 🔶 **1.42×** | 0.50× | PASS | 0 |
| 26 | MICM | vector128, in-place Doolittle | one per thread | 19,672 (19,635–20,508) | 3 | 🔶 **1.41×** | 0.50× | PASS | 1e-13 |
| 27 | MICM | vector4, separate Doolittle | one per thread | 20,241 (20,119–25,346) | 3 | 🔶 **1.37×** | 0.48× | PASS | 1e-13 |
| 28 | MICM | vector4, separate Mozart | one per thread | 20,311 (20,125–27,834) | 3 | 🔶 **1.37×** | 0.48× | PASS | 0 |
| 29 | MICM | vector128, separate Doolittle | one per thread | 20,914 (20,795–20,935) | 3 | 🔶 **1.33×** | 0.47× | PASS | 1e-13 |
| 30 | **GPT-6-Luna Extra High** | the agent's C port (`luna-xhigh`) | one per thread | 21,443 (21,373–21,506) | 5 | 🔶 **1.30×** | 0.46× | PASS | 2e-13 |
| 31 | MICM | vector128, separate Mozart | one per thread | 22,085 (22,015–22,090) | 3 | 🔶 **1.26×** | 0.44× | PASS | 0 |
| 32 | MICM | vector2, in-place Doolittle | one per thread | 22,115 (22,098–23,668) | 3 | 🔶 **1.26×** | 0.44× | PASS | 1e-13 |
| 33 | MICM | vector1, in-place Mozart | 32 cells | 22,242 (22,221–22,894) | 3 | 🔶 **1.25×** | 0.44× | PASS | 0 |
| 34 | MICM | vector2, in-place Mozart | one per thread | 22,380 (22,327–22,543) | 3 | 🔶 **1.24×** | 0.44× | PASS | 0 |
| 35 | MICM | vector1, in-place Doolittle | 32 cells | 22,443 (22,224–23,607) | 3 | 🔶 **1.24×** | 0.44× | PASS | 1e-13 |
| 36 | MICM | standard, in-place Doolittle | 32 cells | 22,759 (22,682–23,522) | 3 | 🔶 **1.22×** | 0.43× | PASS | 1e-13 |
| 37 | MICM | vector2, separate Mozart | one per thread | 22,908 (22,723–28,814) | 3 | 🔶 **1.21×** | 0.43× | PASS | 0 |
| 38 | MICM | standard, in-place Mozart | 32 cells | 23,414 (22,284–24,502) | 26 | 🔶 **1.19×** | 0.42× | PASS | 0 |
| 39 | MICM | vector2, separate Doolittle | one per thread | 23,531 (23,294–30,200) | 3 | 🔶 **1.18×** | 0.42× | PASS | 1e-13 |
| 40 | MICM | vector1, separate Mozart | 32 cells | 24,114 (23,924–24,464) | 3 | 🔶 **1.15×** | 0.41× | PASS | 0 |
| 41 | MICM | vector1, separate Doolittle | 32 cells | 24,145 (23,885–28,070) | 3 | 🔶 **1.15×** | 0.40× | PASS | 1e-13 |
| 42 | MICM | standard, separate Doolittle | 32 cells | 24,231 (23,813–25,295) | 3 | 🔶 **1.15×** | 0.40× | PASS | 1e-13 |
| 43 | MICM | standard, separate Mozart | 32 cells | 24,694 (23,960–24,803) | 3 | 🔶 **1.13×** | 0.40× | PASS | 0 |
| 44 | MICM | vector1, in-place Mozart | one per thread | 27,357 (26,805–27,649) | 3 | 🔶 **1.02×** | 0.36× | PASS | 0 |
| 45 | MICM | vector1, in-place Doolittle | one per thread | 27,358 (27,023–30,968) | 3 | 🔶 **1.02×** | 0.36× | PASS | 1e-13 |
| 46 | MICM | standard, in-place Mozart | one per thread | 27,788 (25,847–28,582) | 26 | 🔶 **1.00×** | 0.35× | PASS | 0 |
| 47 | MICM | standard, in-place Doolittle | one per thread | 27,880 (26,910–28,158) | 3 | 🔶 **1.00×** | 0.35× | PASS | 1e-13 |
| 48 | MICM | vector1, separate Mozart | one per thread | 27,910 (27,788–29,573) | 3 | 🔶 **1.00×** | 0.35× | PASS | 0 |
| 49 | MICM | vector1, separate Doolittle | one per thread | 28,316 (27,890–28,695) | 3 | 🔶 **0.98×** | 0.35× | PASS | 1e-13 |
| 50 | MICM | standard, separate Mozart | one per thread | 28,361 (28,036–29,072) | 3 | 🔶 **0.98×** | 0.34× | PASS | 0 |
| 51 | MICM | standard, separate Doolittle | one per thread | 28,487 (27,647–29,135) | 3 | 🔶 **0.98×** | 0.34× | PASS | 1e-13 |
| 52 | **GPT-6-Luna Max** | the agent's C port (`luna-max`) | one per thread | 31,201 (30,484–31,756) | 10 | 🔶 **0.89×** | 0.31× | PASS | 3e-13 |
| 53 | MICM | vector128, in-place Mozart | 32 cells | 36,531 (35,340–37,389) | 3 | 🔶 **0.76×** | 0.27× | PASS | 0 |
| 54 | MICM | vector128, in-place Doolittle | 32 cells | 36,919 (36,657–38,242) | 3 | 🔶 **0.75×** | 0.26× | PASS | 1e-13 |
| 55 | MICM | vector128, separate Doolittle | 32 cells | 47,392 (46,862–47,914) | 3 | 🔶 **0.59×** | 0.21× | PASS | 1e-13 |
| 56 | MICM | vector128, separate Mozart | 32 cells | 49,635 (49,204–50,984) | 3 | 🔶 **0.56×** | 0.20× | PASS | 0 |
