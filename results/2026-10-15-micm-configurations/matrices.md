Median ms at 10 threads (speed relative to MICM's default), one cell group per thread:

| Matrix layout | in-place Mozart (default) | in-place Doolittle | separate Mozart | separate Doolittle |
|---|---:|---:|---:|---:|
| standard | 27,751 (1.00×) | 27,880 (1.00×) | 28,361 (0.98×) | 28,487 (0.97×) |
| vector1 | 27,357 (1.01×) | 27,358 (1.01×) | 27,910 (0.99×) | 28,316 (0.98×) |
| vector2 | 22,380 (1.24×) | 22,115 (1.25×) | 22,908 (1.21×) | 23,531 (1.18×) |
| vector4 | 19,354 (1.43×) | 19,076 (1.45×) | 20,311 (1.37×) | 20,241 (1.37×) |
| vector8 | 18,192 (1.53×) | 18,014 (1.54×) | 19,123 (1.45×) | 19,027 (1.46×) |
| vector128 | 19,627 (1.41×) | 19,672 (1.41×) | 22,085 (1.26×) | 20,914 (1.33×) |

Median ms at 10 threads (speed relative to MICM's default), 32-cell groups:

| Matrix layout | in-place Mozart (default) | in-place Doolittle | separate Mozart | separate Doolittle |
|---|---:|---:|---:|---:|
| standard | 23,480 (1.18×) | 22,759 (1.22×) | 24,694 (1.12×) | 24,231 (1.15×) |
| vector1 | 22,242 (1.25×) | 22,443 (1.24×) | 24,114 (1.15×) | 24,145 (1.15×) |
| vector2 | 15,311 (1.81×) | 15,361 (1.81×) | 16,640 (1.67×) | 17,056 (1.63×) |
| vector4 | 11,362 (2.44×) | 11,309 (2.45×) | 12,815 (2.17×) | 12,873 (2.16×) |
| vector8 | **9,752 (2.85×)** | 9,819 (2.83×) | 10,890 (2.55×) | 10,990 (2.53×) |
| vector128 | 36,531 (0.76×) | 36,919 (0.75×) | 49,635 (0.56×) | 47,392 (0.59×) |
