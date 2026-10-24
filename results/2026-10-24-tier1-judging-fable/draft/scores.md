## Scores

Totals out of 30; extensibility is twice the mean of the five exercise ratings
(E1 to E5 in brackets).

| Rank | Entrant | Readability | Maintainability | Extensibility (E1–E5) | Total of 30 |
|---:|---|---:|---:|---|---:|
| 1 | [astra-high](Entrant-astra-high) | 8 | 8 | 8.8 (5, 4, 5, 4, 4) | 24.8 |
| 2 | [sol-xhigh](Entrant-sol-xhigh) | 7 | 8 | 8.4 (5, 4, 5, 3, 4) | 23.4 |
| 3 | [astra-xhigh](Entrant-astra-xhigh) | 7 | 7 | 9.2 (5, 5, 5, 4, 4) | 23.2 |
| 4 | [sol-high](Entrant-sol-high) | 7 | 8 | 8.0 (5, 4, 4, 3, 4) | 23.0 |
| 5 | [astra-max](Entrant-astra-max) | 6 | 7 | 8.4 (5, 4, 5, 3, 4) | 21.4 |
| 6 | [sol-max](Entrant-sol-max) | 6 | 6 | 8.0 (5, 4, 4, 3, 4) | 20.0 |
| 7 | [luna-max](Entrant-luna-max) | 5 | 4 | 7.2 (4, 3, 5, 3, 3) | 16.2 |
| 8 | [luna-xhigh](Entrant-luna-xhigh) | 4 | 5 | 6.4 (5, 2, 3, 3, 3) | 15.4 |
| 9 | [luna-high](Entrant-luna-high) | 3 | 2 | 4.4 (2, 2, 2, 3, 2) | 9.4 |

Three groups. The six Sol and Astra entries are one generic design each, with
the same five parameter sets, a sparse LU with a precomputed plan, a hook for
the linear solver and, in six of them, a hook for the controller; they differ
in how much of that is readable and tested, and in whether generated kernels
duplicate the ordering logic. Ranks 2 to 4 are within 0.4 of a point and the
judge would not defend their order against a different weighting of the same
facts; the gap between rank 1 and rank 2 and the gap between rank 6 and rank 7
are real. The three Luna entries port MICM's cases less faithfully, accept a
NaN step, and convert configuration units with fixed exponents; luna-high
prints its test results from a list and generates wrong reactions for any
mechanism whose species order its generator changes.
