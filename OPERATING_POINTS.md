# Declared operating points

Five configurations are declared. Use all available cores by default, or add
`--threads N`; final local measurements use ten threads. All use the complete
mechanism, double precision, group size 32, `atol=1e-18`, fresh Jacobians,
proportional control, maximum-cell RMS error, warm step-size starts and no
intermediate-stage clipping. Rate caching is enabled with four-step blocks,
dynamic scheduling and the OpenMP default chunk. These are runtime parameters.

**Recommended: fourth-order Euler extrapolation, rtol=0.002.** It gives measured
margin above three overall digits and two digits for every scored species,
including the additional 48-hour cases. The accurate point costs less but has
less margin; the strict point provides additional margin at greater cost. The
fast point only narrowly clears the two-digit species criterion; balanced
provides more measured margin within that cheaper tier.

| Point | Method | rtol | Expected level | Minimum overall/species digits, disclosed 1k |
| --- | --- | --- | --- | --- |
| Fast | rodas3 | 0.01 | 2 | 2.041/1.035 |
| Balanced | rodas3 | 0.003 | 2 | 2.645/1.646 |
| Accurate | extrap order 4 | 0.004 | 3 | 3.201/2.390 |
| Recommended | extrap order 4 | 0.002 | 3 | 3.448/2.550 |
| Strict | extrap order 4 | 0.001 | 3 | 4.064/2.919 |

## Fast

```sh
bin/ts1_scenario --method rodas3 \
  --option rtol=0.01 --option atol=1e-18 --option group-size=32 \
  --option rate-cache=1 --option step-block=4 --option schedule=dynamic \
  --option chunk=0 --option jacobian-age=1 --option controller=0 \
  --option norm=1 --option warm-start=1 --option reuse-window=0 \
  --option clip-stages=0 \
  SCENARIO
```

## Balanced

```sh
bin/ts1_scenario --method rodas3 \
  --option rtol=0.003 --option atol=1e-18 --option group-size=32 \
  --option rate-cache=1 --option step-block=4 --option schedule=dynamic \
  --option chunk=0 --option jacobian-age=1 --option controller=0 \
  --option norm=1 --option warm-start=1 --option reuse-window=0 \
  --option clip-stages=0 \
  SCENARIO
```

## Accurate

```sh
bin/ts1_scenario --method extrap \
  --option rtol=0.004 --option atol=1e-18 --option group-size=32 \
  --option order=4 --option rate-cache=1 --option step-block=4 \
  --option schedule=dynamic --option chunk=0 --option jacobian-age=1 \
  --option controller=0 --option norm=1 --option warm-start=1 \
  --option reuse-window=0 --option clip-stages=0 \
  SCENARIO
```

## Recommended

```sh
bin/ts1_scenario --method extrap \
  --option rtol=0.002 --option atol=1e-18 --option group-size=32 \
  --option order=4 --option rate-cache=1 --option step-block=4 \
  --option schedule=dynamic --option chunk=0 --option jacobian-age=1 \
  --option controller=0 --option norm=1 --option warm-start=1 \
  --option reuse-window=0 --option clip-stages=0 \
  SCENARIO
```

## Strict

```sh
bin/ts1_scenario --method extrap \
  --option rtol=0.001 --option atol=1e-18 --option group-size=32 \
  --option order=4 --option rate-cache=1 --option step-block=4 \
  --option schedule=dynamic --option chunk=0 --option jacobian-age=1 \
  --option controller=0 --option norm=1 --option warm-start=1 \
  --option reuse-window=0 --option clip-stages=0 \
  SCENARIO
```

`SCENARIO` is the path to a scenario file. Machine-readable configurations are
in `scans/operating_points.json`; other method defaults are in DESIGN.md.
The numerical levels are expectations supported by measured scenarios, not
proofs for unseen forcing. Every declared point passes its intended level on
all disclosed and added stress cases with both compilers. The 10,000-cell gate
covers all cells, while independent accuracy covers only its first 1,000 cells.
See RESULTS.md for per-scenario accuracy, five-run timing summaries, work–precision
plots and the matched MICM comparison, and ACCURACY.md for scientific limitations.
