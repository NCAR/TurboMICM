# Declared operating points

These five configurations are the complete evaluation list. All use the full TS1 mechanism, double precision, fresh Jacobians, maximum cell RMS error control, and `atol=1e-18` mol m-3. Use all cores by omitting `--threads`; the measured runs used ten workers.

Append the scenario path to each command below. All other parameters retain the defaults in `DESIGN.md`.

## fast

```sh
bin/ts1_scenario --method rodas3 --option rtol=0.01 --option group=16 --option schedule=dynamic SCENARIO
```

## balanced

```sh
bin/ts1_scenario --method rodas3 --option rtol=0.005 --option group=16 --option schedule=dynamic SCENARIO
```

## accurate

```sh
bin/ts1_scenario --method rodas3 --option rtol=0.0005 --option group=16 --option schedule=dynamic SCENARIO
```

## recommended — recommended for an atmospheric model

```sh
bin/ts1_scenario --method rodas3 --option rtol=0.0003 --option group=16 --option schedule=dynamic SCENARIO
```

## precise

```sh
bin/ts1_scenario --method rodas3 --option rtol=0.0001 --option group=16 --option schedule=dynamic SCENARIO
```

| Point | rtol | Lowest disclosed overall digits | Lowest disclosed species digits |
| --- | --- | --- | --- |
| fast | 0.01 | 2.041 | 1.064 |
| balanced | 0.005 | 2.290 | 1.386 |
| accurate | 0.0005 | 3.273 | 2.341 |
| recommended | 0.0003 | 3.571 | 2.554 |
| precise | 0.0001 | 4.250 | 3.069 |

The recommended point supplies margin over the three-digit overall / two-digit species target, including the additional 48-hour validation. `accurate` is the cheapest declared point meeting that target on every disclosed scenario. `fast` is cheapest at level two but has little accuracy margin; `balanced` provides a more conservative level-two alternative. `precise` supports tighter sensitivity work.

These are measured expectations on the disclosed cases, not promises about unknown inputs. The complete raw results, additional scenarios, projection measurements, and scope of the 10,000-cell reference checks are in `RESULTS.md`, `ACCURACY.md`, and `scans/README.md`. The inherited default solver and research alternatives remain callable but are not additional declared evaluation points.
