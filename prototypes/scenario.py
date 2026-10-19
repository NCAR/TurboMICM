"""Run a prototype method on a scenario (or its first cells) and score it.

    python prototypes/scenario.py SCENARIO --method rodas3 --cells 4 [--option KEY=VALUE ...]
        [--reference REF_CSV] [--csv-row FILE --label NAME]

Reads the scenario as bin/ts1_scenario does (../scenarios/README.md), runs the
cells as one group, applies f(t) to the photolysis parameters at each step
(recomputing only the rate constants that read them), records the snapshots,
and prints the operation counts and, with a reference (reference_snapshots.csv
or reference_final.csv; the scenario's own by default), the accuracy in
significant digits over the cells run, with score.py's definitions. With
--csv-row it appends one row to a scan file. The options are the C driver's.
"""
import argparse
import csv
import math
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from methods import Options, Counts, make_solver  # noqa: E402
from tmclib import Library  # noqa: E402


def read_csv(path):
    with open(path, newline="") as file:
        rows = list(csv.reader(file))
    return rows[0], [[float(x) for x in row] for row in rows[1:] if row]


def read_scenario(path):
    path = Path(path)
    values = {}
    for line in path.read_text().splitlines():
        line = line.split("#", 1)[0].strip()
        if not line:
            continue
        key, _, value = line.partition("=")
        values[key.strip()] = value.strip()
    return {
        "cells": path.parent / values["cells"],
        "photolysis": path.parent / values["photolysis"],
        "dt": float(values["dt"]),
        "steps": int(values["steps"]),
        "snapshot_every": int(values.get("snapshot_every", 0)),
    }


class Photolysis:
    def __init__(self, path):
        header, rows = read_csv(path)
        assert header == ["time_s", "factor"], "photolysis file: time_s,factor"
        self.time = np.array([r[0] for r in rows])
        self.factor = np.array([r[1] for r in rows])

    def __call__(self, t):
        if t <= self.time[0]:
            return float(self.factor[0])
        if t >= self.time[-1]:
            return float(self.factor[-1])
        return float(np.interp(t, self.time, self.factor))


def parse_options(items, defaults=None):
    o = defaults or Options()
    kinds = {"rtol": float, "atol": float, "warm_start": lambda v: v in ("1", "yes", "true"), "h_start": float,
             "h_start_factor": float, "h_max": float, "norm": str, "controller": str,
             "clamp": lambda v: v in ("1", "yes", "true"), "quantize": lambda v: v in ("1", "yes", "true"),
             "reuse_age": int, "reuse_change": float, "freeze_band": float, "max_steps": int, "order": int, "sequence": str,
             "safety": float, "fac_min": float, "fac_max": float, "h211_b": float, "h211_k": float,
             "filter_error": lambda v: v in ("1", "yes", "true")}
    for item in items or []:
        key, _, value = item.partition("=")
        if key not in kinds:
            sys.exit(f"unknown option {key}; known: {' '.join(kinds)}")
        setattr(o, key, kinds[key](value))
    return o


def digits(mean_square):
    if mean_square <= 0.0:
        return 16.0
    return min(16.0, -math.log10(math.sqrt(mean_square)))


def score(result_rows, reference_rows, species, ref_species, threshold=1e-14):
    """score.py's accuracy over the given rows: (overall digits, {species: digits}, gate)."""
    col = {name: i for i, name in enumerate(species)}
    sum_sq = count = 0
    per = {}
    nonfinite = below = 0
    for key, ref in reference_rows.items():
        res = result_rows[key]
        for j, name in enumerate(ref_species):
            y_ref = ref[j]
            y = res[col[name]]
            if not math.isfinite(y):
                nonfinite += 1
            elif y < -1e-20:
                below += 1
            if y_ref < threshold:
                continue
            err = abs(y - y_ref) / y_ref if math.isfinite(y) else math.inf
            sum_sq += err * err
            count += 1
            s = per.setdefault(name, [0.0, 0])
            s[0] += err * err
            s[1] += 1
    overall = digits(sum_sq / count) if count else float("nan")
    per_species = {name: digits(t / n) for name, (t, n) in per.items()}
    return overall, per_species, nonfinite == 0 and below == 0


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("scenario")
    parser.add_argument("--method", default="rodas3")
    parser.add_argument("--cells", type=int, default=4, help="the first CELLS cells (default 4; 0: all)")
    parser.add_argument("--option", action="append", default=[])
    parser.add_argument("--reference", help="reference_snapshots.csv or reference_final.csv (default: the scenario's)")
    parser.add_argument("--generic", action="store_true")
    parser.add_argument("--csv-row", type=Path, help="append a row to this scan file")
    parser.add_argument("--label", default="")
    parser.add_argument("--write", type=Path, help="write the snapshots (or final state) to this CSV, like the C driver")
    args = parser.parse_args()

    sc = read_scenario(args.scenario)
    header, rows = read_csv(sc["cells"])
    if args.cells > 0:
        rows = rows[: args.cells]
    n_cells = len(rows)
    photolysis = Photolysis(sc["photolysis"])
    lib = Library(generic=args.generic)
    o = parse_options(args.option)
    if args.method == "micm3":
        o.warm_start = False
    state = lib.state(n_cells)
    kinds = []
    for name in header:
        if name == "temperature":
            kinds.append(("T", 0))
        elif name == "pressure":
            kinds.append(("P", 0))
        elif name in lib.variable_index:
            kinds.append(("v", lib.variable_index[name]))
        elif name in lib.custom_index:
            kinds.append(("p", lib.custom_index[name]))
        else:
            sys.exit(f"column {name} is not a TS1 species or rate parameter")
    for c, row in enumerate(rows):
        T = P = None
        for (kind, index), value in zip(kinds, row):
            if kind == "T":
                T = value
            elif kind == "P":
                P = value
            elif kind == "v":
                state.set_variable(c, index, value)
            else:
                state.set_custom(c, index, value)
        state.set_conditions(c, T, P)
    state.set_tolerances(o.rtol, o.atol)
    state.update_rate_constants()
    base = np.array([[state.get_custom(c, k) for k in lib.photolysis] for c in range(n_cells)])

    snapshot_steps = list(range(sc["snapshot_every"], sc["steps"] + 1, sc["snapshot_every"])) if sc["snapshot_every"] else []
    if not snapshot_steps or snapshot_steps[-1] != sc["steps"]:
        snapshot_steps.append(sc["steps"])
    solver = make_solver(lib, args.method, o)
    y = state.variables()
    counts = Counts()
    snapshots = {}
    last_factor = None
    status = "converged"
    t0 = time.perf_counter()
    for step in range(sc["steps"]):
        factor = photolysis(step * sc["dt"])
        if factor != last_factor:
            for c in range(n_cells):
                for jj, k in enumerate(lib.photolysis):
                    state.set_custom(c, k, base[c, jj] * factor)
            state.update_custom_rate_constants()
            last_factor = factor
        result = solver.solve(state, y, sc["dt"])
        counts.add(result.counts)
        y = np.maximum(result.y, 0.0)  # the library's clamp after every solve
        if result.status != "converged":
            status = f"step {step}: {result.status}"
            break
        if step + 1 in snapshot_steps:
            snapshots[step + 1] = y.copy()
    seconds = time.perf_counter() - t0

    overall = worst = float("nan")
    worst_name = ""
    gate = status == "converged"
    reference = Path(args.reference) if args.reference else None
    if reference is None:
        for name in ("reference_snapshots.csv", "reference_final.csv"):
            if (Path(args.scenario).parent / name).exists():
                reference = Path(args.scenario).parent / name
                break
    if reference and status == "converged":
        ref_header, ref_rows = read_csv(reference)
        keyed = ref_header[:2] == ["step", "cell"]
        start = 4 if keyed else 2
        ref_species = ref_header[start:]
        ref_table, res_table = {}, {}
        if keyed:
            for row in ref_rows:
                step, cell = int(row[0]), int(row[1])
                if cell < n_cells and step in snapshots:
                    ref_table[(step, cell)] = row[start:]
                    res_table[(step, cell)] = snapshots[step][cell]
        else:
            for cell, row in enumerate(ref_rows[:n_cells]):
                ref_table[cell] = row[start:]
                res_table[cell] = y[cell]
        overall, per_species, finite = score(res_table, ref_table, lib.variable_names, ref_species)
        gate = gate and finite
        worst_name, worst = min(per_species.items(), key=lambda kv: kv[1])

    if args.write:
        with open(args.write, "w", newline="") as file:
            w = csv.writer(file)
            if sc["snapshot_every"]:
                w.writerow(["step", "cell", "temperature", "pressure", *lib.variable_names])
                for step in snapshot_steps:
                    for c in range(n_cells):
                        w.writerow([step, c, f"{rows[c][header.index('temperature')]:.17g}", f"{rows[c][header.index('pressure')]:.17g}", *[f"{v:.17g}" for v in snapshots[step][c]]])
            else:
                w.writerow(["temperature", "pressure", *lib.variable_names])
                for c in range(n_cells):
                    w.writerow([f"{rows[c][header.index('temperature')]:.17g}", f"{rows[c][header.index('pressure')]:.17g}", *[f"{v:.17g}" for v in y[c]]])

    fields = {
        "label": args.label, "scenario": Path(args.scenario).parent.name, "method": args.method, "cells": n_cells,
        "options": " ".join(args.option), "status": status, "gate": "PASS" if gate else "FAIL",
        "digits": f"{overall:.3f}", "worst_species": worst_name, "worst_digits": f"{worst:.3f}",
        "steps": counts.steps, "accepted": counts.accepted, "rejected": counts.rejected,
        "function_calls": counts.function_calls, "jacobian_updates": counts.jacobian_updates,
        "decompositions": counts.decompositions, "solves": counts.solves, "matrix_reuses": counts.matrix_reuses,
        "python_s": f"{seconds:.2f}",
    }
    print(" ".join(f"{k}={v}" for k, v in fields.items()))
    if args.csv_row:
        new = not args.csv_row.exists()
        with open(args.csv_row, "a", newline="") as file:
            w = csv.DictWriter(file, fieldnames=list(fields))
            if new:
                w.writeheader()
            w.writerow(fields)
    return 0 if status == "converged" else 1


if __name__ == "__main__":
    sys.exit(main())
