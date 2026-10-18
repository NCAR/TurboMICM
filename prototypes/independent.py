"""Independent dense prototype of the adaptive ROS/Euler-extrapolation formulas.

Only the retained Rosenbrock coefficient tables are read from the C probe.
Chemistry, dense elimination, integration, and control are evaluated in Python.
No Python timings are reported. Run with the checkout's .conda/bin/python.
"""
import csv
import json
import math
import subprocess
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def factor(a):
    a = a.copy()
    for j in range(3):
        inv = 1.0 / a[j, j]
        for i in range(j + 1, 3):
            a[i, j] *= inv
        for k in range(j + 1, 3):
            for i in range(j + 1, 3):
                a[i, k] -= a[i, j] * a[j, k]
    return a


def solve(a, b):
    b = b.copy()
    for i in range(3):
        for j in range(i):
            b[i] -= a[i, j] * b[j]
    for i in reversed(range(3)):
        for j in range(i + 1, 3):
            b[i] -= a[i, j] * b[j]
        b[i] /= a[i, i]
    return b


def integrate(table, method, max_age, controller, window, stiff):
    y = np.array([.8, .15, .05])
    counts = np.zeros(13, dtype=int)
    adaptive_order = method == 10
    if adaptive_order:
        method = 7

    def f(y):
        counts[3] += 1
        a, b, c = y
        # Preserve reaction accumulation order; explicit mechanism independent
        # of the C reaction table/evaluator and its sparse Jacobian jobs.
        out = np.array([-2 * stiff * a, 2 * stiff * a, 0.])
        out[1] -= .7 * b; out[2] += .7 * b
        out[2] -= .1 * c; out[0] += .1 * c
        rate = .3 * a * b
        out[0] -= rate; out[1] -= rate; out[1] += 2 * rate
        return out

    def jac(y):
        counts[4] += 1
        a, b, _ = y
        return np.array([[-2 * stiff - .3 * b, -.3 * a, .1],
                         [2 * stiff + .3 * b, -.7 + .3 * a, 0], [0, .7, -.1]])

    def decompose(h, gamma=1):
        counts[5] += 1
        return factor(np.eye(3) / (h * gamma) - j)

    def linear(a, b):
        counts[6] += 1
        return solve(a, b)

    t, h, age, old_error, factor_h = 0., .1, max_age, 1., 0.
    f0 = None
    rejected = False
    order = table["order"] if method < 5 else method - (3 if method < 8 else 6)
    nominal_order = order
    if method >= 8:
        if order == 2:
            g = 1-1/math.sqrt(2)
            a = np.array([[g, 0], [1-g, g]])
            e = np.array([-g, g])
        else:
            g = .43586652150845899941601945119356
            a = np.array([[g, 0, 0], [(1-g)/2, g, 0],
                          [-1.5*g*g+4*g-.25, 1.5*g*g-5*g+1.25, g]])
            low1 = (.5-g)/((1-g)/2)
            e = a[-1] - [1-low1, low1, 0]
    while t < 3:
        assert counts[0] < 100000
        h = min(h, 3 - t)
        assert h > 0 and t + h != t
        fresh = age == 0 or age >= max_age
        if age >= max_age:
            j = jac(y); age = 0; factor_h = 0
        else:
            counts[7] += 1
        if method < 5:
            if f0 is None:
                f0 = f(y)
            if h != factor_h:
                lu = decompose(h, table["gamma"][0]); factor_h = h
            else:
                counts[8] += 1
            k = []; cached = None
            for stage in range(table["stages"]):
                base = stage * (stage - 1) // 2
                if stage == 0:
                    rhs = f0.copy()
                elif table["new_function"][stage]:
                    stage_y = y.copy()
                    for z in range(stage):
                        stage_y += table["a"][base + z] * k[z]
                    rhs = f(stage_y)
                else:
                    rhs = cached.copy()
                cached = rhs.copy()
                for z in range(stage):
                    rhs += table["c"][base + z] / h * k[z]
                k.append(linear(lu, rhs))
            candidate = y.copy(); error = np.zeros(3)
            for stage in range(table["stages"]):
                candidate += table["m"][stage] * k[stage]
                error += table["e"][stage] * k[stage]
        elif method < 8:
            if f0 is None:
                f0 = f(y)
            tableau = []
            for row in range(nominal_order):
                if row:
                    low = tableau[0].copy()
                lu = decompose(h / (row + 1)); v = y.copy()
                for substep in range(row + 1):
                    v += linear(lu, f(v) if substep else f0)
                tableau.append(v)
                for z in reversed(range(1, row + 1)):
                    tableau[z-1] = tableau[z] + (tableau[z] - tableau[z-1]) / ((row+1)/z - 1)
                if row:
                    error = tableau[0] - low
                    normalized = np.maximum(abs(error), np.maximum(0, -tableau[0])) / (
                        1e-10 + 1e-5*np.maximum(abs(y), abs(tableau[0])))
                    if adaptive_order and math.sqrt(sum(v*v for v in normalized)/3) <= 1:
                        break
            counts[10+row-1] += 1
            candidate = tableau[0]; error = candidate - low
            order = row+1
        else:
            if h != factor_h:
                lu = decompose(h, g); factor_h = h
            else:
                counts[8] += 1
            k = []; failed = False
            for stage in range(order):
                base = y.copy()
                for z in range(stage):
                    base += h*a[stage, z]*k[z]
                v = base.copy(); previous = math.inf
                for iteration in range(8):
                    delta = linear(lu, f(v) + (1/(h*g))*(base-v))
                    counts[9] += 1
                    v += delta
                    normalized = delta / (1e-10 + 1e-5*np.maximum(abs(y), abs(v)))
                    norm = math.sqrt(sum(z*z for z in normalized)/3)
                    if norm <= .03:
                        break
                    if not math.isfinite(norm) or (iteration >= 2 and norm > .9*previous):
                        failed = True; break
                    previous = norm
                else:
                    failed = True
                if failed:
                    break
                k.append(f(v))
            candidate = v; error = np.zeros(3)
            if failed:
                error[:] = math.inf
            else:
                for stage in range(order):
                    error += h*e[stage]*k[stage]
        counts[0] += 1
        normalized = np.maximum(np.abs(error), np.maximum(0, -candidate)) / (
            1e-10 + 1e-5 * np.maximum(np.abs(y), np.abs(candidate)))
        norm = math.sqrt(sum(v*v for v in normalized)/3)
        proposal = .9 * max(norm, 1e-12) ** (-1/order)
        if controller and not rejected and norm <= 1:
            proposal = .9 * max(norm, 1e-12)**(-.7/order) * old_error**(.4/order)
        if norm <= 1:
            y = np.maximum(candidate, 0)
            counts[1] += 1; age += 1; t += h; f0 = None
            old_error = max(norm, 1e-12)
            proposal = min(1 if rejected else 5, max(.2, proposal))
            if abs(proposal - 1) <= window:
                proposal = 1
            h = min(3, h * proposal); rejected = False
        else:
            counts[2] += 1
            if not fresh:
                age = max_age
            else:
                h *= min(.5, max(.1, proposal))
            factor_h = 0; rejected = True
    return y, counts


def main():
    rows = []
    for method in range(11):
        for age, controller, window, stiff in [(1, 0, 0, 1), (1, 1, 0, 100),
                                               (3, 0, .15, 1), (3, 1, .15, 100)]:
            result = subprocess.run([str(ROOT / "bin" / "adaptive_probe"), str(method), str(age),
                                     str(controller), str(window), str(stiff)],
                                    capture_output=True, text=True, check=True)
            c = json.loads(result.stdout)
            y, counts = integrate(c, method, age, controller, window, stiff)
            error = max(abs(y - c["y"]))
            assert error < 128 * np.finfo(float).eps, (method, error, y, c["y"])
            assert list(counts) == c["counts"], (method, counts, c["counts"])
            assert abs(sum(y) - 1) < 1e-12
            rows.append(dict(method=method, age=age, controller=controller, window=window,
                             stiffness=stiff, max_absolute_difference=error,
                             counts=json.dumps(c["counts"]), passed=True))
    with (ROOT / "results" / "prototype_agreement.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    print(f"{len(rows)} prototype comparisons pass; maximum difference {max(r['max_absolute_difference'] for r in rows):.3g}")


if __name__ == "__main__":
    main()
