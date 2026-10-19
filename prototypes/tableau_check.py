"""Check the library's transformed Rosenbrock coefficients against the
standard form, independently.

    python prototypes/tableau_check.py [--h 60] [--methods ros34pw2,...]

For each method, one step of size h from the mini scenario's first cell is
taken two ways: (1) in the standard form of Hairer and Wanner with a dense
Jacobian and numpy's solver, (I - h gamma J) k_i = h f(y0 + sum alpha_ij k_j)
+ h J sum gamma_ij k_j, from the standard-form tables written here from the
literature; (2) in MICM's transformed form with the library's a, c, m, e
(prototypes/methods.py, the library's sparse LU). The two solutions and error
estimates must agree to roundoff. Methods without a standard-form table here
are skipped.
"""
import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from methods import Options, RowSolver, Counts  # noqa: E402
from scenario import read_csv, read_scenario  # noqa: E402
from tmclib import Library  # noqa: E402

# standard-form tables: alpha (strictly lower), gamma_ij (strictly lower), gamma, b, bhat
STANDARD = {
    "ros34pw2": dict(  # Rang and Angermann (2005), as printed by Blom et al. (2016), Table 7
        gamma=4.3586652150845900e-01,
        alpha=[[0, 0, 0, 0], [8.7173304301691801e-01, 0, 0, 0], [8.4457060015369423e-01, -1.1299064236484185e-01, 0, 0], [0, 0, 1, 0]],
        gamma_ij=[[0, 0, 0, 0], [-8.7173304301691801e-01, 0, 0, 0], [-9.0338057013044082e-01, 5.4180672388095326e-02, 0, 0],
                  [2.4212380706095346e-01, -1.2232505839045147e+00, 5.4526025533510214e-01, 0]],
        b=[2.4212380706095346e-01, -1.2232505839045147e+00, 1.5452602553351020e+00, 4.3586652150845900e-01],
        bhat=[3.7810903145819369e-01, -9.6042292212423178e-02, 5.0000000000000000e-01, 2.1793326075422950e-01]),
    "ros2": dict(  # Verwer et al. (1999): alpha_21 = 1, gamma_21 = -2 gamma, b = (1/2, 1/2), embedded first order (1, 0)
        gamma=1.0 + 1.0 / np.sqrt(2.0), alpha=[[0, 0], [1, 0]], gamma_ij=[[0, 0], [-2.0 * (1.0 + 1.0 / np.sqrt(2.0)), 0]],
        b=[0.5, 0.5], bhat=[1.0, 0.0]),
}


def standard_step(state, y0, h, tab):
    """One step in the standard form with a dense Jacobian; returns (y1, error)."""
    s = len(tab["b"])
    n = y0.shape[1]
    J = state.dense_jacobian(0, y0)
    M = np.eye(n) - h * tab["gamma"] * J
    k = []
    for i in range(s):
        arg = y0[0] + sum(tab["alpha"][i][j] * k[j] for j in range(i))
        f = state.forcing(arg[None, :])[0]
        rhs = h * f + h * J @ sum((tab["gamma_ij"][i][j] * k[j] for j in range(i)), np.zeros(n))
        k.append(np.linalg.solve(M, rhs))
    y1 = y0[0] + sum(tab["b"][i] * k[i] for i in range(s))
    err = sum((tab["b"][i] - tab["bhat"][i]) * k[i] for i in range(s))
    return y1, err


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--scenario", type=Path, default=Path(__file__).resolve().parent / "mini" / "scenario.txt")
    parser.add_argument("--h", type=float, default=60.0)
    parser.add_argument("--methods", default=",".join(STANDARD))
    parser.add_argument("--tolerance", type=float, default=1e-6)
    args = parser.parse_args()
    lib = Library()
    sc = read_scenario(args.scenario)
    header, rows = read_csv(sc["cells"])
    state = lib.state(1)
    T = P = None
    for name, value in zip(header, rows[0]):
        if name == "temperature":
            T = value
        elif name == "pressure":
            P = value
        elif name in lib.variable_index:
            state.set_variable(0, lib.variable_index[name], value)
        elif name in lib.custom_index:
            state.set_custom(0, lib.custom_index[name], value)
    state.set_conditions(0, T, P)
    state.update_rate_constants()
    y0 = state.variables()
    failures = 0
    for name in args.methods.split(","):
        if name not in STANDARD:
            print(f"SKIP {name}: no standard-form table here")
            continue
        y_std, e_std = standard_step(state, y0, args.h, STANDARD[name])
        solver = RowSolver(lib, lib.tableau(name), Options())
        state.factor(y0, 1.0 / (args.h * solver.t.gamma))
        F0 = state.forcing(y0)
        y_lib, e_lib = solver.stages(state, y0, F0, args.h, Counts())
        scale = np.maximum(np.abs(y_std), 1e-30)
        dy = np.max(np.abs(y_lib[0] - y_std) / scale)
        de = np.max(np.abs(e_lib[0] - e_std) / np.maximum(np.abs(e_std), np.max(np.abs(e_std)) * 1e-6))
        ok = dy <= args.tolerance and de <= 1e-3
        failures += 0 if ok else 1
        print(f"{'PASS' if ok else 'FAIL'} {name}: h = {args.h:g} s; solution relative difference {dy:.2e}, error-estimate relative difference {de:.2e}; "
              f"library m = {np.array2string(solver.t.m[:solver.t.stages], precision=6)}, e = {np.array2string(solver.t.e[:solver.t.stages], precision=6)}")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
