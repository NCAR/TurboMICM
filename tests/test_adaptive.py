"""Small-case agreement with independent Python stage/Newton implementations."""
import csv
import math
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'prototypes'))
from integrators import rosenbrock, bdf, COUNTERS
from api import C, Model, State, lib, rate
from adaptive_api import AdaptiveStats, counts


def mechanism():
    m = Model(['A', 'B', 'C'])
    m.reaction(['A'], [('B', 1.)], rate(A=2.))
    m.reaction(['B'], [('A', 1.)], rate(A=.3))
    m.reaction(['A', 'B'], [('C', 1.)], rate(A=5.))
    m.reaction(['C'], [('A', 1.), ('B', 1.)], rate(A=.2))
    m.data.default_temperature = 300.; m.data.default_pressure = 101325.
    return m.finalize()


def forcing(y):
    a, b, c = y
    result = [-2*a, 2*a, 0.]
    result[1] -= .3*b; result[0] += .3*b
    result[0] -= 5*a*b; result[1] -= 5*a*b; result[2] += 5*a*b
    result[2] -= .2*c; result[0] += .2*c; result[1] += .2*c
    return result


def jacobian(y):
    a, b, _ = y
    return [[-2-5*b, .3-5*a, .2], [2-5*b, -.3-5*a, .2], [5*b, 5*a, -.2]]


def main():
    records = []
    m = mechanism()
    methods = [(method, 5) for method in [2, 3, 4, 5, 6]] + [(10, order) for order in range(1, 6)]
    for method, order in methods:
        for variant in ['fixed', 'adaptive', 'reuse', 'maxnorm', 'filtered']:
            if method == 10 and variant == 'filtered':
                continue
            p = lib.chem_adaptive_parameters(method)
            p.max_order = order
            p.h_start = .01 if variant == 'fixed' else .08
            p.h_max = .01 if variant == 'fixed' else .2
            p.jacobian_age = 8 if variant == 'reuse' else 1
            p.reuse_ratio = .2
            if variant == 'reuse':
                p.h_start = p.h_max = .01; p.step_hold = 1.2
            if variant == 'maxnorm':
                p.norm = 2
            p.controller = int(method != 10 and variant == 'adaptive')
            if variant == 'filtered':
                p.controller = 2
            settings = {key: getattr(p, key) for key, _ in p._fields_}
            state = State(m); state.set('A', 1.); state.set('B', .1); state.rates()
            state.tolerances(relative=1e-4, absolute=1e-9)
            w = lib.chem_workspace_create(state.p); assert w
            stats = AdaptiveStats()
            code = lib.chem_adaptive_solve(C.byref(p), .75, state.p, w, C.byref(stats))
            assert code == 0, (method, variant, code)
            if method == 10:
                expected, operations = bdf(forcing, jacobian, [1., .1, 0.], .75, settings, 1e-4, 1e-9)
            else:
                expected, operations = rosenbrock(forcing, jacobian, [1., .1, 0.], .75,
                                                 lib.chem_solver(method).parameters, settings, 1e-4, 1e-9)
            actual = [state.get(name) for name in ['A', 'B', 'C']]
            diff = max(abs(a-b) for a, b in zip(actual, expected))
            assert diff < 3e-12, (method, variant, diff, actual, expected)
            assert counts(stats) == operations, (method, variant, counts(stats), operations)
            assert abs(actual[0] + actual[1] + 2*actual[2] - 1.1) < 3e-12
            records.append(dict(method=method, max_order=order, variant=variant, max_absolute_difference=diff, **operations))
            lib.chem_workspace_free(w); state.close()
    # Multiple lanes and a partial group, with a stiff conservative transfer.
    stiff = Model(['X', 'Y']); stiff.reaction(['X'], [('Y', 1.)], rate(A=1e6)); stiff.finalize()
    stiff.data.default_temperature = 300.; stiff.data.default_pressure = 101325.
    for method in [2, 3, 4, 5, 6, 10]:
        state = State(stiff, 3); state.set('X', [1., .1, 5.]); state.rates()
        state.tolerances(relative=1e-5, absolute=1e-12)
        p = lib.chem_adaptive_parameters(method); p.jacobian_age = 10 if method == 10 else 1
        w = lib.chem_workspace_create(state.p); stats = AdaptiveStats()
        assert lib.chem_adaptive_solve(C.byref(p), 1., state.p, w, C.byref(stats)) == 0
        for cell, initial in enumerate([1., .1, 5.]):
            assert math.isfinite(state.get('X', cell)) and state.get('X', cell) >= 0.
            assert state.get('X', cell) <= 2e-5*initial + 1e-12
            assert abs(state.get('Y', cell)-initial) <= 2e-5*initial + 1e-12
        # Invalid settings fail, leaving state unchanged.
        p.h_start = -1.; before = state.get('Y')
        assert lib.chem_adaptive_solve(C.byref(p), 1., state.p, w, C.byref(stats)) == 1
        assert state.get('Y') == before
        lib.chem_workspace_free(w); state.close()
    stiff.close(); m.close()
    dest = ROOT / 'build' / 'prototype-verification.csv'
    with dest.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['method', 'max_order', 'variant', 'max_absolute_difference'] + COUNTERS, lineterminator='\n')
        writer.writeheader(); writer.writerows(records)
    print(f'prototype agreement: {len(records)} configurations, matching counts; stiff/positivity/API checks PASS')


if __name__ == '__main__':
    main()
