"""Generated-versus-generic chemistry, sparse layout and trajectory checks."""
import csv
import ctypes as C
import math
from pathlib import Path
import random
import sys

from chemffi import Data, Stats, Rosenbrock, D, P, DP, lib

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'tools'))
from configure import normalize, reorder
from specialize import symbolic

CONFIGS = {
    'ts1': 'ts1.json',
    'no_user': 'rate_constants_no_user_defined/config.yaml',
    'user': 'rate_constants_user_defined/config.json',
    'robertson': 'robertson/config.json',
    'parameterized': 'parameterized/config.json',
    'checks': 'specialization_checks.json',
    'inert': 'inert.json',
}
owner = C.CDLL(str(ROOT / 'build/specialized_tests.so'))


def close(a, b, context, rtol=1e-10):
    assert math.isfinite(a) and math.isfinite(b), (context, a, b)
    assert abs(a - b) <= rtol * abs(b) + 1e-30, (context, a, b)


def check(prefix, config, cells):
    data = Data.in_dll(owner, prefix + '_data')
    rate = getattr(owner, prefix + '_update_rates')
    rate.argtypes, rate.restype = [P], C.c_int
    solve = getattr(owner, prefix + '_solve')
    solve.argtypes, solve.restype = [P, C.POINTER(Rosenbrock), D, C.POINTER(Stats)], C.c_int
    mech, generic, special = P(), P(), P()
    assert lib.chem_mechanism_create(C.byref(data), C.byref(mech)) == 0
    try:
        # Check EVERY coordinate, including fill and structural zeros.
        normalized = normalize(ROOT / 'configs' / config)
        if [s['name'] for s in reorder(normalized)['species']] == [data.species_names[i].decode() for i in range(data.species_count)]:
            normalized = reorder(normalized)
        plan = symbolic(normalized)
        assert len(plan) == lib.chem_jacobian_size(mech)
        for row in range(data.species_count):
            for col in range(data.species_count):
                assert lib.chem_jacobian_index(mech, row, col) == plan.get((row, col), C.c_size_t(-1).value)
        assert lib.chem_state_create(mech, cells, C.byref(generic)) == 0
        assert lib.chem_state_create(mech, cells, C.byref(special)) == 0
        rng = random.Random(1000 + cells)
        source = list(csv.DictReader((ROOT.parent / 'reference/ts1_cells_100_seed1.csv').open())) if prefix == 'ts1' else None
        for c in range(cells):
            t = float(source[c]['temperature']) if source else rng.uniform(270, 310)
            pressure = float(source[c]['pressure']) if source else rng.uniform(80000, 102000)
            for s in (generic, special):
                lib.chem_temperatures(s)[c] = t
                lib.chem_pressures(s)[c] = pressure
            for i in range(data.species_count):
                name = data.species_names[i].decode()
                value = float(source[c].get(name, data.initial_concentrations[i])) if source else (1.0 if name == 'A' else 0.0)
                for s in (generic, special):
                    lib.chem_concentrations(s)[i * cells + c] = value
            for i in range(data.parameter_count):
                name = data.parameter_names[i].decode()
                if source:
                    value = float(source[c].get(name, data.initial_parameters[i]))
                elif name.endswith('.effective radius [m]'):
                    value = 1e-7
                elif name.endswith('.particle number concentration [# m-3]'):
                    value = 2.5e6
                elif prefix == 'robertson':
                    value = [.04, 3e7, 1e4][i]
                else:
                    value = 1e-8
                for s in (generic, special):
                    lib.chem_parameters(s)[i * cells + c] = value
        for s in (generic, special):
            assert lib.chem_calculate_air_density(s) == 0
        params = lib.chem_rosenbrock_parameters(3)
        params.max_steps = 10000
        # Tight tolerances exercise the error controller and rejected steps as
        # well as the loose default benchmark trajectory.
        for tight in (False, True):
            if tight:
                for s in (generic, special):
                    for i in range(data.species_count): lib.chem_absolute_tolerances(s)[i] = 1e-15
            for step in range(6):
                # Every custom parameter remains runtime data. Vary photolysis,
                # surface, emission and kinetic parameters identically.
                for i in range(data.parameter_count * cells):
                    v = lib.chem_parameters(generic)[i] * (0.5 if step % 2 else 1.5)
                    lib.chem_parameters(generic)[i] = lib.chem_parameters(special)[i] = v
                for s in (generic, special):
                    for c in range(cells):
                        lib.chem_temperatures(s)[c] += .1
                        lib.chem_pressures(s)[c] *= .999
                        lib.chem_air_densities(s)[c] *= 1.001
                assert lib.chem_update_rates(generic) == rate(special) == 0
                for i in range(data.reaction_count * cells):
                    close(lib.chem_rate_constants(special)[i], lib.chem_rate_constants(generic)[i], (prefix, 'rate', i), 1e-13)
                dt = 30 if prefix == 'ts1' else .1
                a, b = Stats(), Stats()
                assert lib.chem_solve(generic, C.byref(params), dt, C.byref(a)) == 0, (prefix, 'generic', a.final_time)
                assert solve(special, C.byref(params), dt, C.byref(b)) == 0, (prefix, 'special', b.final_time)
                for key, _ in Stats._fields_:
                    assert getattr(a, key) == getattr(b, key), (prefix, step, key, getattr(a, key), getattr(b, key))
                for i in range(data.species_count * cells):
                    close(lib.chem_concentrations(special)[i], lib.chem_concentrations(generic)[i], (prefix, cells, tight, step, i))
        # Public argument checking must not dereference a missing state.
        assert solve(None, C.byref(params), 1, C.byref(Stats())) == 1
        assert rate(None) == 1
        print(f'PASS {prefix}: {cells} cells, sparse plan, dynamic rates, default/tight trajectories')
    finally:
        lib.chem_state_free(generic)
        lib.chem_state_free(special)
        lib.chem_mechanism_free(mech)


if __name__ == '__main__':
    for prefix, config in CONFIGS.items():
        for cells in (1, 7, 32, 33):
            check(prefix, config, cells)
    print('All specialized mechanisms agree with the generic library to 1e-10 relative.')
