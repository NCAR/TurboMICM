"""Compare every generated solver with the generic C library at 1e-10 relative.

Run by make: Python only drives validation. Every C source is compiled by make.
Failures in rates, kernels, sparse solves, or integration make this exit nonzero.
"""
import csv
import ctypes as C
import math
import random
import sys
from pathlib import Path

from api import (ROOT, lib, MP, SP, Z, D, DP, I, Stats, Model, State, Conditions,
                 arr, bind)

specialized = C.CDLL(str(ROOT / 'bin' / ('libspecialized.dylib' if sys.platform == 'darwin' else 'libspecialized.so')))


def function(name, result, *arguments):
    fn = getattr(specialized, name)
    fn.restype, fn.argtypes = result, arguments
    return fn


def compare(actual, expected, count, what, rtol=1e-10):
    for i in range(count):
        a, b = actual[i], expected[i]
        assert math.isfinite(a) and math.isfinite(b), (what, i, a, b)
        assert abs(a-b) <= 1e-30+rtol*abs(b), (what, i, a, b, abs(a-b)/(abs(b)+1e-300))


def initialize(state, name):
    s, m = state.data, state.model.data
    if name.startswith('ts1'):
        source = ROOT.parent / 'reference' / 'ts1_cells_100_seed1.csv'
        with source.open() as file:
            rows = list(csv.DictReader(file))
        for c in range(s.cells):
            for key, value in rows[c].items():
                value = float(value)
                if key in ('temperature', 'pressure'):
                    setattr(s.conditions[c],key,value)
                else:
                    v = lib.chem_variable_id(state.model.p,key.encode())
                    if v != Z(-1).value:
                        s.variables[v*s.stride+c] = value
                    else:
                        u = lib.chem_custom_id(state.model.p,key.encode())
                        assert u != Z(-1).value
                        s.custom[u*s.stride+c] = value
    else:
        for c in range(s.cells):
            s.conditions[c] = Conditions(270+c%19*2,85000+c%23*900,0)
            for v in range(m.nvariables):
                initial = m.species[m.variable_species[v]].initial
                s.variables[v*s.stride+c] = (initial if initial else 1e-8)*(1+c*.007)
            for u in range(m.ncustom):
                label = m.custom_names[u].decode()
                value = 1e-7 if 'radius' in label else 2.5e6 if 'number concentration' in label else .0002
                s.custom[u*s.stride+c] = value*(1+c*.003)
    for c in range(s.cells):
        lib.chem_ideal_air_density(C.byref(s.conditions[c]))
    state.tolerances(1e-6,1e-15)


def check_mechanism(name, cells):
    loader_name = 'chem_load_'+name
    if hasattr(lib,loader_name):
        generic_loader = bind(loader_name,MP,I)
    else:
        generic_loader = function(loader_name,MP,I)
    gm = Model(pointer=generic_loader(0))
    sm = Model(pointer=function('chem_load_'+name+'_specialized',MP,I)(1))
    gs, ss = State(gm,cells), State(sm,cells)
    try:
        initialize(gs,name); initialize(ss,name)
        g, s, m = gs.data, ss.data, gm.data
        assert (m.nvariables,m.nnz,m.nreactions,m.ncustom) == (sm.data.nvariables,sm.data.nnz,sm.data.nreactions,sm.data.ncustom)
        for i in range(m.nvariables*m.nvariables):
            assert m.index[i] == sm.data.index[i]
        for i in range(m.nvariables):
            assert m.variable_species[i] == sm.data.variable_species[i]
        rate_fn = function('chem_rates_'+name+'_specialized',None,SP)
        solve_fn = function('chem_solve_'+name+'_specialized',I,D,SP,C.POINTER(Stats))
        factor_fn = function('chem_factor_'+name+'_specialized',None,MP,Z,Z,DP)
        linear_fn = function('chem_linear_'+name+'_specialized',None,MP,Z,Z,DP,DP)
        # Exercise linear algebra with a well-conditioned matrix independent of chemistry.
        rng = random.Random(17)
        matrix = [rng.uniform(-.1,.1) for _ in range(m.nnz*cells)]
        for v in range(m.nvariables):
            for c in range(cells): matrix[m.diagonal[v]*cells+c] += 100.0
        a, b = arr(matrix), arr(matrix)
        lib.chem_lu_mozart(gm.p,cells,cells,a); factor_fn(sm.p,cells,cells,b)
        compare(b,a,len(matrix),'LU')
        rhs = [rng.uniform(-1,1) for _ in range(m.nvariables*cells)]
        x, y = arr(rhs), arr(rhs)
        lib.chem_linear_solve(gm.p,cells,cells,a,x); linear_fn(sm.p,cells,cells,b,y)
        compare(y,x,len(rhs),'triangular solve')
        # Conditions and all custom parameters change between calls. Density is
        # deliberately independent of ideal gas conditions after the first step.
        for step, dt in enumerate((30.0, 60.0, 120.0) if name.startswith('ts1') else (.1,.2,1.0)):
            if step:
                for state in (gs,ss):
                    st = state.data
                    for c in range(cells):
                        st.conditions[c].temperature += 3.5
                        st.conditions[c].pressure *= .99
                        st.conditions[c].air_density *= 1.01
                    for u in range(m.ncustom):
                        for c in range(cells): st.custom[u*cells+c] *= .8 if step == 1 else 1.1
            gs.rates(); rate_fn(ss.p)
            compare(s.rates,g.rates,m.nreactions*cells,'rate constants')
            for kernel, count in (('chem_forcing',m.nvariables),('chem_jacobian',m.nnz)):
                a, b = arr([0]*(count*cells)), arr([0]*(count*cells))
                getattr(lib,kernel)(gs.p,g.variables,a); getattr(lib,kernel)(ss.p,s.variables,b)
                compare(b,a,count*cells,kernel)
            expected = gs.solve(dt)
            actual = Stats(); status = solve_fn(dt,ss.p,C.byref(actual))
            assert status == 0, (name,cells,step,status,actual.final_time)
            assert actual.final_time == expected.final_time == dt
            compare(s.variables,g.variables,m.nvariables*cells,'integrated state')
            assert all(s.variables[i] >= 0 for i in range(m.nvariables*cells))
        print(f'{name} cells={cells} PASS rates, chemistry, LU, solves, adaptive integration',flush=True)
    finally:
        gs.close(); ss.close(); gm.close(); sm.close()


def main():
    for name in ('ts1','example_no_user','example_user','example_robertson','validation','edges','ts1_variant'):
        for cells in (1,7,32,33):
            check_mechanism(name,cells)


if __name__ == '__main__':
    main()
