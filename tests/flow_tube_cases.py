"""Original 3600-step wall-loss regression and its recorded concentration rows."""
import csv
import math
from api import *
from harness import case, equal, near


@case('regression:RosenbrockRegression.FlowTube','regression:BackwardEulerRegression.FlowTube')
def flow_tube(identifier):
    backward = 'BackwardEuler' in identifier; names = ['SOA1','SOA2','a-pinene','O3']; m = Model(names)
    m.reaction(['SOA1'],[],rate(USER,label='r1')); m.reaction(['SOA2'],[],rate(USER,label='r2'))
    m.reaction(['a-pinene','O3'],[('SOA1',.18),('SOA2',.09)],rate(A=8.8e-17*(1e-6*6.02214076e23)))
    m.finalize(); s = State(m)
    try:
        path = ROOT/'tests/flow-tube'/('backward_euler_wall_loss.csv' if backward else 'rosenbrock_wall_loss.csv')
        with path.open() as f: rows = list(csv.reader(f))
        equal(rows[0],['time']+names); equal(len(rows)-1,3601); equal(len(rows[1])-1,4)
        s.set('a-pinene',8e-8); s.set('O3',.00002); s.set('r1',.01,custom=True); s.set('r2',.05,custom=True)
        s.data.conditions[0] = Conditions(298.15,101325,101325/((1.380649e-23*6.02214076e23)*298.15)); s.rates()
        for i in range(3601):
            if i < 10:
                for n,expected in zip(names,map(float,rows[i+1][1:])): near(s.get(n),expected,1e-10)
            if i < 3600: s.advance(1,method=0 if backward else 2)
        for n in names: assert math.isfinite(s.get(n))
    finally: s.close(); m.close()
