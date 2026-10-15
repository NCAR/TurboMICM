"""Original Chapman and flow-tube regression quantities and tolerances."""
from functools import partial
from pathlib import Path
import csv
from chemffi import *
from test_chapman import mechanism,reactions,NAMES,PHOTOS,TEMPERATURES,PRESSURES,condition,POSITIONS
from test_linear import bind
from reference_random import Minstd,MT19937,lognormals
import reference_chapman as reference


def source_permutation(declaration=NAMES,integration=False):
    positions={name:i for i,name in enumerate(declaration)};entries=set()
    for rr,pp,r in reactions(integration):
        for col in rr:entries.update((positions[NAMES[row]],positions[NAMES[col]]) for row in rr+[i for i,v in pp])
    coords=(MatrixIndex*len(entries))(*[MatrixIndex(*ij) for ij in sorted(entries)])
    perm=(Z*9)();assert lib.chem_diagonal_markowitz(9,len(entries),coords,perm)==0
    return [NAMES.index(declaration[old]) for old in perm]


def storage_rows(numbers,columns,group):
    length=3*columns if not group else ((3+group-1)//group)*columns*group
    raw=[next(numbers) for _ in range(length)]
    return [[raw[c*columns+i if not group else (c//group*columns+i)*group+c%group] for i in range(columns)] for c in range(3)]


def regression(kind,group=0):
    perm=source_permutation();numbers=lognormals(4 if kind=='Jacobian' else 2,MT19937(739))
    if kind=='RateConstants':
        for start in range(0,3,group or 3):
            count=min(group or 3,3-start)
            with mechanism(count) as m:
                condition(m,TEMPERATURES[start:],PRESSURES[start:],PHOTOS[start:])
                for cell in range(count):
                    want=reference.rates(TEMPERATURES[start+cell],PHOTOS[start+cell])
                    for i,v in enumerate(want):
                        got=lib.chem_rate_constants(m.state)[i*count+cell]
                        assert got==v,(start+cell,i,got,v)
        return
    raw=storage_rows(numbers,9,group)
    if kind=='Jacobian':raw=storage_rows(numbers,9,group) # source intentionally fills concentrations twice
    variables=[[0.]*9 for _ in range(3)]
    for cell in range(3):
        for index,old in enumerate(perm):variables[cell][old]=raw[cell][index]
    constants=storage_rows(numbers,7,group)
    # Grouping is a runtime tuning choice. Preserve reference input/padding
    # generation and apply every actual cell with independent scratch per group.
    for start in range(0,3,group or 3):
        count=min(group or 3,3-start)
        with mechanism(count) as m:
            m.set_y(variables[start:start+count])
            for c in range(count):
                for i in range(7):lib.chem_rate_constants(m.state)[i*count+c]=constants[start+c][i]
            actual=m.jacobian() if kind=='Jacobian' else m.force()
            for c in range(count):
                expected=reference.jacobian(constants[start+c],variables[start+c]) if kind=='Jacobian' else reference.force(constants[start+c],variables[start+c])
                pairs=[(actual[c][row][col],-expected[row][col]) for row,col in POSITIONS] if kind=='Jacobian' else list(zip(actual[c],expected))
                for got,want in pairs:assert abs(got-want)<=(abs(got)+abs(want))*1e-8+1e-12,(kind,start+c,got,want)


def regression_solve(method,group=0):
    numbers=lognormals(2,Minstd());raw=[[next(numbers) for _ in range(9)] for _ in range(3)]
    perm=source_permutation();variables=[[0.]*9 for _ in range(3)]
    for c in range(3):
        for i,old in enumerate(perm):variables[c][old]=raw[c][i]
    atol=max(v for row in raw for v in row)*1e-12
    expected=[reference.solve(row,reference.rates(TEMPERATURES[c],PHOTOS[c])) for c,row in enumerate(variables)]
    for start in range(0,3,group or 3):
        count=min(group or 3,3-start)
        with mechanism(count) as m:
            m.set_y(variables[start:start+count]);condition(m,TEMPERATURES[start:],PRESSURES[start:],PHOTOS[start:])
            m.solve(500.,method)
            for c,row in enumerate(m.rows()):
                for i,(got,want) in enumerate(zip(row,expected[start+c])):
                    tolerance=(abs(got)+abs(want))*(1e-2 if method==2 else 1e-4)+atol
                    assert abs(got-want)<=tolerance,(method,group,start+c,i,got,want,tolerance)


def flowtube(method):
    r1,r2=rate(9),rate(9);r1.parameter=0;r2.parameter=1
    r3=rate(0);r3.p[0]=8.8e-17*(1e-6*6.02214076e23)
    filename='rosenbrock' if method else 'backward_euler'
    with (Path(__file__).parent/'data'/(filename+'_wall_loss.csv')).open() as f:
        rows=list(csv.reader(f));assert rows[0]==['time','SOA1','SOA2','a-pinene','O3']
        expected=[[float(v) for v in row] for row in rows[1:]]
    with Mechanism(['SOA1','SOA2','a-pinene','O3'],[([0],[],r1),([1],[],r2),([2,3],[(0,.18),(1,.09)],r3)],parameters=['r1','r2']) as m:
        m.set_y([[0,0,8e-8,.00002]]);m.custom[0]=.01;m.custom[1]=.05
        lib.chem_temperatures(m.state)[0]=298.15;lib.chem_pressures(m.state)[0]=101325
        assert lib.chem_calculate_air_density(m.state)==0 and lib.chem_update_rates(m.state)==0
        assert len(expected)==3601
        for t in range(3601):
            if t<10:
                for got,want in zip(m.rows()[0],expected[t][1:]):assert abs(got-want)<=1e-10,(t,got,want)
            if t==3600:break
            elapsed=0
            while elapsed<1:
                if method:status,stats=m.solve(1-elapsed,method)
                else:
                    p=lib.chem_backward_euler_parameters();stats=Stats()
                    status=lib.chem_solve_backward_euler(m.state,C.byref(p),1-elapsed,C.byref(stats));m.rows()
                assert stats.final_time>0,(t,status)
                elapsed+=stats.final_time


def vector(fn,*args):
    for group in [1,2,3,4]:fn(*args,group=group)


CASES={}
for name in ['RateConstants','Forcing','Jacobian']:
    host='regression_test_rosenbrock_dforce_dy' if name=='Jacobian' else 'regression_test_rosenbrock_p_force'
    CASES[host+':RegressionRosenbrock.'+name]=partial(regression,name)
    CASES[host+':RegressionRosenbrock.Vector'+name]=partial(vector,regression,name)
for name,method in [('TwoStageSolve',2),('ThreeStageSolve',3),('FourStageSolve',4),('FourStageDASolve',-4),('SixStageDASolve',6)]:
    CASES['regression_test_solve:RegressionRosenbrock.'+name]=partial(regression_solve,method)
CASES['regression_test_solve:RegressionRosenbrock.VectorSolve']=partial(vector,regression_solve,3)
CASES['regression:RosenbrockRegression.FlowTube']=partial(flowtube,2)
CASES['regression:BackwardEulerRegression.FlowTube']=partial(flowtube,0)
