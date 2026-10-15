"""Chapman fixtures and the legacy fixed-solver regression comparisons.

The independent Python oracle follows the old nine-species solver, including its
1e-4 relative tolerance. The C implementation uses the generic mechanism builder.
Original fixtures: Copyright (C) 2023-2026 University Corporation for Atmospheric Research.
SPDX-License-Identifier: Apache-2.0
"""
import math
import random
from api import *
from harness import case, equal, near
from linear_cases import matrix_model, flatten
from analytical_cases import stiff_model

NAMES = ['M','Ar','CO2','H2O','N2','O1D','O','O2','O3']
PHOTO = ['jO2','jO3a','jO3b']
POSITIONS = [(0,0),(6,0),(7,0),(8,0),(1,1),(2,2),(3,3),(4,4),(5,4),(6,4),(5,5),(6,5),
             (6,6),(7,6),(8,6),(5,7),(6,7),(7,7),(8,7),(5,8),(6,8),(7,8),(8,8)]
TEMPS,PRESSURES = [284.19,215.02,299.31],[101245,100789.2,101398]
PHOTOS = [[1e-4,1e-5,1e-6],[3.2e-4,7.3e-5,3.2e-6],[5.2e-4,8.2e-5,4.6e-6]]


def chapman(integration=False):
    m = Model(NAMES)
    for rs,ps,r in [(['O1D','N2'],[('O',1),('N2',1)],rate(A=2.15e-11,C=110)),
      (['O1D','O2'],[('O',1),('O2',1)],rate(A=3.3e-11,C=55)),(['O','O3'],[('O2',2)],rate(A=8e-12,C=-2060)),
      (['O','O2','M'],[('O3',1),('M',1)],rate(A=6e-34,**({'C':2.4} if integration else {'B':2.4}))),
      (['O2'],[('O',2)],rate(USER,label='jO2')),(['O3'],[('O1D',1),('O2',1)],rate(USER,label='jO3a')),
      (['O3'],[('O',1),('O2',1)],rate(USER,label='jO3b'))]: m.reaction(rs,ps,r)
    return m.finalize()


def fixed_rates(t,photo):
    return [2.15e-11*math.exp(110/t),3.3e-11*math.exp(55/t),8e-12*math.exp(-2060/t),6e-34*(t/300)**2.4]+list(photo)


def fixed_force(k,y):
    f = [0.0]*9
    f[5] = k[5]*y[8]-k[0]*y[4]*y[5]-k[1]*y[5]*y[7]
    f[6] = 2*k[4]*y[7]+k[6]*y[8]+k[0]*y[4]*y[5]+k[1]*y[5]*y[7]-k[2]*y[6]*y[8]-k[3]*y[0]*y[6]*y[7]
    f[7] = -k[4]*y[7]+k[5]*y[8]+k[6]*y[8]+2*k[2]*y[6]*y[8]-k[3]*y[0]*y[6]*y[7]
    f[8] = -k[5]*y[8]-k[6]*y[8]-k[2]*y[6]*y[8]+k[3]*y[0]*y[6]*y[7]
    return f


def fixed_jac(k,y):
    j = [[0.0]*9 for _ in range(9)]
    j[5][4],j[5][5],j[5][7],j[5][8] = -k[0]*y[5],-k[0]*y[4]-k[1]*y[7],-k[1]*y[5],k[5]
    j[6][0],j[6][4],j[6][5] = -k[3]*y[6]*y[7],k[0]*y[5],k[0]*y[4]+k[1]*y[7]
    j[6][6],j[6][7],j[6][8] = -k[2]*y[8]-k[3]*y[0]*y[7],2*k[4]+k[1]*y[5]-k[3]*y[0]*y[6],k[6]-k[2]*y[6]
    j[7][0],j[7][6],j[7][7],j[7][8] = -k[3]*y[6]*y[7],2*k[2]*y[8]-k[3]*y[0]*y[7],-k[4]-k[3]*y[0]*y[6],k[5]+k[6]+2*k[2]*y[6]
    j[8][0],j[8][6],j[8][7],j[8][8] = k[3]*y[6]*y[7],-k[2]*y[8]+k[3]*y[0]*y[7],k[3]*y[0]*y[6],-k[5]-k[6]-k[2]*y[6]
    return j


def dense_lu(a):
    a = [row[:] for row in a]
    for k in range(len(a)):
        for i in range(k+1,len(a)):
            a[i][k] /= a[k][k]
            for j in range(k+1,len(a)): a[i][j] -= a[i][k]*a[k][j]
    return a


def dense_solve(a,b):
    x = list(b)
    for i in range(len(x)):
        for j in range(i): x[i] -= a[i][j]*x[j]
    for i in reversed(range(len(x))):
        for j in range(i+1,len(x)): x[i] -= a[i][j]*x[j]
        x[i] /= a[i][i]
    return x


def fixed_solve(k,y,dt):
    y = list(y); h,time,steps = 1e-6,0,0; reject,last = False,False
    cs = [-1.0156171083877702091975600115545,4.0759956452537699824805835358067,9.2076794298330791242156818474003]
    weights = [1,6.169794704382824559255361568973,-.42772256543218573326238373806514]
    errors = [.5,-2.9079558716805469821718236208017,.22354069897811569627360909276199]
    while time < dt:
        assert steps <= 100
        h = min(h,dt-time); force = fixed_force(k,y)
        while True:
            j = fixed_jac(k,y); alpha = 1/(h*.43586652150845899941601945119356)
            a = dense_lu([[(alpha if i == c else 0)-j[i][c] for c in range(9)] for i in range(9)])
            stages = [dense_solve(a,force)]
            f = fixed_force(k,[y[i]+stages[0][i] for i in range(9)])
            stages.append(dense_solve(a,[f[i]+cs[0]/h*stages[0][i] for i in range(9)]))
            stages.append(dense_solve(a,[f[i]+cs[1]/h*stages[0][i]+cs[2]/h*stages[1][i] for i in range(9)]))
            new_y = y[:]; error = [0.0]*9
            for st in range(3):
                for i in range(9): new_y[i] += weights[st]*stages[st][i]; error[i] += errors[st]*stages[st][i]
            norm = max(1e-10,math.sqrt(sum((error[i]/(1e-3+1e-4*max(abs(y[i]),abs(new_y[i]))))**2 for i in range(9))/9))
            hnew = h*min(6,max(.2,.9/norm**(1/3))); steps += 1
            if norm < 1:
                time += h; y = new_y; hnew = min(hnew,dt)
                h = min(hnew,h) if reject else hnew; reject,last = False,False; break
            if last: hnew = h*.1
            last,reject,h = reject,True,hnew
    return y


@case(*['chapman_ode_solver:ChapmanMechanismHardCodedAndGeneral.'+n for n in ('DefaultConstructor','lin_solve','ReactionNames','PhotolysisNames','SpeciesNames','simple_force','smaller_force','factored_alpha_minus_jac','dforce_dy_times_vector','Solve','solve_10_times_larger','solve_with_random_number_densities')],
      'chapman_ode_solver:RegressionChapmanODESolver.solve_10_times_smaller')
def legacy_cases(identifier):
    name = identifier.split('.')[-1]
    if name in ('lin_solve','factored_alpha_minus_jac','dforce_dy_times_vector'):
        m = matrix_model(9,set(POSITIONS)); s = State(m)
        try:
            a = arr([0]*m.data.nnz)
            for i,j in POSITIONS: a[m.data.index[i*9+j]] = 1
            if name == 'lin_solve':
                b = arr([.5]*9); lib.chem_linear_solve(m.p,1,1,a,b)
                equal(list(b),[.5,.5,.5,.5,.5,-.5,-1,.5,0])
            elif name == 'dforce_dy_times_vector':
                for i in range(5): a[m.data.index[i*9+i]] = 0
                b = arr([0]*9); lib.chem_sparse_multiply(m.p,1,1,a,arr([.5]*9),b)
                for actual,expected in zip(list(b)+[0]*14,[0,0,0,0,0,2,3,2,2]+[0]*14): near(actual,expected,.01)
            else:
                for i,j in POSITIONS: s.data.jacobian[m.data.index[i*9+j]] = -1
                lib.chem_alpha_minus_jacobian(s.p,2); lib.chem_lu_mozart(m.p,1,1,s.data.jacobian)
                expected = [1,-1,-1,-1,1,1,1,1,-1,-1,1,-1,1,-1,-1,-1,-2,-1,3,-1,-2,-3,.125]
                for (i,j),value in zip(POSITIONS,expected):
                    actual = s.data.jacobian[m.data.index[i*9+j]]
                    near(1/actual if i == j else actual,value,.01)
        finally: s.close(); m.close()
        return
    m = chapman(); s = State(m)
    try:
        if name == 'DefaultConstructor': equal(lib.chem_solver(3).parameters.stages,3)
        elif name == 'ReactionNames': equal(m.data.nreactions,7)
        elif name == 'PhotolysisNames': equal(m.data.ncustom,3)
        elif name == 'SpeciesNames': equal([m.data.species[i].name.decode() for i in range(9)],NAMES)
        elif name.endswith('force'):
            small = name == 'smaller_force'; y = 5e-6 if small else 1; k = 3e-8 if small else 1
            for i in range(9): s.data.variables[i] = y
            for i in range(7): s.data.rates[i] = k
            f = arr([0]*9); lib.chem_forcing(s.p,s.data.variables,f)
            for i in range(5): equal(f[i],0)
            for actual,expected in zip(list(f)[5:],[1.49e-13,4.55e-13,1.5e-13,-3e-13] if small else [-1,3,2,-2]):
                near(actual,expected,.01 if small else 0)
        else:
            scale = 10 if name == 'solve_10_times_larger' else .1 if name == 'solve_10_times_smaller' else 1
            initial = [scale*x for x in [1,.392,.0169,0,32.9,0,0,8.84,0]]
            if name == 'solve_with_random_number_densities': initial = [random.Random(i).random() for i in range(9)]
            for n,x in zip(NAMES,initial): s.set(n,x)
            s.data.conditions[0] = Conditions(273.15,100000,2.7e19)
            for n,x in zip(PHOTO,PHOTOS[0]): s.set(n,x,custom=True)
            s.rates(); s.solve(1)
            if name != 'solve_with_random_number_densities':
                expected = [1,.392,.0169,0,32.9,1.8039e-41,.00176789,8.83912,4.5031e-36]
                for n,x in zip(NAMES,expected): near(s.get(n),scale*x,1e-4)
    finally: s.close(); m.close()


@case(*[f'regression_test_rosenbrock_{part}:RegressionRosenbrock.{prefix}{name}' for part,names in [('p_force',('RateConstants','Forcing')),('dforce_dy',('Jacobian',))] for name in names for prefix in ('','Vector')])
def process_regression(identifier):
    name = identifier.split('.')[-1].removeprefix('Vector'); rng = random.Random(739)
    m = chapman(); s = State(m,3)
    try:
        if name == 'RateConstants':
            for c in range(3):
                s.data.conditions[c] = Conditions(TEMPS[c],PRESSURES[c],0)
                for j,p in enumerate(PHOTO): s.data.custom[j*3+c] = PHOTOS[c][j]
            s.rates()
            for c in range(3):
                for j,k in enumerate(fixed_rates(TEMPS[c],PHOTOS[c])): equal(s.data.rates[j*3+c],k)
        else:
            for i in range(27): s.data.variables[i] = rng.lognormvariate(-2,4 if name == 'Jacobian' else 2)
            for i in range(21): s.data.rates[i] = rng.lognormvariate(-2,4 if name == 'Jacobian' else 2)
            f = arr([0]*27)
            if name == 'Forcing': lib.chem_forcing(s.p,s.data.variables,f)
            else: lib.chem_jacobian(s.p,s.data.variables,s.data.jacobian)
            for c in range(3):
                y,k = [s.data.variables[j*3+c] for j in range(9)],[s.data.rates[j*3+c] for j in range(7)]
                if name == 'Forcing':
                    for j,expected in enumerate(fixed_force(k,y)):
                        actual = f[j*3+c]; near(actual,expected,(abs(actual)+abs(expected))*1e-8+1e-12)
                else:
                    expected = fixed_jac(k,y)
                    for i,j in POSITIONS:
                        actual = s.data.jacobian[m.data.index[i*9+j]*3+c]
                        near(actual,-expected[i][j],(abs(actual)+abs(expected[i][j]))*1e-8+1e-12)
    finally: s.close(); m.close()


@case(*['regression_test_solve:RegressionRosenbrock.'+n for n in ('TwoStageSolve','ThreeStageSolve','FourStageSolve','FourStageDASolve','SixStageDASolve','VectorSolve')])
def solve_regression(identifier):
    name = identifier.split('.')[-1]; method = {'TwoStageSolve':2,'ThreeStageSolve':3,'FourStageSolve':4,'FourStageDASolve':5,'SixStageDASolve':6,'VectorSolve':3}[name]
    m = chapman(); s = State(m,3); rng = random.Random(1)
    try:
        initial = [[rng.lognormvariate(-2,2) for _ in range(9)] for _ in range(3)]
        abs_tol = max(max(row) for row in initial)*1e-12
        for c in range(3):
            s.data.conditions[c] = Conditions(TEMPS[c],PRESSURES[c],0)
            for j in range(9): s.data.variables[j*3+c] = initial[c][j]
            for j in range(3): s.data.custom[j*3+c] = PHOTOS[c][j]
        s.rates(); s.solve(500,method)
        for c in range(3):
            reference = fixed_solve(fixed_rates(TEMPS[c],PHOTOS[c]),initial[c],500)
            for j,expected in enumerate(reference):
                actual = s.get(NAMES[j],c)
                near(actual,expected,(abs(actual)+abs(expected))*(1e-2 if method == 2 else 1e-4)+abs_tol)
    finally: s.close(); m.close()


@case('chapman_integration:ChapmanIntegration.CanBuildChapmanSystem')
def chapman_integration(_):
    m = chapman(True); s = State(m)
    try:
        names = ['O','O1D','O2','O3','M','Ar','N2','H2O','CO2']
        for n,x in zip(names,[.1,.1,.1,.2,.2,.2,.3,.3,.3]): s.set(n,x)
        for n,x in zip(PHOTO,[.1,.2,.3]): s.set(n,x,custom=True)
        s.data.conditions[0] = Conditions(2,3,0)
        oxygen = lambda: sum(s.get(n)*w for n,w in [('O',1),('O1D',1),('O2',2),('O3',3)])
        initial = oxygen(); spectators = {n:s.get(n) for n in names[4:]}
        for _ in range(100):
            s.rates(); s.solve(30)
            for n in names: assert math.isfinite(s.get(n)) and s.get(n) >= 0,(n,s.get(n))
            near(oxygen(),initial,abs(initial)*1e-8)
            for n,x in spectators.items(): near(s.get(n),x,abs(x)*1e-8)
        assert s.get('O2') < .5*.1 and s.get('O') > .1 and s.get('O') > .5*initial
    finally: s.close(); m.close()


@case('integrated_reaction_rates:ChapmanIntegration.CanBuildChapmanSystem')
def integrated_rates(_):
    m = Model(['a','b','c','irr_1','irr_2'])
    m.reaction(['a'],[('b',1),('irr_1',1)],rate(A=2.15e-11,C=110))
    m.reaction(['b'],[('c',1),('irr_2',1)],rate(USER,label='r2')); m.finalize(); s = State(m)
    try:
        s.set('a',1); s.set('r2',1,custom=True); s.data.conditions[0] = Conditions(273,1000,0)
        for t in range(100):
            if t > 50: s.set('r2',0,custom=True)
            s.rates(); before = [s.get('irr_1'),s.get('irr_2')]; s.solve(30)
            assert all(s.get(n) >= x for n,x in zip(['irr_1','irr_2'],before))
    finally: s.close(); m.close()


@case('terminator:RosenbrockSolver.Terminator','terminator:RosenbrockSolver.VectorTerminator')
def terminator(_):
    for cells in range(1,5):
        m = Model(['Cl2','Cl']); m.reaction(['Cl2'],[('Cl',2)],rate(USER,label='toy_k1'))
        m.reaction(['Cl','Cl'],[('Cl2',1)],rate(A=1)); m.finalize(); s = State(m,cells); rng = random.Random(1)
        try:
            s.tolerances(1e-8,[1e-20]*2); p = lib.chem_solver(3); p.parameters.max_steps = 100000
            initial = [[rng.lognormvariate(-2,2)*1e-6,rng.lognormvariate(-2,2)*1e-10] for c in range(cells)]
            for step in range(math.floor(2*math.pi/.3)):
                ks = []
                for c in range(cells):
                    lat = math.pi/180*(c*(90/cells)); lat0,lon0 = math.pi*20/180,math.pi*300/180
                    ks.append(max(0,math.sin(lat)*math.sin(lat0)+math.cos(lat)*math.cos(lat0)*math.cos(step*.3-lon0)))
                    for j in range(2): s.data.variables[j*cells+c] = initial[c][j]
                    s.data.conditions[c] = Conditions(298,101300,42)
                s.set('toy_k1',ks,custom=True); s.rates(); s.solve(30,solver=p)
                for c,(cl2,cl) in enumerate(initial):
                    r = ks[c]/4; det = math.sqrt(r*r+2*r*(cl+2*cl2)); e = math.exp(-4*det*30)
                    ell = (1-e)/det/30 if det*30 > 1e-16 else 4
                    f = -ell*(cl-det+r)*(cl+det+r)/(1+e+30*ell*(cl+r))
                    for n,ref in [('Cl',cl+30*f),('Cl2',cl2-30*f/2)]: near(s.get(n,c),ref,abs(ref)*1e-6+1e-14)
        finally: s.close(); m.close()


@case('openmp:OpenMP.OneSolverManyStates')
def shared_solver(_):
    m,s,names = stiff_model('Robertson'); states = [s]+[State(m) for _ in range(7)]
    try:
        for state in states:
            state.tolerances(); state.set('A',1)
            for n,k in [('r1',.04),('r2',3e7),('r3',1e4)]: state.set(n,k,custom=True)
            state.data.conditions[0] = Conditions(272.5,101253.3,1e6); state.rates()
        pointers = (SP*8)(*(state.p for state in states)); stats = (Stats*8)(); p = lib.chem_solver(3)
        for _ in range(10): equal(lib.chem_advance_many(C.cast(lib.chem_step_rosenbrock,P),C.byref(p),200,1,pointers,8,8,stats),0)
        for state in states[1:]: equal([state.get(n) for n in names],[s.get(n) for n in names])
    finally:
        for state in states: state.close()
        m.close()
