"""Rosenbrock building blocks, parameter overrides, and DAE initialization ports."""
import math
from api import *
from harness import case, equal, near, real_equal


def simple_system(constrained=True):
    m = Model(['A','B','C'] if constrained else ['A','B'])
    m.reaction(['A'],[('B',1)],rate(A=.5))
    if constrained: m.constraint('C',[('B',1),('C',1)],2,equilibrium_reactants=1,name='B_C_eq')
    return m.finalize()


def five_species():
    m = Model(['foo','bar','baz','quz','quuz'])
    for rs,ps,r in [(['foo','baz'],[('bar',1),('quuz',2.4)],rate(A=2e-11,C=110)),
                    (['bar'],[('foo',1),('quz',1.4)],rate(A=1e-6)),(['quz'],[],rate(A=3.5e-6))]: m.reaction(rs,ps,r)
    return m.finalize()


@case(*['rosenbrock:RosenbrockSolver.'+name for name in ('StandardAlphaMinusJacobian','VectorAlphaMinusJacobian')])
def alpha_minus_jacobian(_):
    for cells in (1,2,3,4):
        m = five_species(); s = State(m,cells)
        try:
            data = {(0,0):-12.2,(0,2):-42.3,(1,0):-.43,(1,1):-23.4,(2,0):-4.74,(2,2):-6.91,(3,1):-59.1,(3,3):-83.4,(4,0):-78.5,(4,2):-53.6,(4,4):-1}
            expected = []
            for c in range(cells):
                values = dict(data); values[0,1] = -24.3*(c+2); values[1,2] = -83.4/(c+3); expected.append(values)
                for (i,j),v in values.items(): s.data.jacobian[m.data.index[i*5+j]*cells+c] = v
            lib.chem_alpha_minus_jacobian(s.p,42.042)
            for c,values in enumerate(expected):
                for (i,j),v in values.items(): real_equal(s.data.jacobian[m.data.index[i*5+j]*cells+c],v+(42.042 if i == j else 0))
        finally: s.close(); m.close()


@case('rosenbrock:RosenbrockSolver.CanSetTolerances')
def setting_tolerances(_):
    m = Model(); m.species('foo',absolute_tolerance=1e-7); m.species('bar',absolute_tolerance=1e-8)
    m.reaction(['foo'],[('bar',1)],rate(A=2e-11,C=110)); m.finalize()
    try:
        for cells in range(1,11):
            s = State(m,cells)
            try: equal(list(s.data.absolute_tolerance[:2]),[1e-7,1e-8])
            finally: s.close()
    finally: m.close()


@case(*['rosenbrock:RosenbrockSolver.'+name for name in ('StandardNormalizedError','VectorNormalizedError','StandardNormalizedErrorWithConstraints','VectorNormalizedErrorWithConstraints')])
def normalized_error(identifier):
    constrained = 'WithConstraints' in identifier
    for cells in (1,2,3,4,5):
        m = simple_system() if constrained else five_species(); s = State(m,cells)
        n = m.data.nvariables
        old,new,errors = arr([0]*(n*cells)),arr([0]*(n*cells)),arr([0]*(n*cells))
        try:
            total = 0
            for c in range(cells):
                for i in range(n):
                    at = i*cells+c
                    old[at],new[at],errors[at] = ((1+c+.1*i,.8+.5*c+.2*i,.01*(1+c+i)) if constrained else (7.7*c+i,-13.9/(i+1)-c,81.57/(c+7)/(i+3)))
                    scale = s.data.absolute_tolerance[i]+s.data.relative_tolerance*max(abs(old[at]),abs(new[at]))
                    total += errors[at]*errors[at]/(scale*scale)
            expected = max(math.sqrt(total/(cells*n)),1e-10)
            computed = lib.chem_normalized_error(s.p,old,new,errors)
            near(computed,expected,1e-12 if constrained else 1e-11*max(abs(expected),abs(computed)))
        finally: s.close(); m.close()


@case('rosenbrock:RosenbrockSolver.RejectedStepAlphaMatchesInPlaceSolver')
def rejected_step(_):
    results = []
    for doolittle in (True,False):
        m = five_species(); s = State(m)
        try:
            for name in ('foo','bar','baz','quz','quuz'): s.set(name,1e12)
            s.data.conditions[0].temperature = 298.15; s.rates()
            solver = lib.chem_solver(3); solver.parameters.h_start = 1000
            if doolittle: solver.factor = C.cast(lib.chem_lu_doolittle,P)
            stats = s.solve(1000,solver=solver); assert stats.internal_steps-stats.accepted >= 2
            results.append((stats.internal_steps,stats.accepted,stats.rejected))
        finally: s.close(); m.close()
    equal(results[0],results[1])


@case(*['solver_builder:SolverBuilder.'+name for name in ('ThrowsMissingSystem','ThrowsUnusedSpecies','CanBuildBackwardEuler','CanBuildRosenbrock',
    'CanBuildBackwardEulerOverloadedSolverMethod','CanBuildRosenbrockOverloadedSolveMethod','FinalizeReceivesDenseMatrixPolicy','FinalizeSupportsSparseOnlySignature')])
def builder(identifier):
    name = identifier.split('.')[-1]
    if name == 'ThrowsMissingSystem':
        m = Model()
        try: equal(lib.chem_model_finalize(m.p),1)
        finally: m.close()
        return
    m = Model(['A','B','C']+(['D'] if name == 'ThrowsUnusedSpecies' else []))
    m.reaction(['A'],[('B',1)],rate(A=2.15e-11,C=110)); m.reaction(['B'],[('C',1)],rate(A=3.3e-11,C=55))
    if name == 'ThrowsUnusedSpecies':
        m.data.ignore_unused = 0
        try: equal(lib.chem_model_finalize(m.p),3)
        finally: m.close()
        return
    calls = []
    SetupFn = C.CFUNCTYPE(None,MP,P)
    callback = SetupFn(lambda model,ctx:calls.append((model.contents.nvariables,model.contents.nnz)))
    if name.startswith('Finalize'): m.data.setup = C.cast(callback,P)
    m.finalize()
    try:
        if name.startswith('Finalize'): equal(len(calls),1); equal(calls[0][0],3); assert calls[0][1] >= 3
        if name in ('CanBuildBackwardEuler','CanBuildRosenbrock'):
            for cells in (1,4):
                s = State(m,cells)
                try: equal(s.data.cells,cells); equal(s.model.data.nspecies,3)
                finally: s.close()
        elif 'BackwardEulerOverloaded' in name:
            s = State(m)
            try:
                p = lib.chem_euler_parameters(); stats = Stats()
                equal(lib.chem_backward_euler(C.byref(p),5,s.p,C.byref(stats)),0); equal(stats.final_time,5)
                for key in ('function_calls','jacobian_updates','internal_steps','solves'): equal(getattr(stats,key),2)
                p.small = 1; p.max_steps = 1
                equal(lib.chem_backward_euler(C.byref(p),5,s.p,C.byref(stats)),10); equal(stats.final_time,.03125)
                for key in ('function_calls','jacobian_updates','internal_steps','solves'): equal(getattr(stats,key),6)
            finally: s.close()
        elif 'RosenbrockOverloaded' in name:
            s = State(m)
            try:
                s.set('A',1); solver = lib.chem_solver(3); stats = s.solve(5,solver=solver)
                equal(stats.final_time,5)
                for key,value_ in [('function_calls',18),('jacobian_updates',9),('internal_steps',9),('solves',27)]: equal(getattr(stats,key),value_)
                solver.parameters.h_min = 15; solver.parameters.max_steps = 6; s.set('A',1); stats = s.solve(5,solver=solver)
                equal(stats.final_time,5)
                for key,value_ in [('function_calls',2),('jacobian_updates',1),('internal_steps',1),('solves',3)]: equal(getattr(stats,key),value_)
                equal(solver.parameters.h_min,15); equal(solver.parameters.max_steps,6)
            finally: s.close()
    finally: m.close()


INIT_NAMES = ('ConsistentICsUnchanged','MildlyInconsistentICsCorrected','SeverelyInconsistentICsConverge',
    'ScaledConstraintRowsSolveIdenticallyThroughSolve','PureODESystemUnaffected','MultiCellSystems','SubsequentSolveCallsReinitialize',
    'SolverStateToStringNewValue','FailedInitializationRestoresCallerState','SuccessfulInitializationKeepsTheCorrection',
    'WeightedCorrectionIsInvariantToConstraintRowScaling','WeightedCorrectionUsesStateTolerances')


@case(*['constraint_initialization:ConstraintInitialization.'+name for name in INIT_NAMES])
def initialization(identifier):
    name = identifier.split('.')[-1]
    if name == 'SolverStateToStringNewValue': equal(lib.chem_status_string(9),b'Constraint Initialization Failed'); return
    if name == 'ScaledConstraintRowsSolveIdenticallyThroughSolve':
        results = []
        for scale in (1,1e-12,1e12):
            m = Model(['A','B','C']); m.reaction(['A'],[('B',1)],rate(A=.5)); m.constraint('C',[(n,scale) for n in ('A','B','C')],scale*1.3); m.finalize(); s = State(m)
            try:
                s.set('A',.7); s.set('B',.15); s.set('C',5); s.rates(); s.solve(.001,method=5)
                near(s.get('C'),1.3-.7-.15,1e-5*(1.3-.7-.15)); results.append(s.get('C'))
            finally: s.close(); m.close()
        for v in results: near(v,results[0],1e-5*(1.3-.7-.15))
        return
    if name.startswith('WeightedCorrection'):
        for scale in ((1e-12,1,1e12) if 'Invariant' in name else (1,)):
            for tight in ((True,) if 'Invariant' in name else (False,True)):
                m = Model(['X','Z']); m.constraint('Z',[('X',scale),('Z',-scale)],0); m.finalize(); s = State(m)
                try:
                    solver = lib.chem_solver(3); equal(solver.parameters.init_tolerance,.1)
                    s.tolerances(1e-8 if tight else 1e-6,1e-12 if tight else 1e-8)
                    original = 4 if 'Invariant' in name else 1+5e-8
                    s.set('X',1); s.set('Z',original); stats = Stats()
                    equal(lib.chem_initialize_constraints(C.byref(solver),s.p,C.byref(stats)),0)
                    if tight: near(s.get('Z'),1,1e-12)
                    else: equal(s.get('Z'),original)
                    if 'Uses' in name: equal(stats.constraint_init_iterations,2 if tight else 1)
                finally: s.close(); m.close()
        return
    if name in ('FailedInitializationRestoresCallerState','SuccessfulInitializationKeepsTheCorrection'):
        m = Model(['B','C']); m.constraint('C',[('B',1),('C',2)],2,equilibrium_reactants=1,name='squared'); m.finalize(); s = State(m)
        try:
            solver = lib.chem_solver(3); stats = Stats(); s.set('B',.5)
            if name.startswith('Failed'):
                solver.parameters.init_iterations = 1; s.set('C',100); s.rates()
                equal(lib.chem_initialize_constraints(C.byref(solver),s.p,C.byref(stats)),9); equal(s.get('C'),100); equal(s.get('B'),.5)
            else:
                s.set('C',1e-6); s.rates(); equal(lib.chem_initialize_constraints(C.byref(solver),s.p,C.byref(stats)),0); near(s.get('C'),1,1e-3)
                for backtracks in (0,2):
                    solver.parameters.init_backtracks = backtracks; s.set('C',1e-6)
                    equal(lib.chem_initialize_constraints(C.byref(solver),s.p,C.byref(stats)),9); equal(s.get('B'),.5); equal(s.get('C'),1e-6)
        finally: s.close(); m.close()
        return
    pure = name == 'PureODESystemUnaffected'; m = simple_system(not pure)
    rows = [[1,.5,1],[2,.3,5],[.5,.1,100]] if name == 'MultiCellSystems' else [[1,0]] if pure else [[1,.5,5 if name.startswith('Mildly') else 1000 if name.startswith('Severely') else 1]]
    s = State(m,len(rows))
    try:
        for c,row in enumerate(rows):
            for i,v in enumerate(row): s.data.variables[i*s.data.stride+c] = v
            s.data.conditions[c].temperature = 298.15
        s.rates(); stats = s.solve(.01 if pure else .001)
        if pure: equal(stats.constraint_init_iterations,0)
        else:
            if name.startswith('Consistent'): assert stats.constraint_init_iterations <= 1
            if name.startswith('Mildly'): assert stats.constraint_init_iterations > 0
            if name == 'SubsequentSolveCallsReinitialize': s.set('C',999); s.rates(); s.solve(.001)
            for c in range(s.data.cells): near(2*s.get('B',c)-s.get('C',c),0,1e-5 if name == 'MultiCellSystems' else 1e-6)
    finally: s.close(); m.close()
