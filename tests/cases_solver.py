"""Solver unit cases and initialization transaction tests."""
from ffi import *
from harness import case,group


def five_species():
    m = Model(['foo','bar','baz','quz','quuz'])
    for r,p,a,c in [(['foo','baz'],{'bar':1,'quuz':2.4},2e-11,110),
                    (['bar'],{'foo':1,'quz':1.4},1e-6,0),(['quz'],{},3.5e-6,0)]:
        rate = rate_default(0); rate.p[0],rate.p[2] = a,c
        m.reaction(r,p,rate)
    return m.finalize()


@group('rosenbrock:RosenbrockSolver.', '''StandardAlphaMinusJacobian VectorAlphaMinusJacobian CanSetTolerances
StandardNormalizedError VectorNormalizedError StandardNormalizedErrorWithConstraints VectorNormalizedErrorWithConstraints
RejectedStepAlphaMatchesInPlaceSolver''')
def rosenbrock(id_):
    name = id_.split('.')[-1]
    if name == 'RejectedStepAlphaMatchesInPlaceSolver':
        results = []
        with five_species() as m:
            for algo in (0,1):
                with State(m) as s:
                    require(state_set_lu(s.handle,algo))
                    s.y[:] = 1e12; s.conditions(298.15,0,0); s.update()
                    p = rosenbrock_default(3); p.h_start = 1e3
                    _,stats = s.solve(1e3,params=p)
                    assert stats.internal_steps-stats.accepted >= 2
                    results.append((stats.internal_steps,stats.accepted,stats.rejected))
        assert results[0] == results[1]
        return
    if name == 'CanSetTolerances':
        with Model() as m:
            m.species('foo',atol=1e-7); m.species('bar',atol=1e-8)
            r = rate_default(0); r.p[0],r.p[2] = 2e-11,110.
            m.reaction(['foo'],{'bar':1},r); m.finalize()
            for nc in range(1,11):
                with State(m,nc) as s:
                    zero = np.zeros((2,nc)); error = np.ones((2,nc))
                    actual = normalized_error(s.handle,ptr(zero),ptr(zero),ptr(error))
                    near(actual,np.sqrt((1e14+1e16)/2),rtol=1e-15)
        return
    constrained = name.endswith('WithConstraints')
    m = Model(['A','B','C']) if constrained else five_species()
    with m:
        if constrained:
            r = rate_default(0); r.p[0] = .5
            m.reaction(['A'],{'B':1},r); m.equilibrium('C',{'B':1},{'C':1},10.,-2400.)
            m.finalize()
        for nc in range(1,6):
            with State(m,nc) as s:
                if 'Alpha' in name:
                    values = np.zeros((model_jacobian_size(m.handle),nc))
                    for c in range(nc):
                        entries = {(0,0):-12.2,(0,1):-24.3*(c+2),(0,2):-42.3,(1,0):-.43,(1,1):-23.4,
                            (1,2):-83.4/(c+3),(2,0):-4.74,(2,2):-6.91,(3,1):-59.1,(3,3):-83.4,
                            (4,0):-78.5,(4,2):-53.6,(4,4):-1.}
                        for (i,j),v in entries.items():
                            idx = model_jacobian_index(m.handle,i,j); assert idx >= 0
                            values[idx,c] = v
                    expected = values.copy()
                    for i in range(5): expected[model_jacobian_index(m.handle,i,i)] += 42.042
                    alpha_minus_jacobian(s.handle,42.042,ptr(values))
                    near(values,expected,rtol=0)
                else:
                    i,j = np.meshgrid(np.arange(nc),np.arange(len(m.names)))
                    if constrained:
                        old,new,errors = 1+i+.1*j,.8+.5*i+.2*j,.01*(1+i+j)
                    else:
                        old,new,errors = 7.7*i+j,-13.9/(j+1)-i,81.57/(i+7)/(j+3)
                    old,new,errors = [np.ascontiguousarray(a,dtype=float) for a in (old,new,errors)]
                    expected = max(np.sqrt(np.mean(errors**2/(.001+1e-6*np.maximum(abs(old),abs(new)))**2)),1e-10)
                    actual = normalized_error(s.handle,ptr(old),ptr(new),ptr(errors))
                    near(actual,expected,rtol=0 if constrained else 1e-11,atol=1e-12 if constrained else 0)


INIT_NAMES = '''ConsistentICsUnchanged MildlyInconsistentICsCorrected SeverelyInconsistentICsConverge
ScaledConstraintRowsSolveIdenticallyThroughSolve PureODESystemUnaffected MultiCellSystems
SubsequentSolveCallsReinitialize SolverStateToStringNewValue FailedInitializationRestoresCallerState
SuccessfulInitializationKeepsTheCorrection WeightedCorrectionIsInvariantToConstraintRowScaling
WeightedCorrectionUsesStateTolerances'''


@group('constraint_initialization:ConstraintInitialization.',INIT_NAMES)
def initialize(id_):
    name = id_.split('.')[-1]
    p,stats = rosenbrock_default(3),Stats()
    if name == 'SolverStateToStringNewValue':
        assert status_string(-11) == b'Constraint Initialization Failed'
        return
    if name in ('FailedInitializationRestoresCallerState','SuccessfulInitializationKeepsTheCorrection'):
        with Model(['B','C']) as m:
            m.equilibrium('C',{'B':1},{'C':2},2.,name='squared'); m.finalize()
            with State(m) as s:
                s.conditions(298.15); s.update(); s.set('B',.5)
                if name.startswith('Failed'):
                    p.init_max_iterations = 1; s.set('C',100.)
                    assert initialize_constraints(s.handle,C.byref(p),C.byref(stats)) == -11
                    near(s.y,[[.5],[100.]],rtol=0)
                else:
                    s.set('C',1e-6)
                    require(initialize_constraints(s.handle,C.byref(p),C.byref(stats)))
                    near(s.y[1,0],1.,rtol=0,atol=1e-3)
                    for backtracks in (0,2):
                        p.init_max_backtracks = backtracks; s.set('C',1e-6)
                        assert initialize_constraints(s.handle,C.byref(p),C.byref(stats)) == -11
                        near(s.y,[[.5],[1e-6]],rtol=0)
        return
    if name.startswith('WeightedCorrection'):
        for scale in (1e-12,1.,1e12) if 'Invariant' in name else (1.,):
            with Model(['X','Z']) as m:
                m.linear('Z',{'X':scale,'Z':-scale},0); m.finalize()
                with State(m) as s:
                    if 'Invariant' in name:
                        s.y[:] = [[1],[4]]; s.tolerance(1e-8,1e-12)
                        require(initialize_constraints(s.handle,C.byref(p),C.byref(stats)))
                        near(s.y[1,0],1.,rtol=0,atol=1e-12)
                    else:
                        assert p.init_tolerance == .1
                        for rtol,atol,iters in [(1e-6,1e-8,1),(1e-8,1e-12,2)]:
                            s.y[:] = [[1],[1+5e-8]]; s.tolerance(rtol,atol); stats = Stats()
                            require(initialize_constraints(s.handle,C.byref(p),C.byref(stats)))
                            assert stats.constraint_init_iterations == iters
                            near(s.y[1,0],1+5e-8 if iters == 1 else 1.,rtol=0,atol=0 if iters == 1 else atol)
        return
    if name == 'ScaledConstraintRowsSolveIdenticallyThroughSolve':
        values = []
        for scale in (1.,1e-12,1e12):
            with Model(['A','B','C']) as m:
                r = rate_default(0); r.p[0] = .5
                m.reaction(['A'],{'B':1},r)
                m.linear('C',{'A':scale,'B':scale,'C':scale},1.3*scale); m.finalize()
                with State(m) as s:
                    s.y[:] = [[.7],[.15],[5]]; s.conditions(298.15); s.update(); s.solve(.001,43)
                    values.append(s.y[2,0]); near(values[-1],.45,rtol=0,atol=1e-5*.45)
        near(np.array(values),np.full(3,values[0]),rtol=0,atol=1e-5*.45)
        return
    with Model(['A','B','C']) as m:
        r = rate_default(0); r.p[0] = .5
        m.reaction(['A'],{'B':1},r)
        if name != 'PureODESystemUnaffected': m.equilibrium('C',{'B':1},{'C':1},2.)
        m.finalize()
        with State(m,3 if name == 'MultiCellSystems' else 1) as s:
            s.conditions(298.15); s.update()
            s.y[:] = np.array([[1,.5,1],[2,.3,5],[.5,.1,100]]).T if s.cells == 3 else [[1],[.5],[1]]
            if name.startswith('Mildly'): s.set('C',5.)
            if name.startswith('Severely'): s.set('C',1000.)
            _,stats = s.solve(.001)
            if name == 'ConsistentICsUnchanged': assert stats.constraint_init_iterations <= 1
            if name == 'PureODESystemUnaffected': assert stats.constraint_init_iterations == 0
            else:
                assert stats.constraint_init_iterations > 0
                near(s.y[2],2*s.y[1],rtol=0,atol=1e-6)
            if name == 'SubsequentSolveCallsReinitialize':
                s.set('C',999.); s.update(); s.solve(.001)
                near(s.y[2],2*s.y[1],rtol=0,atol=1e-6)


@group('backward_euler:BackwardEuler.','CanCallSolve IsConverged')
def backward_euler(id_):
    with Model(['A','B','C']) as m:
        for r,p,a,c in [(['A'],{'B':1},2.15e-11,110),(['B'],{'C':1},3.3e-11,55)]:
            rate = rate_default(0); rate.p[0],rate.p[2] = a,c
            m.reaction(r,p,rate)
        m.finalize()
        with State(m,4) as s:
            s.conditions(272.5,101253.3,1e6); s.set('A',1.); s.tolerance(atol=1e-6); s.update()
            if id_.endswith('CanCallSolve'):
                s.solve(1,'be'); return
            p = backward_euler_default(); p.small = 1e-6; s.tolerance(.001,1e-6)
            residual,y = np.zeros_like(s.y),np.zeros_like(s.y)
            check = lambda: backward_euler_converged(s.handle,ptr(residual),ptr(y),C.byref(p))
            assert check()
            residual[1,0] = 1e-5; assert not check()
            p.small = 1e-4; assert check()
            residual[2,3] = 1e-3; assert not check()
            y[2,3] = 10.; assert check()
            residual[2,3] = .1; assert not check()
            s.tolerance(.001,[1e-6,1e-6,1.]); assert check()


@case('solver_lambda_rate_constant:Solver.GetLambdaRateConstantByNameCanOverrideLambda')
def override_lambda(_):
    fun_type = C.CFUNCTYPE(D,C.POINTER(Conditions),DP,P)
    callback = fun_type(lambda c,p,x: 1e-3*c.contents.temperature)
    replacement = fun_type(lambda c,p,x: 2e-3*c.contents.temperature)
    with Model(['A','B']) as m:
        r = rate_default(10); r.function = C.cast(callback,P)
        index = m.reaction(['A'],{'B':1},r); m.finalize()
        model_rate(m.handle,index).contents.function = C.cast(replacement,P)
        with State(m) as s:
            s.tolerance(atol=1e-6); s.set('A',1.); s.conditions(300.,101253.3,1e6); s.update()
            near(s.rates[0,0],.6,rtol=0,atol=1e-12)
            _,stats = s.solve(1.,'be'); assert stats.final_time == 1.
