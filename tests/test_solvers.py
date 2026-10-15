"""Direct ports of the solver unit-test quantities."""
# Copyright (C) 2023-2026 University Corporation for Atmospheric Research
# SPDX-License-Identifier: Apache-2.0
from chemffi import *
from test_constraints import near
import math

CASES={}

def case(suite,*names):
    def decorate(fn):
        for name in names: CASES[f'{suite}.{name}']=fn
        return fn
    return decorate


def five_species(cells):
    return Mechanism(['foo','bar','baz','quz','quuz'],reactions=[
        ([0,2],[(1,1),(4,2.4)],rate(0,[2e-11,0,110,300,0])),
        ([1],[(0,1),(3,1.4)],rate(0,[1e-6,0,0,300,0])),
        ([3],[],rate(0,[3.5e-6,0,0,300,0]))],cells=cells)


@case('RosenbrockSolver','StandardNormalizedError','VectorNormalizedError')
def normalized():
    for cells in [1,2,3,4,5]:
        with five_species(cells) as m:
            old=(D*(m.n*cells))();new=(D*(m.n*cells))();errors=(D*(m.n*cells))()
            expected=0
            for c in range(cells):
                for j in range(m.n):
                    index=j*cells+c
                    old[index]=7.7*c+j;new[index]=-13.9/(j+1)-c;errors[index]=81.57/(c+7)/(j+3)
                    scale=1e-3+1e-6*max(abs(old[index]),abs(new[index]))
                    expected+=errors[index]*errors[index]/(scale*scale)
            expected=max(math.sqrt(expected/(m.n*cells)),1e-10)
            computed=lib.chem_normalized_error(m.state,old,new,errors)
            assert abs(computed-expected)/max(abs(computed),abs(expected)) <= 1e-11


@case('RosenbrockSolver','StandardNormalizedErrorWithConstraints','VectorNormalizedErrorWithConstraints')
def normalized_constraints():
    for cells in [1,2,3,4,5]:
        with Mechanism(['A','B','C'],reactions=[([0],[(1,1)],rate(0,[.5,0,0,300,0]))],
                constraints=[equilibrium(2,[(1,1)],[(2,1)],10,-2400)],cells=cells) as m:
            assert lib.chem_mechanism_data(m.handle).contents.constraint_count==1
            old=(D*(3*cells))();new=(D*(3*cells))();errors=(D*(3*cells))()
            expected=0
            for c in range(cells):
                for j in range(3):
                    k=j*cells+c
                    old[k]=1+c+.1*j;new[k]=.8+.5*c+.2*j;errors[k]=.01*(1+c+j)
                    scale=1e-3+1e-6*max(abs(old[k]),abs(new[k]))
                    expected+=errors[k]*errors[k]/(scale*scale)
            expected=math.sqrt(expected/(3*cells))
            near(lib.chem_normalized_error(m.state,old,new,errors),expected,1e-12)


@case('BackwardEuler','CanCallSolve')
def backward_solve():
    rx=[([0],[(1,1)],rate(0,[2.15e-11,0,110,300,0])),([1],[(2,1)],rate(0,[3.3e-11,0,55,300,0]))]
    with Mechanism(['A','B','C'],reactions=rx) as m:
        m.set_y([[1,0,0]])
        lib.chem_temperatures(m.state)[0]=272.5;lib.chem_pressures(m.state)[0]=101253.3
        for j in range(3): lib.chem_absolute_tolerances(m.state)[j]=1e-6
        assert lib.chem_update_rates(m.state)==0
        p=lib.chem_backward_euler_parameters();stats=Stats()
        assert lib.chem_solve_backward_euler(m.state,C.byref(p),1,C.byref(stats))==0


@case('BackwardEuler','IsConverged')
def backward_converged():
    with Mechanism(['A','B','C'],cells=4) as m:
        p=lib.chem_backward_euler_parameters();p.small=1e-6
        lib.chem_relative_tolerance(m.state,1e-3)
        for j in range(3):lib.chem_absolute_tolerances(m.state)[j]=1e-6
        residual=(D*12)();y=(D*12)()
        def converged():return bool(lib.chem_backward_euler_converged(m.state,C.byref(p),residual,y))
        assert converged()
        residual[1*4+0]=1e-5;assert not converged()
        p.small=1e-4;assert converged()
        residual[2*4+3]=1e-3;assert not converged()
        y[2*4+3]=10;assert converged()
        residual[2*4+3]=1e-1;assert not converged()
        lib.chem_absolute_tolerances(m.state)[2]=1;assert converged()

@case('RosenbrockSolver','CanSetTolerances')
def initial_tolerances():
    for cells in range(1,11):
        with Mechanism(['foo','bar'],reactions=[([0],[(1,1)],rate(0,[2e-11,0,110,300,0]))],
                atols=[1e-7,1e-8],cells=cells) as m:
            assert list(lib.chem_absolute_tolerances(m.state)[:2])==[1e-7,1e-8]


@case('RosenbrockSolver','StandardAlphaMinusJacobian','VectorAlphaMinusJacobian')
def alpha_minus_jacobian():
    for cells in range(1,5):
        with five_species(cells) as m:
            count=lib.chem_jacobian_size(m.handle)
            assert count>=13
            values=(D*(count*cells))()
            for cell in range(cells):
                entries={(0,0):-12.2,(0,1):-24.3*(cell+2),(0,2):-42.3,(1,0):-.43,(1,1):-23.4,
                    (1,2):-83.4/(cell+3),(2,0):-4.74,(2,2):-6.91,(3,1):-59.1,(3,3):-83.4,
                    (4,0):-78.5,(4,2):-53.6,(4,4):-1}
                for (row,col),v in entries.items():
                    k=lib.chem_jacobian_index(m.handle,row,col)
                    assert k<Z(-1).value
                    values[k*cells+cell]=v
            lib.chem_alpha_minus_jacobian(m.handle,cells,values,42.042)
            for cell in range(cells):
                entries={(0,0):-12.2,(0,1):-24.3*(cell+2),(0,2):-42.3,(1,0):-.43,(1,1):-23.4,
                    (1,2):-83.4/(cell+3),(2,0):-4.74,(2,2):-6.91,(3,1):-59.1,(3,3):-83.4,
                    (4,0):-78.5,(4,2):-53.6,(4,4):-1}
                for (row,col),v in entries.items():
                    expected=42.042+v if row==col else v
                    assert values[lib.chem_jacobian_index(m.handle,row,col)*cells+cell]==expected


@case('RosenbrockSolver','RejectedStepAlphaMatchesInPlaceSolver')
def rejected_alpha():
    results=[]
    for linear_solver in [lib.chem_linear_mozart(),lib.chem_linear_doolittle()]:
        with five_species(1) as m:
            m.set_y([[1e12]*5]);assert lib.chem_update_rates(m.state)==0
            assert lib.chem_set_linear_solver(m.state,linear_solver)==0
            p=lib.chem_rosenbrock_parameters(3);p.h_start=1e3
            status,stats=m.solve(1e3,parameters=p)
            assert status==0,status
            assert stats.internal_steps-stats.accepted>=2
            results.append((stats.internal_steps,stats.accepted,stats.rejected))
    assert results[0]==results[1],results


def simple_constrained(cells=1):
    return Mechanism(['A','B','C'],reactions=[([0],[(1,1)],rate(0,[.5,0,0,300,0]))],
        constraints=[equilibrium(2,[(1,1)],[(2,1)],2)],cells=cells)


def initialize(m,p=None):
    p=p or lib.chem_rosenbrock_parameters(3)
    stats=Stats()
    status=lib.chem_initialize_constraints(m.state,C.byref(p),C.byref(stats))
    m.rows()
    return status,stats


@case('ConstraintInitialization','ConsistentICsUnchanged')
def init_consistent():
    with simple_constrained() as m:
        m.set_y([[1,.5,1]]);assert lib.chem_update_rates(m.state)==0
        status,stats=m.solve(.001)
        assert status==0 and stats.constraint_init_iterations<=1


@case('ConstraintInitialization','MildlyInconsistentICsCorrected')
def init_mild():
    with simple_constrained() as m:
        m.set_y([[1,.5,5]]);assert lib.chem_update_rates(m.state)==0
        status,stats=m.solve(.001)
        assert status==0 and stats.constraint_init_iterations>0
        row=m.rows()[0];near(2*row[1]-row[2],0,1e-6)


@case('ConstraintInitialization','SeverelyInconsistentICsConverge')
def init_severe():
    with simple_constrained() as m:
        m.set_y([[1,.5,1000]]);assert lib.chem_update_rates(m.state)==0
        status,stats=m.solve(.001)
        assert status==0
        row=m.rows()[0];near(2*row[1]-row[2],0,1e-6)


@case('ConstraintInitialization','ScaledConstraintRowsSolveIdenticallyThroughSolve')
def init_row_scaled():
    values=[]
    for scale in [1,1e-12,1e12]:
        with Mechanism(['A','B','C'],reactions=[([0],[(1,1)],rate(0,[.5,0,0,300,0]))],
                constraints=[linear(2,[(0,scale),(1,scale),(2,scale)],1.3*scale)]) as m:
            m.set_y([[.7,.15,5]]);assert lib.chem_update_rates(m.state)==0
            status,stats=m.solve(.001,method=-4)
            assert status==0,status
            values.append(m.rows()[0][2])
    consistent=1.3-.7-.15
    near(values[0],consistent,1e-5*consistent)
    for value in values[1:]:near(value,values[0],1e-5*consistent)


@case('ConstraintInitialization','PureODESystemUnaffected')
def init_ode():
    with Mechanism(['A','B'],reactions=[([0],[(1,1)],rate(0,[.5,0,0,300,0]))]) as m:
        m.set_y([[1,0]]);assert lib.chem_update_rates(m.state)==0
        status,stats=m.solve(.01)
        assert status==0 and stats.constraint_init_iterations==0


@case('ConstraintInitialization','MultiCellSystems')
def init_multicell():
    with simple_constrained(3) as m:
        m.set_y([[1,.5,1],[2,.3,5],[.5,.1,100]])
        assert lib.chem_update_rates(m.state)==0
        status,stats=m.solve(.001);assert status==0
        for row in m.rows():near(2*row[1]-row[2],0,1e-5)


@case('ConstraintInitialization','SubsequentSolveCallsReinitialize')
def init_repeat():
    with simple_constrained() as m:
        m.set_y([[1,.5,1]]);assert lib.chem_update_rates(m.state)==0
        assert m.solve(.001)[0]==0
        m.y[2]=999
        assert lib.chem_update_rates(m.state)==0
        assert m.solve(.001)[0]==0
        row=m.rows()[0];near(2*row[1]-row[2],0,1e-6)


@case('ConstraintInitialization','SolverStateToStringNewValue')
def init_string():
    assert lib.chem_status_string(8)==b'Constraint Initialization Failed'


@case('ConstraintInitialization','FailedInitializationRestoresCallerState')
def init_rollback():
    with Mechanism(['B','C'],constraints=[equilibrium(1,[(0,1)],[(1,2)],2)]) as m:
        m.set_y([[.5,100]])
        p=lib.chem_rosenbrock_parameters(3);p.constraint_init_max_iterations=1
        assert initialize(m,p)[0]==8
        assert m.rows()==[[.5,100]]


@case('ConstraintInitialization','SuccessfulInitializationKeepsTheCorrection')
def init_cold_start():
    with Mechanism(['B','C'],constraints=[equilibrium(1,[(0,1)],[(1,2)],2)]) as m:
        m.set_y([[.5,1e-6]])
        p=lib.chem_rosenbrock_parameters(3)
        assert initialize(m,p)[0]==0
        near(m.rows()[0][1],1,1e-3)
        for backtracks in [0,2]:
            p.constraint_init_max_backtracks=backtracks
            m.set_y([[.5,1e-6]])
            assert initialize(m,p)[0]==8
            assert m.rows()==[[.5,1e-6]]


@case('ConstraintInitialization','WeightedCorrectionIsInvariantToConstraintRowScaling')
def init_weighted_scale():
    for scale in [1e-12,1,1e12]:
        with Mechanism(['X','Z'],constraints=[linear(1,[(0,scale),(1,-scale)],0)]) as m:
            lib.chem_relative_tolerance(m.state,1e-8)
            for i in range(2):lib.chem_absolute_tolerances(m.state)[i]=1e-12
            m.set_y([[1,4]])
            assert initialize(m)[0]==0
            near(m.rows()[0][1],1,1e-12)


@case('ConstraintInitialization','WeightedCorrectionUsesStateTolerances')
def init_tolerances():
    assert lib.chem_rosenbrock_parameters(3).constraint_init_tolerance==.1
    for tight in [False,True]:
        with Mechanism(['X','Z'],constraints=[linear(1,[(0,1),(1,-1)],0)]) as m:
            lib.chem_relative_tolerance(m.state,1e-8 if tight else 1e-6)
            for i in range(2):lib.chem_absolute_tolerances(m.state)[i]=1e-12 if tight else 1e-8
            m.set_y([[1,1+5e-8]])
            status,stats=initialize(m)
            assert status==0
            if tight:
                near(m.rows()[0][1],1,1e-12);assert stats.constraint_init_iterations==2
            else:
                assert m.rows()[0][1]==1+5e-8 and stats.constraint_init_iterations==1

@case('Solver','GetLambdaRateConstantByNameCanOverrideLambda')
def lambda_override():
    signature=C.CFUNCTYPE(D,D,D,D,DP,P)
    original=signature(lambda t,p,air,custom,context:1e-3*t)
    replacement=signature(lambda t,p,air,custom,context:2e-3*t)
    r=rate(11);r.name=b'lambda_rc';r.function=C.cast(original,P)
    with Mechanism(['A','B'],reactions=[([0],[(1,1)],r)],atols=[1e-6,1e-6]) as m:
        handle=lib.chem_rate_by_name(m.handle,b'lambda_rc')
        assert bool(handle)
        handle.contents.function=C.cast(replacement,P)
        m.set_y([[1,0]])
        lib.chem_temperatures(m.state)[0]=300
        lib.chem_pressures(m.state)[0]=101253.3
        lib.chem_air_densities(m.state)[0]=1e6
        assert lib.chem_update_rates(m.state)==0
        near(lib.chem_rate_constants(m.state)[0],.6,1e-12)
        p=lib.chem_backward_euler_parameters();stats=Stats()
        assert lib.chem_solve_backward_euler(m.state,C.byref(p),1,C.byref(stats))==0
        assert stats.final_time==1
