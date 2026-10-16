"""Chapman unit and regression cases through mechanism-generic C operations."""
import math
from functools import lru_cache
from ffi import *
from harness import case,group
from cases_analytical import add_user,advance
import chapman_oracle as oracle


def model(reorder=False):
    m = Model(oracle.NAMES)
    for rs,ps,a,b,c in [(['O1D','N2'],{'O':1,'N2':1},2.15e-11,0,110),
        (['O1D','O2'],{'O':1,'O2':1},3.3e-11,0,55),(['O','O3'],{'O2':2},8e-12,0,-2060),
        (['O','O2','M'],{'O3':1,'M':1},6e-34,2.4,0)]:
        r = rate_default(0); r.p[:3] = a,b,c; m.reaction(rs,ps,r)
    add_user(m,['O2'],{'O':2},1e-4,'jO2'); add_user(m,['O3'],{'O1D':1,'O2':1},1e-5,'jO3a')
    add_user(m,['O3'],{'O':1,'O2':1},1e-6,'jO3b')
    return m.finalize(reorder=reorder)


@group('chapman_ode_solver:ChapmanMechanismHardCodedAndGeneral.','''DefaultConstructor lin_solve ReactionNames
PhotolysisNames SpeciesNames simple_force smaller_force factored_alpha_minus_jac dforce_dy_times_vector
Solve solve_10_times_larger solve_with_random_number_densities''')
@case('chapman_ode_solver:RegressionChapmanODESolver.solve_10_times_smaller')
def unit(id_):
    name = id_.split('.')[-1]
    if name == 'DefaultConstructor': assert rosenbrock_default(3).stages == 3; return
    if name in ('lin_solve','factored_alpha_minus_jac','dforce_dy_times_vector'):
        pat = np.zeros((9,9),np.uint8)
        for i,j in oracle.POSITIONS: pat[i,j] = 1
        sp = sparse_create(9,pat.ctypes.data_as(BP)); assert sp
        try:
            a = np.ones(sparse_size(sp)); mapping = [sparse_index(sp,i,j) for i,j in oracle.POSITIONS]
            assert len(a) == 23
            if name == 'lin_solve':
                b = np.full(9,.5); lu_solve(sp,1,ptr(a),ptr(b)); near(b,[.5,.5,.5,.5,.5,-.5,-1,.5,0],rtol=0)
            elif name == 'factored_alpha_minus_jac':
                a[:] = -1
                for i in range(9): a[sparse_index(sp,i,i)] += 2
                require(lu_factor(sp,1,ptr(a),0))
                # The reference stores reciprocal U diagonals; the generic API stores U.
                packed = a[mapping].copy()
                for q,(i,j) in enumerate(oracle.POSITIONS):
                    if i == j: packed[q] = 1/packed[q]
                near(packed,[1,-1,-1,-1,1,1,1,1,-1,-1,1,-1,1,-1,-1,-1,-2,-1,3,-1,-2,-3,.125],rtol=0,atol=.01)
            else:
                for i in range(5): a[sparse_index(sp,i,i)] = 0
                b = np.full(23,.5); result = np.zeros(23)
                sparse_multiply(sp,1,ptr(a),ptr(b),ptr(result))
                near(result,[0,0,0,0,0,2,3,2,2]+[0]*14,rtol=0,atol=.01)
        finally: sparse_free(sp)
        return
    with model() as m:
        if name == 'ReactionNames': assert model_reaction_count(m.handle) == 7; return
        if name == 'PhotolysisNames': assert model_parameter_count(m.handle) == 3; return
        if name == 'SpeciesNames': assert m.names == oracle.NAMES; return
        with State(m) as s:
            if name.endswith('force'):
                small = name == 'smaller_force'; s.y[:] = 5e-6 if small else 1.; s.rates[:] = 3e-8 if small else 1.
                s.conditions(273.15,1e5,6e-14 if small else 0)
                target = [0]*5+([1.49e-13,4.55e-13,1.5e-13,-3e-13] if small else [-1,3,2,-2])
                near(s.forcing()[:,0],target,rtol=0,atol=.01 if small else 0); return
            scale = 10 if name.endswith('larger') else .1 if name.endswith('smaller') else 1
            s.y[:,0] = np.array([1,.392,.0169,0,32.9,0,0,8.84,0])*scale
            if 'random' in name: s.y[:,0] = np.random.default_rng(3).random(9).astype(np.float32)
            s.conditions(273.15,1e5,2.7e19); s.update(); s.tolerance(1e-4,.001)
            p = rosenbrock_default(3); p.h_start = 1e-6; p.max_steps = 100
            s.solve(1.,params=p)
            if 'random' not in name:
                target = np.array([1,.392,.0169,0,32.9,0,.00176789,8.83912,0])*scale
                target[5],target[8] = 1.8039e-41*scale**3,4.5031e-36*scale**3
                near(s.y[:,0],target,rtol=0,atol=1e-4)


@group('regression_test_rosenbrock_p_force:RegressionRosenbrock.','RateConstants VectorRateConstants Forcing VectorForcing')
@group('regression_test_rosenbrock_dforce_dy:RegressionRosenbrock.','Jacobian VectorJacobian')
def regression_functions(id_):
    name = id_.split('.')[-1]
    with model(reorder=True) as m, State(m,3) as s:
        if name.endswith('RateConstants'):
            t,p = [284.19,215.02,299.31],[101245.,100789.2,101398.]
            photo = np.array([[1e-4,1e-5,1e-6],[3.2e-4,7.3e-5,3.2e-6],[5.2e-4,8.2e-5,4.6e-6]])
            s.conditions(t,p); s.custom[:] = photo.T; s.update()
            near(s.rates,np.array([oracle.rates(v,r) for v,r in zip(t,photo)]).T,rtol=0)
            return
        rng = np.random.default_rng(739); sigma = 4 if name.endswith('Jacobian') else 2
        s.y[:] = rng.lognormal(-2,sigma,s.y.shape); s.rates[:] = rng.lognormal(-2,sigma,s.rates.shape)
        if name.endswith('Forcing'):
            actual = s.forcing()
            for c in range(3):
                expected = oracle.forcing(s.rates[:,c],s.y[[m.indices[n] for n in oracle.NAMES],c])
                v = actual[[m.indices[n] for n in oracle.NAMES],c]
                assert np.all(abs(v-expected) <= (abs(v)+abs(expected))*1e-8+1e-12)
        else:
            actual = s.jacobian()
            for c in range(3):
                expected = oracle.jacobian(s.rates[:,c],s.y[[m.indices[n] for n in oracle.NAMES],c])
                for i,j in oracle.POSITIONS:
                    a,b = actual[c,m.indices[oracle.NAMES[i]],m.indices[oracle.NAMES[j]]],expected[i,j]
                    assert abs(a-b) <= (abs(a)+abs(b))*1e-8+1e-12


@group('regression_test_solve:RegressionRosenbrock.','TwoStageSolve ThreeStageSolve FourStageSolve FourStageDASolve SixStageDASolve VectorSolve')
def regression_solve(id_):
    method = {'TwoStageSolve':2,'ThreeStageSolve':3,'FourStageSolve':4,'FourStageDASolve':43,'SixStageDASolve':6,'VectorSolve':3}[id_.split('.')[-1]]
    with model(reorder=True) as m:
        initial = np.random.default_rng(193).lognormal(-2,2,(3,9))
        photo = np.array([[1e-4,1e-5,1e-6],[3.2e-4,7.3e-5,3.2e-6],[5.2e-4,8.2e-5,4.6e-6]])
        temperature,pressure = [284.19,215.02,299.31],[101245.,100789.2,101398.]
        expected = [oracle.solve(initial[c,[m.indices[n] for n in oracle.NAMES]],oracle.rates(temperature[c],photo[c]),500.) for c in range(3)]
        for lu in (0,1):
            with State(m,3) as s:
                require(state_set_lu(s.handle,lu)); s.y[:] = initial.T; s.custom[:] = photo.T
                s.conditions(temperature,pressure); s.update(); s.solve(500.,method)
                for c in range(3):
                    actual = s.y[[m.indices[n] for n in oracle.NAMES],c]; target = expected[c]
                    tol = .01 if method == 2 else 1e-4
                    assert np.all(abs(actual-target) <= (abs(actual)+abs(target))*tol+initial.max()*1e-12), (method,actual,target)
