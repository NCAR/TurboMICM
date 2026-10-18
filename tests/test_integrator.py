"""C/prototype equivalence, generic kinetics, analytic and failure tests."""
import csv
from dataclasses import asdict
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'prototypes'))
from adaptive import Integrator, Options, rodas31, COUNTERS
from ffi import *


class OptionsC(C.Structure):
    _fields_ = [(name,I if field.type is int else D)
                for name,field in Options.__dataclass_fields__.items()]


class Extra(C.Structure):
    _fields_ = [(name,Z) for name in COUNTERS[7:]]


integrator_default = bind('integrator_default',OptionsC,I)
integrator_create = bind('integrator_create',P,P,C.POINTER(OptionsC))
integrator_free = bind('integrator_free',None,P)
integrator_reset = bind('integrator_reset',None,P)
integrator_step = bind('integrator_step',I,P,D,C.POINTER(Stats),C.POINTER(Extra))


def retained_table(method):
    if method == 431:
        return rodas31()
    if method == 102:
        return None
    p = rosenbrock_default(method)
    return dict(stages=p.stages,gamma=p.gamma[0],order=p.estimator_order,
                new=list(p.new_function),a=list(p.a),c=list(p.c),m=list(p.m),e=list(p.e))


def model():
    m = Model(['A','B','C'])
    j = m.parameter('light',.4)
    photo = rate_default(7); photo.custom_index = j
    m.reaction(['A'],{'B':1},photo)
    rate = rate_default(0); rate.p[0] = 3
    m.reaction(['B','B'],{'B':1,'C':1},rate)
    rate.p[0] = 4
    m.reaction(['B','C'],{'A':1,'C':1},rate)
    return m.finalize()


def equivalence():
    rows = []
    for method in (2,3,4,43,6,102,431):
        for variant in ('exact','pi','reuse','filter','maxnorm'):
            o = Options(method=method)
            if variant == 'pi': o.controller = 1
            if variant == 'reuse': o.jacobian_steps,o.lu_band = 3,2.
            if variant == 'filter': o.controller,o.filter_b,o.filter_k = 2,2.,2.
            if variant == 'maxnorm': o.norm = 2
            oc = OptionsC(**asdict(o))
            with model() as m, State(m,2) as s:
                s.y[:] = [[1,.7],[.02,.03],[.01,.1]]
                initial_mass = s.y.sum(axis=0)
                s.conditions(); s.tolerance(rtol=1e-4,atol=1e-10)
                light = .4

                def f(y):
                    a,b,c = y
                    r1,r2,r3 = light*a,3*b*b,4*b*c
                    return np.array([-r1+r3,r1-r2-r3,r2])

                def jac(y):
                    a,b,c = y
                    return np.array([[[-light,4*c[k],4*b[k]],
                                      [light,-6*b[k]-4*c[k],-4*b[k]],
                                      [0,6*b[k],0]] for k in range(y.shape[1])])

                oracle = Integrator(f,jac,s.y,o,retained_table(method),rtol=1e-4,atol=1e-10)
                work = integrator_create(s.handle,C.byref(oc)); assert work
                try:
                    for light,dt in ((.4,.5),(.2,.4),(.7,.3)):
                        s.parameter('light',light); s.update()
                        st,ex = Stats(),Extra()
                        require(integrator_step(work,dt,C.byref(st),C.byref(ex)))
                        expected,counts = oracle.step(dt)
                        # Accumulated double precision roundoff; not a loose ODE error tolerance.
                        near(s.y,expected,rtol=8e-13,atol=3e-15)
                        near(s.y.sum(axis=0),initial_mass,rtol=2e-14)
                        assert np.all(s.y >= 0)
                        actual = {k:getattr(st,k) if hasattr(st,k) else getattr(ex,k) for k in COUNTERS}
                        assert actual == counts, (method,variant,actual,counts)
                        rows.append(dict(method=method,variant=variant,dt=dt,
                                         max_difference=float(np.max(abs(s.y-expected))),**actual))
                    before = s.y.copy()
                    require(integrator_step(work,0,C.byref(st),C.byref(ex)))
                    near(s.y,before,rtol=0)
                    assert integrator_step(work,-1,C.byref(st),C.byref(ex)) == -1
                    assert integrator_step(work,float('nan'),C.byref(st),C.byref(ex)) == -1
                finally:
                    integrator_free(work)
    return rows


def analytical():
    # Different mechanism and extreme stiffness: irreversible linear chain.
    for method in (2,3,4,43,6,102,431):
        with Model(['X','Y','Z']) as m:
            rate = rate_default(0); rate.p[0] = 1e5
            m.reaction(['X'],{'Y':1},rate)
            rate.p[0] = .1; m.reaction(['Y'],{'Z':1},rate); m.finalize()
            with State(m,3) as s:
                s.conditions(); s.y[:] = [[1,.2,.8],[0,0,0],[0,0,0]]; s.update()
                s.tolerance(rtol=1e-5,atol=1e-12)
                o = integrator_default(method)
                work = integrator_create(s.handle,C.byref(o)); assert work
                st,ex = Stats(),Extra()
                try:
                    require(integrator_step(work,2,C.byref(st),C.byref(ex)))
                    y = 1e5/(1e5-.1)*np.exp(-.2)
                    expected = np.array([0,y,1-y])[:,None]*np.array([1,.2,.8])
                    near(s.y,expected,rtol=2e-4,atol=2e-11)
                    near(s.y.sum(axis=0),[1,.2,.8],rtol=2e-13)
                    assert st.final_time == 2 and np.all(s.y >= 0)
                    # A changed state is legal after resetting the predictor.
                    integrator_reset(work)
                finally:
                    integrator_free(work)
    # Invalid settings fail before they can change caller-owned concentrations.
    with model() as m, State(m) as s:
        for key,value in [('method',999),('jacobian_steps',0),('factor_min',0),('h_start',float('nan')),
                          ('norm',3),('max_steps',0),('lu_band',.5),('newton_tolerance',1)]:
            o = integrator_default(43); setattr(o,key,value)
            assert not integrator_create(s.handle,C.byref(o))


if __name__ == '__main__':
    rows = equivalence()
    analytical()
    if len(sys.argv) == 2:
        with open(sys.argv[1],'w',newline='') as file:
            writer = csv.DictWriter(file,fieldnames=list(rows[0]))
            writer.writeheader(); writer.writerows(rows)
    print(f'PASS: {len(rows)} C/prototype comparisons with identical counters; seven stiff analytical cases; invalid options')
