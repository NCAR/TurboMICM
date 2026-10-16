"""Conservation, integrated flux, and shared-model OpenMP integration."""
import subprocess
from ffi import *
from harness import case
from cases_analytical import add_user,robertson,advance


def chapman_model():
    m = Model(['O','O1D','O2','O3','M','Ar','N2','H2O','CO2'])
    for rs,ps,a,c in [(['O1D','N2'],{'O':1,'N2':1},2.15e-11,110),
                     (['O1D','O2'],{'O':1,'O2':1},3.3e-11,55),
                     (['O','O3'],{'O2':2},8e-12,-2060),
                     (['O','O2','M'],{'O3':1,'M':1},6e-34,2.4)]:
        r = rate_default(0); r.p[0],r.p[2] = a,c; m.reaction(rs,ps,r)
    add_user(m,['O2'],{'O':2},.1,'jO2')
    add_user(m,['O3'],{'O1D':1,'O2':1},.2,'jO3a')
    add_user(m,['O3'],{'O':1,'O2':1},.3,'jO3b')
    return m.finalize(reorder=True)


@case('chapman_integration:ChapmanIntegration.CanBuildChapmanSystem')
def chapman(_):
    with chapman_model() as m:
        for lu in (0,1):
            with State(m) as s:
                require(state_set_lu(s.handle,lu)); s.y[:,0] = [.1,.1,.1,.2,.2,.2,.3,.3,.3]
                s.conditions(2,3); s.update()
                initial = s.y.copy()
                oxygen = lambda: sum(w*s.y[m.indices[n],0] for n,w in [('O',1),('O1D',1),('O2',2),('O3',3)])
                total = oxygen()
                for _ in range(100):
                    s.update(); s.solve(30.)
                    assert np.all(np.isfinite(s.y)) and np.all(s.y >= 0)
                    near(oxygen(),total,rtol=1e-8)
                    for name in ('M','Ar','N2','H2O','CO2'):
                        near(s.y[m.indices[name]],initial[m.indices[name]],rtol=1e-8)
                assert s.y[m.indices['O2'],0] < .5*initial[m.indices['O2'],0]
                assert s.y[m.indices['O'],0] > max(initial[m.indices['O'],0],.5*total)


@case('integrated_reaction_rates:ChapmanIntegration.CanBuildChapmanSystem')
def integrated_rates(_):
    with Model(['a','b','c','irr_1','irr_2']) as m:
        r = rate_default(0); r.p[0],r.p[2] = 2.15e-11,110
        m.reaction(['a'],{'b':1,'irr_1':1},r)
        add_user(m,['b'],{'c':1,'irr_2':1},1,'r2'); m.finalize(reorder=True)
        with State(m) as s:
            s.set('a',1); s.conditions(273,1000)
            for t in range(100):
                if t > 50: s.parameter('r2',0)
                s.update(); old = s.y.copy(); s.solve(30)
                for n in ('irr_1','irr_2'): assert s.y[m.indices[n],0] >= old[m.indices[n],0]


@case('openmp:OpenMP.OneSolverManyStates')
def openmp(_):
    result = subprocess.run([str(ROOT/'bin/test_native_cpu'),'openmp'],capture_output=True,text=True)
    assert result.returncode == 0, result.stdout+result.stderr
