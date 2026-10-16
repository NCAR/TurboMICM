"""Terminator mechanism at twenty solar longitudes and one to four cells."""
from ffi import *
from harness import group
from cases_analytical import add_user


@group('terminator:RosenbrockSolver.','Terminator VectorTerminator')
def terminator(_):
    rng = np.random.default_rng(735)
    for nc in range(1,5):
        with Model() as m:
            m.species('Cl2',atol=1e-20); m.species('Cl',atol=1e-20)
            m.reaction(['Cl','Cl'],{'Cl2':1},rate_default(0)); add_user(m,['Cl2'],{'Cl':2},0,'toy_k1')
            m.finalize(reorder=True)
            with State(m,nc) as s:
                initial = rng.lognormal(-2,2,(nc,2)).T*np.array([[1e-6],[1e-10]])
                s.tolerance(1e-8,1e-20); s.conditions(298,101300,42)
                params = rosenbrock_default(3); params.max_steps = 100000
                for step in range(int(np.floor(2*np.pi/.3))):
                    lat = np.pi/180*np.arange(nc)*(90/nc)
                    k = np.maximum(0,np.sin(lat)*np.sin(np.pi*20/180)+np.cos(lat)*np.cos(np.pi*20/180)*np.cos(step*.3-np.pi*300/180))
                    s.set('Cl2',initial[0]); s.set('Cl',initial[1]); s.parameter('toy_k1',k); s.update(); s.solve(30,params=params)
                    r = k/4; ci,di = initial[1],initial[0]; total = ci+2*di
                    det = np.sqrt(r*r+2*r*total); e = np.exp(-4*det*30)
                    l = np.full(nc,4.); active = det*30 > 1e-16; l[active] = (1-e[active])/det[active]/30
                    cf = -l*(ci-det+r)*(ci+det+r)/(1+e+30*l*(ci+r))
                    target = np.stack((di-15*cf,ci+30*cf))
                    actual = s.y[[m.indices['Cl2'],m.indices['Cl']]]
                    near(actual,target,rtol=1e-6,atol=1e-14)
