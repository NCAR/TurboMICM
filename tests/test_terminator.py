"""Latitude/longitude-dependent chlorine mechanism closed-form assertions."""
from functools import partial
import math
from chemffi import *
from reference_random import lognormals


def terminator(vector=False):
    for repeat in range(1 if vector else 2):
        for nc in [1,2,3,4]:
            photo=rate(9);photo.parameter=0
            with Mechanism(['Cl2','Cl'],[([0],[(1,2.)],photo),([1,1],[(0,1.)],rate(0))],parameters=['toy_k1'],cells=nc,atols=[1e-20]*2) as m:
                lib.chem_relative_tolerance(m.state,1e-8)
                numbers=lognormals()
                initial=[[next(numbers)*1e-6,next(numbers)*1e-10] for _ in range(nc)]
                parameters=lib.chem_rosenbrock_parameters(3);parameters.max_steps=100000
                for step in range(math.floor(2*math.pi/.3)):
                    m.set_y(initial)
                    for c in range(nc):
                        lat=math.pi/180*(c*(90./nc))
                        k1=max(0.,math.sin(lat)*math.sin(math.pi*20/180)+math.cos(lat)*math.cos(math.pi*20/180)*math.cos(step*.3-math.pi*300/180))
                        m.custom[c]=k1
                        lib.chem_temperatures(m.state)[c]=298.
                        lib.chem_pressures(m.state)[c]=101300.
                        lib.chem_air_densities(m.state)[c]=42.
                    assert lib.chem_update_rates(m.state)==0
                    status,stats=m.solve(30.,parameters=parameters)
                    assert status==0,(vector,nc,step,status,stats.final_time)
                    for c,(cl2_i,cl_i) in enumerate(initial):
                        r=m.custom[c]/4.
                        cly=cl_i+2*cl2_i
                        det=math.sqrt(r*r+2*r*cly)
                        e=math.exp(-4*det*30.)
                        l=(1-e)/det/30. if det*30.>1e-16 else 4.
                        cl_f=-l*(cl_i-det+r)*(cl_i+det+r)/(1+e+30*l*(cl_i+r))
                        reference=[cl2_i-30*cl_f/2,cl_i+30*cl_f]
                        for j,(a,b) in enumerate(zip(m.rows()[c],reference)):
                            assert abs(a-b)<=abs(b)*1e-6+1e-14,(vector,nc,step,c,j,a,b)


CASES={'RosenbrockSolver.Terminator':terminator,'RosenbrockSolver.VectorTerminator':partial(terminator,True)}
