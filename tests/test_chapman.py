"""Chapman case ports through generic reaction, matrix, and solver APIs."""
from functools import partial
import math
import struct
from chemffi import *
from test_linear import Pattern,FACTOR,SOLVE
from reference_random import Minstd

NAMES=['M','Ar','CO2','H2O','N2','O1D','O','O2','O3']
POSITIONS=[(0,0),(6,0),(7,0),(8,0),(1,1),(2,2),(3,3),(4,4),(5,4),(6,4),(5,5),(6,5),(6,6),(7,6),(8,6),(5,7),(6,7),(7,7),(8,7),(5,8),(6,8),(7,8),(8,8)]
PHOTOS=[[1e-4,1e-5,1e-6],[3.2e-4,7.3e-5,3.2e-6],[5.2e-4,8.2e-5,4.6e-6]]
TEMPERATURES=[284.19,215.02,299.31]
PRESSURES=[101245.,100789.2,101398.]
bind('chem_matrix_vector_product',None,P,Z,DP,DP,DP)


def reactions(integration=False):
    result=[]
    for rr,pp,a,b,c in [([5,4],[(6,1),(4,1)],2.15e-11,0,110),
                        ([5,7],[(6,1),(7,1)],3.3e-11,0,55),
                        ([6,8],[(7,2)],8e-12,0,-2060),
                        ([6,7,0],[(8,1),(0,1)],6e-34,0 if integration else 2.4,2.4 if integration else 0)]:
        r=rate(0);r.p[0],r.p[1],r.p[2]=a,b,c;result.append((rr,pp,r))
    for parameter,(rr,pp) in enumerate([([7],[(6,2)]),([8],[(5,1),(7,1)]),([8],[(6,1),(7,1)])]):
        r=rate(9);r.parameter=parameter;result.append((rr,pp,r))
    return result


def mechanism(cells=1,integration=False):
    return Mechanism(NAMES,reactions(integration),parameters=['jO2','jO3a','jO3b'],cells=cells)


def condition(m,temperatures=TEMPERATURES,pressures=PRESSURES,photos=PHOTOS):
    for c in range(m.nc):
        lib.chem_temperatures(m.state)[c]=temperatures[c]
        lib.chem_pressures(m.state)[c]=pressures[c]
        for j in range(3):m.custom[j*m.nc+c]=photos[c][j]
    assert lib.chem_update_rates(m.state)==0


def simple(kind):
    with mechanism() as m:
        if kind=='DefaultConstructor':assert lib.chem_rosenbrock_parameters(3).stages==3
        elif kind=='ReactionNames':assert m.data.reaction_count==7
        elif kind=='PhotolysisNames':assert m.data.parameter_count==3
        elif kind=='SpeciesNames':assert m.data.species_count==9
        else:
            small=kind=='smaller_force';m.set_y([[5e-6 if small else 1.]*9])
            for i in range(7):lib.chem_rate_constants(m.state)[i]=3e-8 if small else 1.
            lib.chem_air_densities(m.state)[0]=6e-14 if small else 0
            expected=[0]*5+([1.49e-13,4.55e-13,1.5e-13,-3e-13] if small else [-1,3,2,-2])
            for i,(got,want) in enumerate(zip(m.force()[0],expected)):
                assert abs(got-want)<=(.01 if small and i>=5 else 0),(i,got,want)


def matrix(kind):
    with Pattern(9,POSITIONS) as m:
        dense=[[0.]*9 for _ in range(9)]
        for i,j in POSITIONS:dense[i][j]=1.
        packed=m.pack([dense],0);linear=lib.chem_linear_mozart()
        if kind=='lin_solve':
            x=(D*9)(*([.5]*9));SOLVE(linear.solve)(m.handle,1,packed,x)
            assert list(x)==[.5,.5,.5,.5,.5,-.5,-1,.5,0]
        elif kind=='factored_alpha_minus_jac':
            for k in range(m.nnz):packed[k]=-packed[k]
            lib.chem_alpha_minus_jacobian(m.handle,1,packed,2.)
            assert FACTOR(linear.factor)(m.handle,1,packed)==0
            expected=[1,-1,-1,-1,1,1,1,1,-1,-1,1,-1,1,-1,-1,-1,-2,-1,3,-1,-2,-3,.125]
            for (i,j),want in zip(POSITIONS,expected):
                got=packed[m.indices[i][j]]
                if i==j:got=1/got # legacy storage uses reciprocal pivots
                assert abs(got-want)<=.01,(i,j,got,want)
        else:
            for i in range(5):packed[m.indices[i][i]]=0
            x=(D*23)(*([.5]*23));out=(D*23)()
            lib.chem_matrix_vector_product(m.handle,1,packed,x,out)
            for got,want in zip(out,[0]*5+[2,3,2,2]+[0]*14):assert abs(got-want)<=.01,(got,want)


def solve(scale=1,random=False):
    with mechanism() as m:
        values=[1,.392,.0169,0,32.9,0,0,8.84,0]
        if random:
            engine=Minstd()
            values=[struct.unpack('f',struct.pack('f',(engine()-engine.minimum)/engine.span))[0] for _ in range(9)]
        else:values=[v*scale for v in values]
        m.set_y([values]);condition(m,[273.15],[100000],[[1e-4,1e-5,1e-6]])
        lib.chem_air_densities(m.state)[0]=2.7e19
        lib.chem_relative_tolerance(m.state,1e-4)
        p=lib.chem_rosenbrock_parameters(3);p.max_steps=100;p.h_start=1e-6
        status,_=m.solve(1.,parameters=p)
        assert status==0,status
        if not random:
            expected=[1,.392,.0169,0,32.9,1.8039e-41,.00176789,8.83912,4.5031e-36]
            for i,(got,want) in enumerate(zip(m.rows()[0],expected)):
                want*=scale**3 if i in (5,8) else scale
                assert abs(got-want)<=1e-4,(got,want)


def integration():
    # The source assigns a raw row after Markowitz reordering. Use that same
    # permutation so every named initial concentration matches the source.
    from test_regression import source_permutation
    permutation=source_permutation(['O','O1D','O2','O3','M','Ar','N2','H2O','CO2'],True)
    initial=[0.]*9
    for old,value in zip(permutation,[.1,.1,.1,.2,.2,.2,.3,.3,.3]):initial[old]=value
    with mechanism(integration=True) as m:
        m.set_y([initial]);condition(m,[2.],[3.],[[.1,.2,.3]])
        def oxygen(row):return row[6]+row[5]+2*row[7]+3*row[8]
        total=oxygen(initial)
        for t in range(100):
            assert lib.chem_update_rates(m.state)==0
            status,_=m.solve(30.);assert status==0,(t,status)
            row=m.rows()[0]
            assert all(math.isfinite(v) and v>=0 for v in row)
            assert abs(oxygen(row)-total)<=1e-8*abs(total)
            for i in range(5):assert abs(row[i]-initial[i])<=1e-8*abs(initial[i])
        assert row[7]<.5*initial[7] and row[6]>initial[6] and row[6]>.5*total


def integrated_rates():
    r=rate(0);r.p[0],r.p[2]=2.15e-11,110
    user=rate(9);user.parameter=0
    with Mechanism(['a','b','c','irr_1','irr_2'],[([0],[(1,1),(3,1)],r),([1],[(2,1),(4,1)],user)],parameters=['r2']) as m:
        m.set_y([[1,0,0,0,0]]);lib.chem_temperatures(m.state)[0]=273;lib.chem_pressures(m.state)[0]=1000
        for t in range(100):
            m.custom[0]=0. if t>50 else 1.
            assert lib.chem_update_rates(m.state)==0
            previous=m.rows()[0];m.solve(30.);current=m.rows()[0]
            assert current[3]>=previous[3] and current[4]>=previous[4]


CASES={}
for name in ['DefaultConstructor','ReactionNames','PhotolysisNames','SpeciesNames','simple_force','smaller_force']:
    CASES['chapman_ode_solver:ChapmanMechanismHardCodedAndGeneral.'+name]=partial(simple,name)
for name in ['lin_solve','factored_alpha_minus_jac','dforce_dy_times_vector']:
    CASES['chapman_ode_solver:ChapmanMechanismHardCodedAndGeneral.'+name]=partial(matrix,name)
for name,scale,random in [('Solve',1,False),('solve_10_times_larger',10,False),('solve_with_random_number_densities',1,True)]:
    CASES['chapman_ode_solver:ChapmanMechanismHardCodedAndGeneral.'+name]=partial(solve,scale,random)
CASES['chapman_ode_solver:RegressionChapmanODESolver.solve_10_times_smaller']=partial(solve,.1)
CASES['chapman_integration:ChapmanIntegration.CanBuildChapmanSystem']=integration
CASES['integrated_reaction_rates:ChapmanIntegration.CanBuildChapmanSystem']=integrated_rates
