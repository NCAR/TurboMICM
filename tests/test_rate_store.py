"""Reaction-rate storage metadata and condition-dependent multipliers."""
from functools import partial
import math
from systemffi import *
from reference_random import lognormals


class Derived(C.Structure):
    _fields_ = [('branch_k0',D),('branch_z',D),('surface_speed_factor',D)]


bind('chem_rate_derived',Derived,C.POINTER(Rate))
for suffix in ['count','offset']:bind('chem_rate_family_'+suffix,Z,P,I)
bind('chem_rate_family_index',Z,P,I,Z)
bind('chem_rate_specification',C.POINTER(Rate),P,Z)


def multiplier(r,p,name='parameterized'):
    r.multiplier_count = 1
    r.multipliers = (Parameterization*1)(p)
    r.multiplier_names = (C.c_char_p*1)(name.encode())
    return r


def surface(diff=1e-5,mw=.025,probability=1.,index=0,status=0):
    with Species('c',[] if mw is None else [('molecular weight [kg mol-1]',0,mw)]) as s:
        member = PhaseSpecies(s.handle,diff or 0,diff is not None)
        r = Rate()
        assert lib.chem_surface_rate(C.byref(member),probability,index,C.byref(r)) == status
        return r


def mechanism(rates,parameters=(),cells=1):
    return Mechanism(['a','b','c'],[([0],[(1,1)],r) for r in rates],parameters=parameters,cells=cells)


def counts():
    ud = rate(9);ud.parameter = 0
    with mechanism([rate(0),rate(0),rate(1),rate(4),ud,surface(index=1)],['ud','radius','number']) as m:
        expected = [2,1,0,0,0,1,0,0,1,1,0]
        assert [lib.chem_rate_family_count(m.handle,i) for i in range(11)] == expected
        assert [lib.chem_rate_family_offset(m.handle,i) for i in range(1,11)] == [2,3,3,3,3,4,4,4,5,6]
        assert lib.chem_rate_family_offset(m.handle,10)+lib.chem_rate_family_count(m.handle,10) == 6


def arrhenius():
    r = rate(0,[2.15e-4,1.2,110.,300.,0.])
    with mechanism([r]) as m:
        assert lib.chem_rate_family_count(m.handle,0) == 1
        stored = lib.chem_rate_specification(m.handle,0).contents
        assert list(stored.p[:4]) == [2.15e-4,1.2,110.,300.]


def branch():
    r = rate(5,[1.,0.,.5,3.])
    with mechanism([r]) as m:
        assert lib.chem_rate_family_count(m.handle,4) == 1
        derived = lib.chem_rate_derived(lib.chem_rate_specification(m.handle,0))
        k0 = 2e-22 * 6.02214076e23 * 1e-6 * math.exp(3.)
        a = k0 * (2.45e19 / 6.02214076e23 * 1e6)
        b = .43 * (293./298.)**-8
        z = a/(1+a/b) * .41**(1/(1+math.log10(a/b)**2))
        assert abs(derived.branch_k0-k0) <= 1e-10*k0
        assert abs(derived.branch_z-z) <= 1e-10*abs(z)


def user_indices():
    r1,r2 = rate(9,[2.]),rate(9,[.5]);r1.parameter,r2.parameter = 0,1
    with mechanism([r1,r2],['rate1','rate2']) as m:
        assert lib.chem_rate_family_count(m.handle,8) == 2
        for i,scale in [(0,2.),(1,.5)]:
            r = lib.chem_rate_specification(m.handle,i).contents
            assert r.parameter == i and r.p[0] == scale


def surface_data(after_user=False):
    r = surface(diff=1e-5 if after_user else 2.3e2,probability=1. if after_user else .74,index=int(after_user))
    rates = [rate(9),r] if after_user else [r]
    with mechanism(rates,['ud','radius','number'] if after_user else ['radius','number']) as m:
        assert lib.chem_rate_family_count(m.handle,9) == 1
        if after_user:
            assert lib.chem_rate_family_count(m.handle,8) == 1
            assert lib.chem_rate_specification(m.handle,0).contents.parameter == 0
        stored = lib.chem_rate_specification(m.handle,int(after_user))
        assert stored.contents.parameter == int(after_user)
        if not after_user:
            assert stored.contents.p[0] == 2.3e2 and stored.contents.p[2] == .74
            speed = lib.chem_rate_derived(stored).surface_speed_factor
            assert abs(speed - 8*(1.380649e-23*6.02214076e23)/(math.pi*.025)) <= 1e-14


def callback():
    fn = C.CFUNCTYPE(D,D,D,D,DP,P)(lambda t,p,air,custom,ctx:t*1e-3)
    r = rate(11);r.function = C.cast(fn,P);r.name = b'lam'
    with mechanism([rate(0),r]) as m:
        assert lib.chem_rate_family_count(m.handle,0) == 1 and lib.chem_rate_family_count(m.handle,10) == 1
        index = lib.chem_rate_family_index(m.handle,10,0)
        assert index == 1
        stored = lib.chem_rate_specification(m.handle,index)
        assert stored and stored.contents.function
        assert abs(lib.chem_rate_evaluate(stored,300,0,0,None)-.3) <= 1e-14


def multipliers():
    r = multiplier(rate(0),Parameterization(air_density=2.,has_value=1),'b')
    with mechanism([rate(0),r]) as m:
        assert [lib.chem_rate_specification(m.handle,i).contents.multiplier_count for i in range(2)] == [0,1]
        stored = lib.chem_rate_specification(m.handle,1)
        assert abs(lib.chem_rate_evaluate(stored,300,0,5,None)-10.) <= 1e-14
        # Verify the full state-update path too.
        lib.chem_air_densities(m.state)[0] = 5.
        assert lib.chem_update_rates(m.state) == 0
        assert abs(lib.chem_rate_constants(m.state)[1]-10.) <= 1e-14


def total():
    fn = C.CFUNCTYPE(D,D,D,D,DP,P)(lambda t,p,air,custom,ctx:1.)
    cb = rate(11);cb.function = C.cast(fn,P)
    with mechanism([rate(0),rate(1),rate(9),surface(index=1),cb],['ud','radius','number']) as m:
        assert lib.chem_rate_family_offset(m.handle,10)+lib.chem_rate_family_count(m.handle,10) == 5


def process(vector=False):
    # Source's lognormal(0,.01) sequence; changing the mean of the existing
    # generator retains its polar-normal sequence.
    numbers = lognormals(.01,mean=0.)
    conditions = [[next(numbers)*scale for scale in [285.,101100.,10.,1e-2,1e-8,1e5]] for _ in range(5)]
    for group in ([1,2,3,4] if vector else [5]):
        for start in range(0,5,group):
            rows = conditions[start:start+group]
            arr = multiplier(rate(0,[12.2,0,300,300,0]),Parameterization(air_density=.82,has_value=1),'bar')
            ud = multiplier(rate(9),Parameterization(air_density=.82,has_value=1),'bar')
            with mechanism([arr,ud,surface(diff=2.3e2,index=1)],['bar_user','radius','number'],len(rows)) as m:
                nc = len(rows)
                for c,(t,p,air,user,radius,number) in enumerate(rows):
                    lib.chem_temperatures(m.state)[c] = t
                    lib.chem_pressures(m.state)[c] = p
                    lib.chem_air_densities(m.state)[c] = air
                    m.custom[c],m.custom[nc+c],m.custom[2*nc+c] = user,radius,number
                assert lib.chem_update_rates(m.state) == 0
                for c,(t,p,air,user,radius,number) in enumerate(rows):
                    expected = [12.2*math.exp(300/t)*(air*.82),user*(air*.82),
                        4*number*math.pi*radius*radius/(radius/2.3e2+4/math.sqrt(8*(1.380649e-23*6.02214076e23)/(math.pi*.025)*t))]
                    for j,reference in enumerate(expected):
                        tol = (1e-10 if j==2 else 1e-12)*reference
                        assert abs(lib.chem_rate_constants(m.state)[j*nc+c]-reference) <= tol


def reaction(copy=False):
    names = ['foo','bar'] if copy else ['O3','NO','NO2','O2']
    r = rate(0,[12.2,0,300,300,0]) if copy else rate(0)
    with Mechanism(names,[([0,1],[] if copy else [(2,1.),(3,1.)],r)]) as m:
        stored = m.data.reactions[0]
        assert stored.reactant_count == 2
        if copy:
            spec = lib.chem_rate_specification(m.handle,0).contents
            assert spec.type == r.type and spec.p[0] == r.p[0]
            assert stored.product_count == m.data.reactions[0].product_count == 0
            assert m.data.species_names[stored.reactants[0]].decode() == 'foo'
        else:assert m.data.species_names[stored.products[0].species].decode() == 'NO2'
        with Species(names[0]) as s, Phase('gas',[PhaseSpecies(s.handle,0,0)]) as phase:
            assert phase.name == 'gas'


CASES = {
    'ReactionRateConstantStore.OffsetsAreContiguousCumulativeSizes':counts,
    'ReactionRateConstantStore.ArrheniusParametersPreserved':arrhenius,
    'ReactionRateConstantStore.BranchedDerivedFieldsComputed':branch,
    'ReactionRateConstantStore.UserDefinedCustomParamIndex':user_indices,
    'ReactionRateConstantStore.SurfaceDataFieldsAndCustomParamIndex':surface_data,
    'ReactionRateConstantStore.SurfaceCustomParamIndexAfterUserDefined':partial(surface_data,True),
    'ReactionRateConstantStore.LambdaEntriesRcIndex':callback,
    'ReactionRateConstantStore.ParameterizedMultipliers':multipliers,
    'ReactionRateConstantStore.TotalRateConstantsMatchesProcessCount':total,
    'ReactionRateConstantStore.SurfaceMissingDiffusionCoefficientThrows':partial(surface,diff=None,status=1),
    'ReactionRateConstantStore.SurfaceMissingMolecularWeightThrows':partial(surface,mw=None,status=9),
    'Process.Matrix':process,
    'Process.VectorMatrix':partial(process,True),
    'Process.BuildsChemicalReaction':reaction,
    'Process.ChemicalReactionCopyAssignmentSucceeds':partial(reaction,True),
}
