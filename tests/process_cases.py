"""Reaction assembly, rate storage, parameterizations, forcing, and Jacobian ports."""
import math
import random
from api import *
from harness import case, equal, near, real_equal
from constraint_cases import verify_jacobian, jac_value

RateFn = C.CFUNCTYPE(D,C.POINTER(Conditions),DP,P)


@case(*['reaction_rate_store:ReactionRateConstantStore.'+name for name in (
    'OffsetsAreContiguousCumulativeSizes','ArrheniusParametersPreserved','BranchedDerivedFieldsComputed',
    'UserDefinedCustomParamIndex','SurfaceDataFieldsAndCustomParamIndex','SurfaceCustomParamIndexAfterUserDefined',
    'LambdaEntriesRcIndex','ParameterizedMultipliers','TotalRateConstantsMatchesProcessCount')])
def rate_store(identifier):
    name = identifier.split('.')[-1]
    callback = RateFn(lambda c,u,ctx:c.contents.temperature*.001)
    multiplier = RateFn(lambda c,u,ctx:2*c.contents.air_density)
    m = Model(['a','b','c'])
    if name == 'ParameterizedMultipliers':
        m.close(); m = Model(['a']); index = m.species('b',parameterized=2); m.data.species[index].parameterization = C.cast(multiplier,P)
        m.reaction(['a'],[('b',1)]); m.reaction(['b'],[('a',1)])
    elif name == 'ArrheniusParametersPreserved': m.reaction(['a'],[('b',1)],rate(A=2.15e-4,B=1.2,C=110,D=300))
    elif name == 'BranchedDerivedFieldsComputed': m.reaction(['a'],[('b',1)],rate(BRANCHED,X=1,Y=0,a0=.5,n=3))
    elif name == 'UserDefinedCustomParamIndex':
        m.reaction(['a'],[('b',1)],rate(USER,label='rate1',A=2)); m.reaction(['a'],[('b',1)],rate(USER,label='rate2',A=.5))
    elif name == 'SurfaceDataFieldsAndCustomParamIndex': m.reaction(['c'],[('b',1)],rate(SURFACE,label='surf',diffusion=230,molecular_weight=.025,probability=.74))
    elif name == 'SurfaceCustomParamIndexAfterUserDefined':
        m.reaction(['a'],[('b',1)],rate(USER,label='ud')); m.reaction(['c'],[('b',1)],rate(SURFACE,label='surf',diffusion=1e-5,molecular_weight=.025))
    elif name == 'LambdaEntriesRcIndex':
        m.reaction(['a'],[('b',1)]); m.reaction(['a'],[('b',1)],rate(LAMBDA,label='lam',function=C.cast(callback,P)))
    else:
        for r in ([rate(),rate(),rate(TROE),rate(TUNNELING),rate(USER,label='ud1'),rate(SURFACE,label='surf1',diffusion=1e-5,molecular_weight=.025)] if name.startswith('Offsets') else
                  [rate(),rate(TROE),rate(USER,label='ud'),rate(SURFACE,label='surf',diffusion=1e-5,molecular_weight=.025),rate(LAMBDA,label='lam',function=C.cast(callback,P))]):
            m.reaction(['a'],[('b',1)],r)
    m.finalize()
    try:
        if name.startswith('Offsets'):
            equal(m.data.nreactions,6)
            for kind,expected in [(TROE,2),(TERNARY,3),(JPL19,3),(BRANCHED,3),(TUNNELING,3),(TAYLOR,4),(REVERSIBLE,4),(USER,4),(SURFACE,5),(LAMBDA,6)]: equal(lib.chem_rate_offset(m.p,kind),expected)
            equal([m.data.reactions[i].rate.type for i in range(6)],[ARRHENIUS,ARRHENIUS,TROE,TUNNELING,USER,SURFACE])
        elif name == 'ArrheniusParametersPreserved':
            r = m.data.reactions[0].rate
            for key,v in [('A',2.15e-4),('B',1.2),('C',110),('D',300)]: real_equal(getattr(r,key),v)
        elif name == 'BranchedDerivedFieldsComputed':
            r = m.data.reactions[0].rate; k0 = 2e-22*6.02214076e23*1e-6*math.exp(3)
            near(r.k0,k0,1e-10*k0)
            a,b = k0*(2.45e19/6.02214076e23*1e6),.43*(293/298)**-8
            z = a/(1+a/b)*.41**(1/(1+math.log10(a/b)**2))
            near(r.z,z,1e-10*abs(z))
        elif name == 'UserDefinedCustomParamIndex':
            for i,scale in enumerate([2,.5]): equal(m.data.reactions[i].rate.custom_index,i); real_equal(m.data.reactions[i].rate.A,scale)
        elif name == 'SurfaceDataFieldsAndCustomParamIndex':
            r = m.data.reactions[0].rate; real_equal(r.diffusion,230); real_equal(r.probability,.74); equal(r.custom_index,0)
            near(r.mean_free_speed_factor,8*(1.380649e-23*6.02214076e23)/(math.pi*.025),1e-14)
        elif name == 'SurfaceCustomParamIndexAfterUserDefined': equal(m.data.reactions[0].rate.custom_index,0); equal(m.data.reactions[1].rate.custom_index,1)
        elif name == 'LambdaEntriesRcIndex':
            equal(lib.chem_rate_offset(m.p,LAMBDA),1); r = m.data.reactions[1].rate; assert r.function
            near(value(r,temperature=300),.3,1e-14)
        elif name == 'ParameterizedMultipliers':
            s = State(m)
            try:
                s.data.conditions[0] = Conditions(300,0,5); s.rates(); near(s.data.rates[1],10,1e-14)
            finally: s.close()
        else: equal(lib.chem_rate_offset(m.p,LAMBDA)+sum(m.data.reactions[i].rate.type == LAMBDA for i in range(m.data.nreactions)),m.data.nreactions)
    finally: m.close()


@case('process:Process.Matrix','process:Process.VectorMatrix')
def process_rates(_):
    callback = RateFn(lambda c,u,ctx:.82*c.contents.air_density)
    m = Model(['foo']); index = m.species('bar',parameterized=2); m.data.species[index].parameterization = C.cast(callback,P)
    m.reaction(['foo','bar'],[],rate(A=12.2,C=300))
    m.reaction(['foo'],[],rate(SURFACE,label='foo_surf',molecular_weight=.025,diffusion=230))
    m.reaction(['bar'],[],rate(USER,label='bar_user')); m.finalize(); s = State(m,5)
    try:
        rng = random.Random(1)
        for c in range(5):
            t,p,d = [rng.lognormvariate(0,.01)*v for v in (285,101100,10)]
            user,radius,number = [rng.lognormvariate(0,.01)*v for v in (1e-2,1e-8,1e5)]
            s.data.conditions[c] = Conditions(t,p,d)
            for label,val in [('bar_user',user),('foo_surf.effective radius [m]',radius),('foo_surf.particle number concentration [# m-3]',number)]:
                index = lib.chem_custom_id(m.p,label.encode()); s.data.custom[index*5+c] = val
        s.rates()
        for c in range(5):
            cond = s.data.conditions[c]; user = s.data.custom[c]; radius,number = s.data.custom[5+c],s.data.custom[10+c]
            expected = 12.2*math.exp(300/cond.temperature)*cond.air_density*.82
            near(s.data.rates[c],expected,1e-12*expected)
            expected = user*(cond.air_density*.82); near(s.data.rates[5+c],expected,1e-12*expected)
            speed = math.sqrt(8*(1.380649e-23*6.02214076e23)/(math.pi*.025)*cond.temperature)
            expected = 4*number*math.pi*radius*radius/(radius/230+4/speed)
            near(s.data.rates[10+c],expected,1e-10*expected)
    finally: s.close(); m.close()


@case('process:Process.BuildsChemicalReaction','process:Process.ChemicalReactionCopyAssignmentSucceeds')
def reaction_assembly(identifier):
    m = Model(['O3','NO','NO2','O2'])
    try:
        source = rate(A=12.2,C=300); m.reaction(['O3','NO'],[('NO2',1),('O2',1)],source)
        r = m.data.reactions[0]; equal(r.nreactants,2); equal(m.data.species[r.products[0]].name,b'NO2')
        m.reaction(list(r.reactants[:r.nreactants]),list(zip(r.products[:r.nproducts],r.yields[:r.nproducts])),r.rate)
        r = m.data.reactions[0]
        copy = m.data.reactions[1]; equal(copy.nreactants,r.nreactants); equal(copy.nproducts,r.nproducts); equal(copy.rate.A,12.2)
        source.A = 9; equal(copy.rate.A,12.2)
    finally: m.close()


def basic_process_model():
    m = Model(['foo','bar','baz','quz','quuz','corge']); index = m.species('qux',parameterized=2)
    callback = RateFn(lambda c,u,ctx:.72*c.contents.air_density); m.callbacks.append(callback); m.data.species[index].parameterization = C.cast(callback,P)
    for rs,ps in [(['foo','baz'],[('bar',1),('quuz',2.4)]),(['bar','qux'],[('foo',1),('quz',1.4)]),(['quz'],[]),(['baz','qux'],[('bar',1),('quz',2.5)])]:
        m.reaction(rs,ps,rate(A=12.2,C=300))
    return m.finalize()


@case('process_set:ProcessSet.Matrix','process_set:ProcessSet.VectorMatrix')
def basic_process(_):
    m = basic_process_model(); s = State(m,2)
    try:
        rows = [[.1,.2,.3,.4,.5,0],[1.1,1.2,1.3,1.4,1.5,0]]
        constants = [[10,20*70*.72,30,40*70*.72],[110,120*80*.72,130,140*80*.72]]
        for c,row in enumerate(rows):
            for i,v in enumerate(row): s.data.variables[i*2+c] = v
            for i,v in enumerate(constants[c]): s.data.rates[i*2+c] = v
        f = arr([1000]*12); lib.chem_generic_forcing(m.p,2,2,s.data.variables,s.data.rates,f)
        jac = arr([100]*(m.data.nnz*2)); lib.chem_generic_jacobian(m.p,2,2,s.data.variables,s.data.rates,jac)
        expected_pairs = {(0,0),(0,1),(0,2),(1,0),(1,1),(1,2),(2,0),(2,2),(3,1),(3,2),(3,3),(4,0),(4,2)}
        used = {m.data.species[x].name.decode() for rx in m.data.reactions[:m.data.nreactions] for x in list(rx.reactants[:rx.nreactants])+list(rx.products[:rx.nproducts])}
        equal(used,{'foo','bar','baz','quz','quuz','qux'})
        for c,(y,k) in enumerate(zip(rows,constants)):
            r1,r2,r3,r4 = k[0]*y[0]*y[2],k[1]*y[1],k[2]*y[3],k[3]*y[2]
            expected = [1000-r1+r2,1000+r1-r2+r4,1000-r1-r4,1000+1.4*r2-r3+2.5*r4,1000+2.4*r1,1000]
            for i,v in enumerate(expected): real_equal(f[i*2+c],v)
            derivatives = {(0,0):k[0]*y[2],(0,1):-k[1],(0,2):k[0]*y[0],(1,0):-k[0]*y[2],(1,1):k[1],(1,2):-k[0]*y[0]-k[3],
                           (2,0):k[0]*y[2],(2,2):k[0]*y[0]+k[3],(3,1):-1.4*k[1],(3,2):-2.5*k[3],(3,3):k[2],(4,0):-2.4*k[0]*y[2],(4,2):-2.4*k[0]*y[0]}
            equal(set(derivatives),expected_pairs)
            for (i,j),v in derivatives.items(): real_equal(jac_value(m,s,jac,i,j,c),100+v)
    finally: s.close(); m.close()


@case('process_set:RandomProcessSet.Matrix')
def random_process(_):
    for cells,reactions,species in [(200,50,40),(300,30,20),(400,100,80)]:
        rng = random.Random(1); m = Model([str(i) for i in range(species)])
        for _ in range(reactions):
            rs = [rng.randrange(species) for _ in range(rng.randrange(4))]
            ps = [(rng.randrange(species),1.2) for _ in range(rng.randrange(11))]; m.reaction(rs,ps,rate(A=12.2,C=300))
        m.finalize(); s = State(m,cells)
        try:
            for i in range(species*cells): s.data.variables[i] = rng.lognormvariate(-2,4)
            for i in range(reactions*cells): s.data.rates[i] = rng.lognormvariate(-2,4)
            f = arr([1000]*(species*cells)); lib.chem_generic_forcing(m.p,cells,cells,s.data.variables,s.data.rates,f)
            assert all(math.isfinite(v) for v in f)
        finally: s.close(); m.close()


@case('process_set:ProcessSetAlgebraicVariables.CudaMatrix')
def masking(_):
    m = Model(['A','B','C','D']); m.reaction(['A','B'],[('C',1),('D',1)]); m.finalize(); s = State(m,2)
    try:
        m.data.algebraic[1] = m.data.algebraic[3] = 1
        for c,values in enumerate([[1,2,3,4],[5,6,7,8]]):
            for i,v in enumerate(values): s.data.variables[i*2+c] = v
            s.data.rates[c] = 10 if c == 0 else 20
        f = arr([1000]*8); jac = arr([500]*(m.data.nnz*2))
        lib.chem_generic_forcing(m.p,2,2,s.data.variables,s.data.rates,f)
        lib.chem_generic_jacobian(m.p,2,2,s.data.variables,s.data.rates,jac)
        for c in range(2):
            equal(f[c],1000-[20,600][c]); equal(f[4+c],1000+[20,600][c])
            equal(f[2+c],1000); equal(f[6+c],1000)
            for row in (1,3):
                for col in (0,1): equal(jac_value(m,s,jac,row,col,c),500)
            for row,sign in [(0,1),(2,-1)]:
                for col,derivatives in [(0,[20,120]),(1,[10,100])]: equal(jac_value(m,s,jac,row,col,c),500+sign*derivatives[c])
    finally: s.close(); m.close()


@case('process_set:ProcessSetFiniteDifferenceJacobian.Matrix')
def process_fd(_):
    m = Model(['A','B','C'])
    m.reaction(['A','B'],[('C',2)],rate(A=1)); m.reaction(['C'],[('A',1)],rate(A=1)); m.finalize()
    s = State(m,2)
    try:
        for name,values in [('A',[1,.5]),('B',[2,4]),('C',[3,1.5])]: s.set(name,values)
        for i,values in enumerate([[5,8],[10,12]]):
            for c,v in enumerate(values): s.data.rates[i*2+c] = v
        verify_jacobian(s)
    finally: s.close(); m.close()
