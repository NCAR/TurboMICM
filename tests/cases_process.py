"""Reaction construction, rate storage, forcing, and structural Jacobian checks."""
import math
from ffi import *
from harness import case,group


def parameterized(m,name,coefficients):
    s=species_create(name.encode());assert s
    try:
        require(species_parameterize(s,1,*coefficients))
        m.ids[name]=require(model_add_species_spec(m.handle,s))
    finally:species_free(s)


@group('process:Process.','Matrix VectorMatrix')
def process_rates(_):
    with Model() as m:
        m.species('foo',**{'molecular weight [kg mol-1]':.025,'diffusion coefficient [m2 s-1]':230})
        parameterized(m,'bar',[0,0,0,.82])
        user=m.parameter('bar_user');radius=m.parameter('foo_surf.effective radius [m]')
        m.parameter('foo_surf.particle number concentration [# m-3]')
        a=rate_default(0);a.p[0],a.p[2]=12.2,300
        u=rate_default(7);u.custom_index=user
        surf=rate_default(6);surf.custom_index=radius
        m.reaction(['foo','bar'],{},a);m.reaction(['bar'],{},u);m.reaction(['foo'],{},surf);m.finalize()
        with State(m,5) as s:
            rng=np.random.default_rng(0)
            t,p,air=rng.lognormal(0,.01,(3,5))*np.array([285,101100,10])[:,None]
            uv,rv,nv=rng.lognormal(0,.01,(3,5))*np.array([1e-2,1e-8,1e5])[:,None]
            s.conditions(t,p,air);s.parameter('bar_user',uv)
            s.parameter('foo_surf.effective radius [m]',rv);s.parameter('foo_surf.particle number concentration [# m-3]',nv)
            s.update()
            near(s.rates[0],12.2*np.exp(300/t)*.82*air,rtol=1e-12)
            near(s.rates[1],uv*.82*air,rtol=1e-12)
            expected=4*nv*math.pi*rv*rv/(rv/230+4/np.sqrt(8*(1.380649e-23*6.02214076e23)*t/(math.pi*.025)))
            near(s.rates[2],expected,rtol=1e-10)


@group('process:Process.','BuildsChemicalReaction ChemicalReactionCopyAssignmentSucceeds')
def process_build(_):
    with Model(['O3','NO','NO2','O2']) as m:
        rate=rate_default(0);rate.p[0],rate.p[2]=12.2,300
        reactants=ints([0,1]);products=ints([2,3]);yields=doubles([1,1])
        require(model_add_reaction(m.handle,2,reactants,2,products,yields,C.byref(rate)))
        reactants[0]=3;products[0]=0;rate.p[0]=999
        copy=model_clone_builder(m.handle);assert copy
        try:
            for handle in (m.handle,copy):
                info=ReactionInfo();require(model_reaction_info(handle,0,C.byref(info)))
                assert info.reactants==2 and info.products==2
                assert list(info.reactant_ids[:2])==[0,1] and list(info.product_ids[:2])==[2,3]
                assert species_name(model_species_info(handle,info.product_ids[0]))==b'NO2'
                assert info.rate.contents.p[0]==12.2
        finally:model_free(copy)


@group('reaction_rate_store:ReactionRateConstantStore.',
       'OffsetsAreContiguousCumulativeSizes TotalRateConstantsMatchesProcessCount')
def rate_groups(id_):
    with Model(['a','b']) as m:
        m.species('c',**{'molecular weight [kg mol-1]':.025,'diffusion coefficient [m2 s-1]':1e-5})
        u=m.parameter('ud1');radius=m.parameter('surf1.effective radius [m]');m.parameter('surf1.particle number concentration [# m-3]')
        types=[0,0,1,4,7,6] if 'Offsets' in id_ else [0,1,7,6,10]
        @C.CFUNCTYPE(D,C.POINTER(Conditions),DP,P)
        def lam(c,p,ctx):return 1.
        for kind in types:
            r=rate_default(kind)
            if kind==1:r.p[0]=r.p[3]=1
            if kind in (6,7):r.custom_index=radius if kind==6 else u
            if kind==10:r.function=C.cast(lam,P)
            m.reaction(['c' if kind==6 else 'a'],{'b':1},r)
        offsets=[];counts=[]
        for kind in (0,1,2,3,5,4,8,9,7,6,10):
            offset,count=Z(),Z();require(model_rate_group(m.handle,kind,C.byref(offset),C.byref(count)))
            offsets.append(offset.value);counts.append(count.value)
            for i in range(count.value):assert model_rate_in_group(m.handle,kind,i).contents.type==kind
            assert not model_rate_in_group(m.handle,kind,count.value)
        if 'Offsets' in id_:
            assert counts==[2,1,0,0,0,1,0,0,1,1,0]
            assert offsets==[0,2,3,3,3,3,4,4,4,5,6]
        assert offsets[-1]+counts[-1]==model_reaction_count(m.handle)==len(types)


@case('reaction_rate_store:ReactionRateConstantStore.ArrheniusParametersPreserved')
def arr_parameters(_):
    with Model(['a','b']) as m:
        rate=rate_default(0);rate.p[:5]=[2.15e-4,1.2,110,300,0]
        m.reaction(['a'],{'b':1},rate)
        near(list(model_rate(m.handle,0).contents.p[:5]),[2.15e-4,1.2,110,300,0],rtol=0)


@case('reaction_rate_store:ReactionRateConstantStore.BranchedDerivedFieldsComputed')
def branched_parameters(_):
    with Model(['a','b']) as m:
        rate=rate_default(5);rate.p[:5]=[1,0,.5,3,0]
        m.reaction(['a'],{'b':1},rate);stored=model_rate(m.handle,0).contents
        k0=2e-22*6.02214076e23*1e-6*math.exp(3)
        a=k0*2.45e19/6.02214076e23*1e6;b=.43*(293/298)**-8
        z=a/(1+a/b)*.41**(1/(1+math.log10(a/b)**2))
        near(list(stored.p[5:7]),[k0,z],rtol=1e-10)


@case('reaction_rate_store:ReactionRateConstantStore.UserDefinedCustomParamIndex')
def user_parameters(_):
    with Model(['a','b']) as m:
        for name,scale in (('rate1',2),('rate2',.5)):
            r=rate_default(7);r.custom_index=m.parameter(name);r.p[0]=scale
            m.reaction(['a'],{'b':1},r)
        for i,scale in enumerate((2,.5)):
            r=model_rate(m.handle,i).contents
            assert r.custom_index==i and r.p[0]==scale


@group('reaction_rate_store:ReactionRateConstantStore.','SurfaceDataFieldsAndCustomParamIndex SurfaceCustomParamIndexAfterUserDefined')
def surface_parameters(_):
    for has_user in (False,True):
        with Model(['a','b']) as m:
            diffusion=1e-5 if has_user else 230
            m.species('c',**{'molecular weight [kg mol-1]':.025,'diffusion coefficient [m2 s-1]':diffusion})
            if has_user:
                r=rate_default(7);r.custom_index=m.parameter('ud');m.reaction(['a'],{'b':1},r)
            r=rate_default(6);r.custom_index=m.parameter('surf.effective radius [m]');r.p[2]=.74
            m.parameter('surf.particle number concentration [# m-3]');index=m.reaction(['c'],{'b':1},r)
            stored=model_rate(m.handle,index).contents
            assert stored.p[0]==diffusion and stored.p[2]==.74 and stored.custom_index==int(has_user)
            near(np.array(stored.p[1]),np.array(8*(1.380649e-23*6.02214076e23)/(math.pi*.025)),rtol=0,atol=1e-14)


@case('reaction_rate_store:ReactionRateConstantStore.LambdaEntriesRcIndex')
def lambda_index(_):
    @C.CFUNCTYPE(D,C.POINTER(Conditions),DP,P)
    def callback(c,p,ctx):return c.contents.temperature*1e-3
    with Model(['a','b']) as m:
        m.reaction(['a'],{'b':1},rate_default(0))
        r=rate_default(10);r.function=C.cast(callback,P)
        assert m.reaction(['a'],{'b':1},r)==1
        r=model_rate(m.handle,1).contents
        assert r.function
        near(np.array(rate_calculate(C.byref(r),C.byref(Conditions(300,0,0)),None)),np.array(.3),rtol=0,atol=1e-14)


@case('reaction_rate_store:ReactionRateConstantStore.ParameterizedMultipliers')
def parameter_multipliers(_):
    with Model(['a']) as m:
        parameterized(m,'b',[0,0,0,2])
        r=rate_default(0);m.reaction(['a'],{'b':1},r);m.reaction(['b'],{'a':1},r);m.finalize()
        for i in (0,1):
            info=ReactionInfo();require(model_reaction_info(m.handle,i,C.byref(info)))
            assert info.parameterized_reactants==i
        with State(m) as s:
            s.conditions(300,0,5);s.update();near(s.rates[:,0],[1,10],rtol=0,atol=1e-14)


@group('reaction_rate_store:ReactionRateConstantStore.','SurfaceMissingDiffusionCoefficientThrows SurfaceMissingMolecularWeightThrows')
def missing_surface(_):
    for property_ in ('molecular weight [kg mol-1]','diffusion coefficient [m2 s-1]'):
        with Model(['b']) as m:
            m.species('c',**{property_:.025})
            r=rate_default(6);r.custom_index=m.parameter('surf.effective radius [m]')
            m.parameter('surf.particle number concentration [# m-3]')
            assert model_add_reaction(m.handle,1,ints([m.ids['c']]),1,ints([m.ids['b']]),doubles([1]),C.byref(r))==-3


@group('process_set:ProcessSet.','Matrix VectorMatrix')
def process_set(_):
    names=['foo','bar','baz','quz','quuz','corge']
    with Model(names) as m:
        parameterized(m,'qux',[0,0,0,.72])
        r=rate_default(0);r.p[0],r.p[2]=12.2,300
        reactions=[(['foo','baz'],{'bar':1,'quuz':2.4}),(['bar','qux'],{'foo':1,'quz':1.4}),
                   (['quz'],{}),(['baz','qux'],{'bar':1,'quz':2.5})]
        for reactants,products in reactions:m.reaction(reactants,products,r)
        m.finalize()
        assert [model_species_used(m.handle,m.ids[n]) for n in names+['qux']]==[1,1,1,1,1,0,1]
        pattern=np.zeros((6,6),np.uint8)
        require(model_kinetic_pattern(m.handle,pattern.ctypes.data_as(BP),36))
        expected_pairs=[(0,0),(0,1),(0,2),(1,0),(1,1),(1,2),(2,0),(2,2),(3,1),(3,2),(3,3),(4,0),(4,2)]
        assert list(zip(*np.nonzero(pattern)))==expected_pairs
        with State(m,2) as s:
            s.y[:]=np.array([[.1,.2,.3,.4,.5,0],[1.1,1.2,1.3,1.4,1.5,0]]).T
            s.rates[:]=[[10,110],[20*70*.72,120*80*.72],[30,130],[40*70*.72,140*80*.72]]
            expected_f=np.full((6,2),1000.)
            expected_j=np.full((2,6,6),100.)
            # Independent analytical mass-action derivatives, including duplicate reactants.
            for ri,(reactants,products) in enumerate(reactions):
                ids=[m.indices[x] for x in reactants if x!='qux']
                for c in range(2):
                    v=s.rates[ri,c]
                    for i in ids:v*=s.y[i,c]
                    for i in ids:expected_f[i,c]-=v
                    for p,yield_ in products.items():expected_f[m.indices[p],c]+=yield_*v
                    for j,col in enumerate(ids):
                        derivative=s.rates[ri,c]
                        for k,i in enumerate(ids):
                            if k!=j:derivative*=s.y[i,c]
                        for i in ids:expected_j[c,i,col]+=derivative
                        for p,yield_ in products.items():expected_j[c,m.indices[p],col]-=yield_*derivative
            forcing=np.full_like(s.y,1000.)
            require(add_force(s.handle,ptr(s.y),ptr(forcing)))
            near(forcing,expected_f,rtol=4*np.finfo(float).eps)
            values=np.full((model_jacobian_size(m.handle),2),100.)
            require(add_negative_jacobian(s.handle,ptr(s.y),ptr(values)))
            for row,col in expected_pairs:
                actual=values[model_jacobian_index(m.handle,row,col)]
                near(actual,expected_j[:,row,col],rtol=4*np.finfo(float).eps)


@case('process_set:RandomProcessSet.Matrix')
def random_processes(_):
    rng=np.random.default_rng(23)
    for cells,nreactions,nspecies in ((200,50,40),(300,30,20),(400,100,80)):
        with Model([str(i) for i in range(nspecies)]) as m:
            reactions=[]
            for ri in range(nreactions):
                reactants=list(rng.integers(nspecies,size=rng.integers(4)))
                products=list(rng.integers(nspecies,size=rng.integers(11)))
                m.reaction([str(i) for i in reactants],[(str(i),1.2) for i in products],rate_default(0))
                reactions.append((reactants,products))
            m.finalize()
            with State(m,cells) as s:
                s.y[:]=rng.lognormal(-2,4,s.y.shape);s.rates[:]=rng.lognormal(-2,4,s.rates.shape)
                expected=np.full_like(s.y,1000.)
                for ri,(reactants,products) in enumerate(reactions):
                    rate=s.rates[ri].copy()
                    for i in reactants:rate*=s.y[i]
                    for i in reactants:expected[i]-=rate
                    for i in products:expected[i]+=rate*1.2
                actual=np.full_like(s.y,1000.);require(add_force(s.handle,ptr(s.y),ptr(actual)))
                # Source checks successful evaluation only; add a numerical oracle.
                near(actual,expected,rtol=1e-11,atol=1e-8)


@case('process_set:ProcessSetAlgebraicVariables.CudaMatrix')
def algebraic_masking(_):
    with Model(['A','B','C','D']) as m:
        m.reaction(['A','B'],{'C':1,'D':1},rate_default(0))
        m.linear('B',{'B':1},0,'maskB');m.linear('D',{'D':1},0,'maskD');m.finalize()
        with State(m,2) as s:
            s.y[:]=np.array([[1,2,3,4],[5,6,7,8]]).T;s.rates[:]=[[10,20]]
            f=np.full_like(s.y,1000.);require(add_force(s.handle,ptr(s.y),ptr(f)))
            near(f.T,[[980,1000,1020,1000],[400,1000,1600,1000]],rtol=0)
            j=np.full((model_jacobian_size(m.handle),2),500.)
            require(add_negative_jacobian(s.handle,ptr(s.y),ptr(j)))
            expected=[[[520,510],[500,500],[480,490],[500,500]],[[620,600],[500,500],[380,400],[500,500]]]
            for row in range(4):
                for col in range(2):near(j[model_jacobian_index(m.handle,row,col)],np.array(expected)[:,row,col],rtol=0)


def verify_state_jacobian(s,atol=1e-7,rtol=1e-7):
    n,nc=s.y.shape
    @ForcingFunction
    def function(n_,c_,y,f,_):force(s.handle,y,f)
    fd=np.zeros((n,n,nc));require(finite_difference(n,nc,ptr(s.y),1e-8,function,None,ptr(fd)))
    negative=-np.moveaxis(s.jacobian(),0,-1).copy()
    pattern=np.zeros((n,n),np.uint8);require(model_kinetic_pattern(s.model.handle,pattern.ctypes.data_as(BP),n*n))
    comparison=compare_jacobian(n,nc,pattern.ctypes.data_as(BP),ptr(negative),ptr(fd),atol,rtol)
    assert comparison.passed,(comparison.max_absolute_error,comparison.worst_row,comparison.worst_column)
    assert check_sparsity(n,nc,pattern.ctypes.data_as(BP),ptr(fd),1e-6).passed


@case('process_set:ProcessSetFiniteDifferenceJacobian.Matrix')
def finite_jacobian(_):
    with Model(['A','B','C']) as m:
        m.reaction(['A','B'],{'C':2},rate_default(0));m.reaction(['C'],{'A':1},rate_default(0));m.finalize()
        with State(m,2) as s:
            s.y[:]=np.array([[1,2,3],[.5,4,1.5]]).T;s.rates[:]=[[5,8],[10,12]]
            verify_state_jacobian(s)
