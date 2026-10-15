"""External constraint composition, filtering, parameters, and FD case ports."""
from functools import partial
import math
from builderffi import *
from test_dae_integration import eq,arr
from test_jacobian import finite_difference,Comparison,Sparse

CONSTRAINT=C.CFUNCTYPE(D,D,DP,DP,DP,P)


def callback_constraint(owner,data,row,terms,constant=0.,parameter=None):
    indices=[(lookup(data,n),v) for n,v in terms]
    c=lib.chem_constraint_default(2);c.algebraic_species=lookup(data,row)
    c.reactant_count=len(indices);c.reactants=owner.array(Product,[Product(i,1.) for i,v in indices])
    def evaluate(t,y,custom,gradient,ctx):
        value=-constant
        for position,(i,coefficient) in enumerate(indices):
            a=custom[parameter] if parameter is not None and position==0 else coefficient
            value+=a*y[i]
            if gradient:gradient[i]+=a
        return value
    fn=CONSTRAINT(evaluate);owner.keep.append(fn);c.function=C.cast(fn,P)
    return c


class Equilibria(Definition):
    def __init__(self,ratios=(),mass=None,temperature=False):
        self.ratios,self.mass,self.temperature=ratios,mass,temperature
        super().__init__(parameters=['C_K_eq'] if temperature else [],constraints=len(ratios)+(mass is not None))
    def configure(self,data,process,cs):
        parameter=lookup(data,'C_K_eq',True) if self.temperature else None
        for j,(r,p,k) in enumerate(self.ratios):cs[j]=callback_constraint(self,data,p,[(r,k),(p,-1.)],parameter=parameter)
        if self.mass:
            row,names,total=self.mass
            cs[len(self.ratios)]=callback_constraint(self,data,row,[(n,1.) for n in names],total)
        def update(state,ctx):
            if parameter is not None:
                nc=lib.chem_cell_count(state);custom=lib.chem_parameters(state);temp=lib.chem_temperatures(state)
                for c in range(nc):custom[parameter*nc+c]=2.*math.exp(3000.*(1/298.15-1/temp[c]))
        self.attach(process,update=update)


class Aerosol(Definition):
    def __init__(self,k,total=None,solvent=False):
        self.k,self.total,self.solvent=k,total,solvent
        super().__init__(variables=['AEROSOL.A_AQ'],constraints=int(total is not None),used=['A_GAS'])
    def configure(self,data,process,cs):
        gas,aq=lookup(data,'A_GAS'),lookup(data,'AEROSOL.A_AQ')
        solvent=lookup(data,'S') if self.solvent else None
        coordinates=[(gas,gas),(aq,gas)]+([(gas,solvent),(aq,solvent)] if self.solvent else [])
        coords=self.array(MatrixIndex,[MatrixIndex(*ij) for ij in coordinates])
        process.nonzero_count=len(coords);process.nonzeros=coords
        flats={}
        def finalize(m,ctx):
            for r,c in coordinates:
                flats[r,c]=lib.chem_jacobian_index(m,r,c);assert flats[r,c]!=Z(-1).value
        def forcing(state,y,f,ctx):
            nc=lib.chem_cell_count(state)
            for c in range(nc):
                value=self.k*y[gas*nc+c]
                if solvent is not None:value*=y[solvent*nc+c]
                f[gas*nc+c]-=value;f[aq*nc+c]+=value
        def jacobian(state,y,jac,ctx):
            nc=lib.chem_cell_count(state)
            for c in range(nc):
                derivative=self.k*(y[solvent*nc+c] if solvent is not None else 1.)
                jac[flats[gas,gas]*nc+c]+=derivative;jac[flats[aq,gas]*nc+c]-=derivative
                if solvent is not None:
                    jac[flats[gas,solvent]*nc+c]+=self.k*y[gas*nc+c]
                    jac[flats[aq,solvent]*nc+c]-=self.k*y[gas*nc+c]
        self.attach(process,forcing,jacobian,finalize=finalize)
        if self.total is not None:cs[0]=callback_constraint(self,data,'AEROSOL.A_AQ',[('A_GAS',1.),('AEROSOL.A_AQ',1.)],self.total)


class NoOp(Definition):
    def configure(self,data,process,cs):self.attach(process,lambda *args:None,lambda *args:None)


def setup(m,rows,temp=298.):
    m.set_y(rows)
    for c in range(m.nc):lib.chem_temperatures(m.state)[c]=temp;lib.chem_pressures(m.state)[c]=101325.


def step(m,b,dt=1.,method=-4):
    assert lib.chem_update_rates(m.state)==0
    status,stats=m.solve(dt,method)
    b.check();assert status==0,(status,stats.final_time)
    return m.rows()


def metadata(total=1.):
    model=Aerosol(.01,total)
    b=Builder(['A_GAS'],models=[model])
    with b.build() as m:
        assert m.n==2 and m.data.constraint_count==int(total is not None)
        assert [lib.chem_is_algebraic(m.handle,j) for j in range(2)]==[0,int(total is not None)]
        b.check()


def conservation(k=.1,initial=.9,steps=20,dt=10.,cells=1,solvent=False):
    b=Builder(['A_GAS','S'] if solvent else ['A_GAS'],models=[Aerosol(k,1.,solvent)])
    with b.build(cells) as m:
        rows=[]
        for c in range(cells):
            gas=initial-.2*c
            rows.append([gas,1.,1-gas] if solvent else [gas,1-gas])
        setup(m,rows)
        aq=lib.chem_species_index(m.handle,b'AEROSOL.A_AQ')
        assert lib.chem_is_algebraic(m.handle,aq)==1
        if solvent:assert lib.chem_jacobian_index(m.handle,aq,1)!=Z(-1).value
        for _ in range(steps):
            for row in step(m,b,dt):assert abs(row[0]+row[aq]-1.)<=1e-4


def algebraic_terms():
    b=Builder(['A_GAS'],models=[Aerosol(10.,1.)])
    with b.build() as m:
        setup(m,[[.7,.2]])
        f=(D*2)();lib.chem_forcing(m.state,m.y,f);lib.chem_constraint_forcing(m.state,m.y,f)
        assert abs(f[0]+7.)<=1e-12 and abs(f[1]-(.7+.2-1.))<=1e-12
        jac=(D*lib.chem_jacobian_size(m.handle))()
        lib.chem_negative_jacobian(m.state,m.y,jac);lib.chem_constraint_negative_jacobian(m.state,m.y,jac)
        for ij,value in [((0,0),10.),((1,0),-1.),((1,1),-1.)]:assert abs(jac[lib.chem_jacobian_index(m.handle,*ij)]-value)<=1e-12
        b.check()


def combined(builtin_row=False):
    if builtin_row:
        constraints=[eq(2,[(1,1.)],[(2,1.)],5.)]
        b=Builder(['A_GAS','B'],constraints=constraints,models=[Aerosol(.1)])
    else:
        constraints=[eq(2,[(1,1.)],[(2,1.)],5.)]
        b=Builder(['A_GAS','B','C'],[([0],[(1,1.)],arr(.05))],constraints,models=[Aerosol(.01,1.)])
    with b.build() as m:
        aq=lib.chem_species_index(m.handle,b'AEROSOL.A_AQ')
        assert m.data.constraint_count==(1 if builtin_row else 2)
        assert lib.chem_is_algebraic(m.handle,aq)==1
        if builtin_row:assert lib.chem_jacobian_index(m.handle,aq,0)!=Z(-1).value
        else:assert lib.chem_is_algebraic(m.handle,2)==1
        b.check()


def trajectory(kind,steps=200,dt=1.):
    reactions=[([0],[(1,1.)],arr(.1))]
    models=[]
    if kind=='kinetic':reactions += [([1],[(2,1.)],arr(10.)),([2],[(1,1.)],arr(2.))]
    elif kind=='conservative':models=[Equilibria([('B','C',5.)],('B',list('ABC'),1.))]
    else:models=[Equilibria([('B','C',5.)])]
    b=Builder(list('ABC'),reactions,models=models)
    with b.build() as m:
        setup(m,[[1.,0.,0.]])
        for _ in range(steps):row=step(m,b,dt,3 if kind=='kinetic' else -4)[0]
        return row


def converged(ratio_only=False):
    kin=trajectory('kinetic');other=trajectory('simple' if ratio_only else 'conservative')
    if ratio_only:
        for row in [kin,other]:assert abs(row[2]/row[1]-5.)<=1e-3
        assert abs(sum(kin)-1.)<=1e-6 and abs(sum(other)-6.)<=.1
    else:
        for row in [kin,other]:
            for a,b,tol in zip(row,[0.,1/6,5/6],[1e-6,1e-3,1e-3]):assert abs(a-b)<=tol
        for a,b,tol in zip(kin,other,[1e-6,1e-3,1e-3]):assert abs(a-b)<=tol
        assert abs(sum(kin)-1.)<=1e-6 and abs(sum(other)-1.)<=1e-3


def equivalent():
    rx=[([0],[(1,1.)],arr(.1))]
    bi=Builder(list('ABC'),rx,[eq(2,[(1,1.)],[(2,1.)],5.)]);ex=Builder(list('ABC'),rx,models=[Equilibria([('B','C',5.)])])
    with bi.build() as a,ex.build() as b:
        setup(a,[[1.,0.,0.]]);setup(b,[[1.,0.,0.]])
        for _ in range(100):
            aa=step(a,bi,.5)[0];bb=step(b,ex,.5)[0]
            for x,y in zip(aa,bb):assert abs(x-y)<=1e-8
            for row in [aa,bb]:assert abs(5*row[1]-row[2])<=1e-6


def composed():
    rx=[([0],[(1,1.)],arr(.1))]
    kinetic=rx+[([1],[(2,1.)],arr(10.)),([2],[(1,1.)],arr(2.)),([1],[(3,1.)],arr(8.)),([3],[(1,1.)],arr(8/3))]
    models=[Equilibria([('B','C',5.)]),Equilibria([('B','D',3.)]),Equilibria(mass=('B',list('ABCD'),1.))]
    ba=Builder(list('ABCD'),kinetic);bb=Builder(list('ABCD'),rx,models=models)
    with ba.build() as a,bb.build() as b:
        setup(a,[[1.,0.,0.,0.]]);setup(b,[[1.,0.,0.,0.]])
        for _ in range(200):aa=step(a,ba,method=3)[0];cc=step(b,bb)[0]
        for row in [aa,cc]:
            for x,y,tol in zip(row,[0.,1/9,5/9,3/9],[1e-6,1e-3,1e-3,1e-3]):assert abs(x-y)<=tol
        for x,y in zip(aa[1:],cc[1:]):assert abs(x-y)<=1e-3
        assert abs(sum(aa)-1.)<=1e-6 and abs(sum(cc)-1.)<=1e-3
        assert abs(5*cc[1]-cc[2])<=1e-6 and abs(3*cc[1]-cc[3])<=1e-6


def fd(mode='process'):
    if mode=='equilibrium':b=Builder(list('AB'),models=[Equilibria([('A','B',2.5)])]);rows=[[3.,5.]]
    else:b=Builder(['A_GAS'],models=[Aerosol(3.5 if mode=='process' else .1,None if mode=='process' else 1.)]);rows=[[.8,.2],[.3,.7]] if mode=='process' else [[.6,.4],[.2,.8]]
    with b.build(len(rows)) as m:
        setup(m,rows)
        function=lib.chem_forcing if mode=='process' else lib.chem_constraint_forcing
        derivative=lib.chem_negative_jacobian if mode=='process' else lib.chem_constraint_negative_jacobian
        finite=finite_difference(rows,lambda nc,n,y,f:function(m.state,y,f))
        values=(D*(lib.chem_jacobian_size(m.handle)*m.nc))();derivative(m.state,m.y,values)
        if mode=='process':pattern=[(0,0),(1,0)]
        else:pattern=[(1,0),(1,1)]
        with Sparse(2,blocks=m.nc,entries=pattern) as a:
            for c in range(m.nc):
                for i,j in pattern:a[c,i,j]=values[lib.chem_jacobian_index(m.handle,i,j)*m.nc+c]
            comparison=Comparison()
            assert lib.chem_compare_jacobian(a.handle,finite,1e-7,1e-7,C.byref(comparison))==0 and comparison.passed
            if mode!='constraint':assert lib.chem_check_jacobian_sparsity(a.handle,finite,1e-6,C.byref(comparison))==0 and comparison.passed
        b.check()


def temperature():
    b=Builder(list('ABC'),[([0],[(1,1.)],arr(.1))],models=[Equilibria([('B','C',2.)],temperature=True)])
    for temp in [298.15,350.]:
        with b.build() as m:
            setup(m,[[1.,0.,0.]],temp)
            for _ in range(200):row=step(m,b)[0]
            k=2.*math.exp(3000.*(1/298.15-1/temp))
            assert row[1]>0 and abs(row[2]/row[1]-k)<=1e-4
            if temp==350.:assert k>2.


def order():
    results=[]
    for reverse in [False,True]:
        models=[Equilibria([('B','C',5.)]),NoOp()]
        if reverse:models.reverse()
        b=Builder(list('ABC'),[([0],[(1,1.)],arr(.1))],models=models)
        with b.build() as m:
            assert m.data.constraint_count==1 and lib.chem_is_algebraic(m.handle,2)==1
            setup(m,[[1.,0.,0.]])
            for _ in range(20):row=step(m,b)[0]
            assert row[1]>0 and abs(5*row[1]-row[2])<=1e-6
            results.append(row[1:])
    assert results[0]==results[1]


def reuse():
    b=Builder(list('ABCD'),[([0],[(1,1.)],arr(.1))],[eq(3,[(0,1.)],[(3,1.)],2.)],models=[Equilibria([('B','C',5.)])])
    with b.build() as first,b.build() as second:
        results=[]
        for m in [first,second]:
            assert m.data.constraint_count==2
            setup(m,[[1.,0.,0.,2.]])
            for _ in range(20):row=step(m,b)[0]
            results.append(row)
        a,bvalue,c,d=results[0]
        assert abs(2*a-d)<=1e-6 and abs(5*bvalue-c)<=1e-6
        assert results[0]==results[1]


CASES={
'ExternalModelConstraints.AddExternalModelWithConstraints':metadata,
'ExternalModelConstraints.AddExternalModelProcessOnly':partial(metadata,None),
'ExternalModelConstraints.DAESolveEnforcesConservation':conservation,
'ExternalModelConstraints.AlgebraicRowHoldsOnlyConstraintTerms':algebraic_terms,
'ExternalModelConstraints.CombinedBuiltInAndExternalConstraints':combined,
'ExternalModelConstraints.AddExternalModelOnlyStandardRosenbrock':partial(metadata,None),
'ExternalModelConstraints.AddExternalModelConstraintsOnly':partial(conservation,.01,.8,1,50.),
'ExternalModelConstraints.MultiGridCell':partial(conservation,cells=3,steps=1,dt=50.),
'ExternalModelConstraints.KineticVsConservativeConstraintConvergence':converged,
'ExternalModelConstraints.SimpleConstraintPreservesRatioNotMass':partial(converged,True),
'ExternalModelConstraints.BuiltInVsExternalModelConstraintStepByStep':equivalent,
'ExternalModelConstraints.MultiEquilibriumKineticVsComposedConstraints':composed,
'ExternalModelConstraints.ProcessJacobianElementInAlgebraicRowSurvivesFiltering':partial(conservation,solvent=True,steps=10),
'ExternalModelConstraints.ProcessJacobianElementInBuiltInAlgebraicRowSurvivesFiltering':partial(combined,True),
'ExternalModelConstraints.TemperatureDependentConstraintParameter':temperature,
'ExternalModelConstraints.ConstraintActiveWhenProcessOnlyModelAddedFirst':order,
'ExternalModelConstraints.BuilderCanBuildMultipleSolvers':reuse,
'ExternalModelFiniteDifferenceJacobian.ProcessForcingJacobian':fd,
'ExternalModelFiniteDifferenceJacobian.ConstraintResidualJacobian':partial(fd,'constraint'),
'ExternalModelFiniteDifferenceJacobian.EquilibriumConstraintModelJacobian':partial(fd,'equilibrium'),
}
