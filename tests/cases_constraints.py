"""Built-in algebraic constraints; all residual/Jacobian checks are numerical."""
from ffi import *
from harness import group


def finite_difference(state, rtol=1e-6, atol=1e-8):
    analytic = state.jacobian(constraints=True)
    original = state.y.copy()
    nc,n = state.cells,len(state.model.names)
    numerical = np.zeros_like(analytic)
    for i in range(n):
        step = 1e-5*np.maximum(np.abs(original[i]),1e-3)
        state.y[i] = original[i]+step
        plus = state.forcing(constraints=True)
        state.y[i] = original[i]-step
        minus = state.forcing(constraints=True)
        state.y[:] = original
        numerical[:,:,i] = ((plus-minus)/(2*step)).T
    near(analytic,numerical,rtol=rtol,atol=atol)


LINEAR_NAMES = '''Construction AlgebraicSpecies WeightedTerms ZeroConstant FractionalCoefficients
ResidualComputationThroughConstraintSet JacobianComputationThroughConstraintSet WeightedSumResidualAndJacobian
ThreeSpeciesConservationResidual ZeroConstantResidual FractionalCoefficientsResidualAndJacobian
JacobianIndependentOfConcentrations FiniteDifferenceJacobianSimpleConservation FiniteDifferenceJacobianWeightedSum'''


@group('linear:LinearConstraint.', LINEAR_NAMES)
def linear(id_):
    name = id_.split('.')[-1]
    if name in ('AlgebraicSpecies','ThreeSpeciesConservationResidual'):
        terms,constant,states = {'A':1.,'B':1.,'C':1.},10.,[[2,3,5],[2,3,6]]
    elif 'Weighted' in name:
        terms,constant,states = {'A':2.,'B':3.,'C':-1.},5.,[[1,2,3],[.2,.5,.8]]
    elif 'Fractional' in name:
        terms,constant,states = {'A':.5,'B':1.5},2.,[[1,1],[2,.5]]
    elif 'ZeroConstant' in name:
        terms,constant,states = {'A':1.,'B':-1.},0.,[[.5,.5],[.6,.5]]
    else:
        terms,constant,states = {'A':1.,'B':1.},1.,[[.3,.7],[.5,.6]]
    with Model(terms) as m:
        row = list(terms)[-1]
        m.linear(row,terms,constant)
        m.finalize()
        assert model_constraint_count(m.handle) == 1
        assert model_is_algebraic(m.handle,m.indices[row]) == 1
        with State(m,len(states)) as s:
            s.y[:] = np.array(states).T
            expected = np.zeros_like(s.y)
            expected[-1] = np.array(states)@np.array(list(terms.values()))-constant
            near(s.forcing(True),expected,rtol=0,atol=1e-10)
            jac = np.zeros((len(states),len(terms),len(terms)))
            jac[:,-1,:] = list(terms.values())
            near(s.jacobian(True),jac,rtol=0,atol=1e-10)
            if name == 'JacobianIndependentOfConcentrations':
                s.y[:] *= 10
                near(s.jacobian(True),jac,rtol=0,atol=1e-10)
            if name.startswith('FiniteDifference'):
                finite_difference(s)


EQUILIBRIUM_NAMES = '''Construction AlgebraicSpecies SingleReactantSingleProduct MultipleReactantsAndProducts
InvalidEquilibriumConstant EmptyReactantsThrows EmptyProductsThrows InvalidStoichiometryThrows
ResidualComputationThroughConstraintSet JacobianComputationThroughConstraintSet ComplexStoichiometryResidual
FiniteDifferenceJacobianSimple FiniteDifferenceJacobianComplexStoichiometry'''


@group('equilibrium:EquilibriumConstraint.', EQUILIBRIUM_NAMES)
def equilibrium(id_):
    name = id_.split('.')[-1]
    if name in ('InvalidEquilibriumConstant','EmptyReactantsThrows','EmptyProductsThrows','InvalidStoichiometryThrows'):
        with Model(['A','B']) as m:
            configurations = [(1,[0],[1.],1,[1],[1.],1.)]
            if name == 'InvalidEquilibriumConstant':
                configurations = [(1,[0],[1.],1,[1],[1.],k) for k in (-1,0)]
            elif name == 'EmptyReactantsThrows': configurations = [(0,[],[],1,[1],[1.],1.)]
            elif name == 'EmptyProductsThrows': configurations = [(1,[0],[1.],0,[],[],1.)]
            else:
                configurations = [(1,[0],[k],1,[1],[1.],1.) for k in (0,-1)] + [(1,[0],[1.],1,[1],[k],1.) for k in (0,-2)]
            for nr,r,cr,np_,p,cp,k in configurations:
                assert model_add_equilibrium(m.handle,b'invalid',1,nr,ints(r),doubles(cr),np_,ints(p),doubles(cp),k,-2400.) == -1
        return
    complex_ = name in ('MultipleReactantsAndProducts','ComplexStoichiometryResidual','FiniteDifferenceJacobianComplexStoichiometry')
    simple = name == 'SingleReactantSingleProduct'
    names = ['A','B'] if simple else (['A','B','C'] if complex_ else ['A','B','AB'])
    with Model(names) as m:
        if complex_:
            m.equilibrium('B',{'A':2.},{'B':1.,'C':1.},100.)
            rows = [[.1,.5,2.],[.15,.8,1.5]]
        elif simple:
            m.equilibrium('B',{'A':1.},{'B':1.},10.)
            rows = [[.1,1.],[.2,1.]]
        else:
            m.equilibrium('AB',{'A':1.,'B':1.},{'AB':1.},1000.)
            rows = [[.001,.001,.001],[.02,.03,.05],[.01,.02,.05],[.1,.2,5.]]
        m.finalize()
        assert model_constraint_count(m.handle) == 1
        assert model_parameter_name(m.handle,0) == b'equilibrium'
        with State(m,len(rows)) as s:
            s.y[:] = np.array(rows).T
            f,j = np.zeros_like(s.y),np.zeros((len(rows),len(names),len(names)))
            for c,y in enumerate(rows):
                if complex_:
                    f[1,c] = 100*y[0]**2-y[1]*y[2]
                    j[c,1,:] = [200*y[0],-y[2],-y[1]]
                elif simple:
                    f[1,c] = 10*y[0]-y[1]; j[c,1,:] = [10,-1]
                else:
                    f[2,c] = 1000*y[0]*y[1]-y[2]; j[c,2,:] = [1000*y[1],1000*y[0],-1]
            near(s.forcing(True),f,rtol=0,atol=1e-10)
            near(s.jacobian(True),j,rtol=0,atol=1e-10)
            if name.startswith('FiniteDifference'):
                finite_difference(s)


SET_NAMES = '''Construction ReplaceStateRowsMapsToAlgebraicSpecies NonZeroJacobianElements MultipleConstraints
AddForcingTerms SubtractJacobianTerms EmptyConstraintSet UnknownSpeciesThrows ThreeDStateOneConstraint
FourDStateTwoConstraints CoupledConstraintsSharedSpecies VectorizedMatricesRespectGridCellIndexing
VectorMatrix1 VectorMatrix2 VectorMatrix3'''


@group('constraint_set:ConstraintSet.',SET_NAMES)
def constraint_set(id_):
    name = id_.split('.')[-1]
    if name == 'UnknownSpeciesThrows':
        with Model(['A','B']) as m:
            assert model_add_equilibrium(m.handle,b'unknown',2,2,ints([3,4]),doubles([1,1]),1,ints([2]),doubles([1]),.033,-24000) < 0
        return
    if name == 'EmptyConstraintSet':
        with Model(['A','B','C']) as m:
            m.finalize()
            assert model_constraint_count(m.handle) == 0
            with State(m) as s:
                s.y[:] = [[.1],[.2],[.3]]
                values = np.ones((4,1))
                constraint_force(s.handle,ptr(s.y),ptr(values))
                near(values,np.ones((4,1)),rtol=0)
        return
    specs,rows = [],[]
    if name == 'ReplaceStateRowsMapsToAlgebraicSpecies':
        names = ['A','B','C']; specs = [('C',{'B':1},{'C':1})]; rows = [[1,2,3]]
    elif name == 'ThreeDStateOneConstraint':
        names = ['X','Y','Z']; specs = [('Y',{'X':1},{'Y':1})]; rows = [[10,.2,.5],[10,.33,.3]]
    elif name == 'FourDStateTwoConstraints':
        names = ['A','B','C','D']; specs = [('B',{'A':1},{'B':1}),('A',{'C':1,'D':1},{'A':1})]
        rows = [[.33,.01089,10,1],[.33,.01089,5,1],[.33,.02,5,1]]
    elif name == 'CoupledConstraintsSharedSpecies':
        names = ['A','B','C']; specs = [('B',{'A':1},{'B':1}),('C',{'A':1},{'C':1})]
        rows = [[10,.33,.33],[10,.2,.4]]
    elif name == 'MultipleConstraints':
        names = ['A','B','AB','C','D']; specs = [('AB',{'A':1,'B':1},{'AB':1}),('D',{'C':1},{'D':1})]
        rows = [[.01,.02,.05,.1,.2]]
    else:
        names = ['A','B','AB']; specs = [('AB',{'A':1,'B':1},{'AB':1})]
        rows = [[.2,.4,.6],[.1,.3,.7],[.01,.02,.05]]
        if name == 'VectorizedMatricesRespectGridCellIndexing':
            rows = [[.1*(i+1),.2*(i+2),.3*(i+3)] for i in range(7)]
    with Model(names) as m:
        for i,(row,r,p) in enumerate(specs):
            m.equilibrium(row,r,p,.033,name='eq'+str(i))
        m.finalize()
        assert model_constraint_count(m.handle) == len(specs)
        for row,r,p in specs:
            ri = m.indices[row]
            assert model_is_algebraic(m.handle,ri) == 1
            dependencies = {m.indices[x] for x in (r|p)}
            present = {col for col in range(len(names)) if model_jacobian_index(m.handle,ri,col) >= 0}
            assert dependencies <= present
        with State(m,len(rows)) as s:
            s.y[:] = np.array(rows).T
            expected = np.zeros_like(s.y)
            jac = np.zeros((len(rows),len(names),len(names)))
            for row,r,p in specs:
                ri = m.indices[row]
                for cell in range(s.cells):
                    def product(terms):
                        v = 1.
                        for name,power in terms.items(): v *= s.y[m.indices[name],cell]**power
                        return v
                    expected[ri,cell] = .033*product(r)-product(p)
                    for terms,scale in [(r,.033),(p,-1.)]:
                        for name,power in terms.items():
                            jac[cell,ri,m.indices[name]] += scale*product(terms)*power/s.y[m.indices[name],cell]
            near(s.forcing(True),expected,rtol=0,atol=1e-10)
            near(s.jacobian(True),jac,rtol=0,atol=1e-10)
